











##########################################################################
set.seed(1234)

library(ggplot2)
library(ggstar)
library(rELA)
library(vegan)
library(foreach)
library(doParallel)
library(tidyr)
library(compositions)


source("packages/01_1_function.R")
#########################################################################

#setwd("/Users/ngch/Desktop/Now_using/data/sugaFB_ELA")

save.dir <- "02_01_ELA_prep_abundance_threshold"
dir.create(save.dir)

#####################
#parameters
tg_th <- 100
#qth <- 10^-5 # do not change!!



taxa <- "Family"

##select target taxa
unident_th <- 0.5

########################################################################
tx_f <- readRDS("Base_data/Fungi/taxa_list_mod.rds")
tx_b <- as.data.frame(readRDS("Base_data/OTU_basedata_set/taxonomy_list_dada2.rds"))


#Genus abundance in all dataset

inputdir <- "From_localcomputer/01_caverage_rarefaction"
tb_al <- list(Fungi=readRDS(sprintf("%s/covrarefy_sqtb_fungi.rds",inputdir))$table,
             Prokaryota=readRDS(sprintf("%s/covrarefy_sqtb_bacteria.rds",inputdir))$table)


if(taxa=="OTU"){
  tb_m1 <- tb_al
}else{
  #Genus
  tb_m0 <-list(Fungi=Taxa.mat(tb_al$Fungi, tx_f, taxa),
               Prokaryota=Taxa.mat(tb_al$Prokaryota, tx_b, taxa))
  tm <- list(Fungi=tb_m0$Fungi/rowSums(tb_m0$Fungi),
             Prokaryota=tb_m0$Prokaryota/rowSums(tb_m0$Prokaryota))
  tb_m1 <- list(Fungi=tm$Fungi[which(tm$Fungi[,"Unidentified"]<unident_th),],
                Prokaryota=tm$Prokaryota[which(tm$Prokaryota[,"Unidentified"]<unident_th),])
}


inf1 <- readRDS("Base_data/comp_sample_info.rds")

inff <- inf1[which(inf1$ID %in% rownames(tb_m1$Fungi)),]
plnamf <- names(table(inff[inff$target=="fungi","plant"])[table(inff[inff$target=="fungi","plant"])>50])

infb <- inf1[which(inf1$ID %in% rownames(tb_m1$Prokaryota)),]
plnamb <- names(table(infb[infb$target=="fungi","plant"])[table(infb[infb$target=="fungi","plant"])>50])

plnam <- intersect(plnamf,plnamb)
plnm <- setdiff(plnam,c("-","Unidentified"))
info <- inf1[which(inf1$plant%in%plnm),]


saveRDS(info,sprintf("%s/comp_sample_info_plant2.rds",save.dir))


sp_info <- spread(data.frame(id=info$ID,plant=info$plant,val=1),key=plant,value=val)
rownames(sp_info) <- sp_info$id
sp_info[is.na(sp_info)] <- 0


saveRDS(sp_info,sprintf("%s/ELA_input_plant.rds",save.dir))


################
tb_all <- list(Fungi=tb_al$Fungi[which(rownames(tb_al$Fungi) %in% rownames(sp_info)),],
               Prokaryota=tb_al$Prokaryota[which(rownames(tb_al$Prokaryota) %in% rownames(sp_info)),])
tb_m <- list(Fungi=tb_m1$Fungi[which(rownames(tb_m1$Fungi) %in% rownames(sp_info)),],
             Prokaryota=tb_m1$Prokaryota[which(rownames(tb_m1$Prokaryota) %in% rownames(sp_info)),])
saveRDS(tb_m,sprintf("%s/Df_RA_%s.rds",save.dir,taxa))
##########################################

##select target taxa

tbgf <-Taxa.mat(tb_all$Fungi, tx_f, "Genus")
tbgb <- Taxa.mat(tb_all$Prokaryota, tx_b, "Genus")

clrf <- clr((tbgf+1)/rowSums(tbgf+1))
clrb <- clr((tbgb+1)/rowSums(tbgb+1))

tgf <- clrf[,which(colSums(tbgf>0)>tg_th)]
tgb <- clrb[,which(colSums(tbgb>0)>tg_th)]


saveRDS(list(Fungi=tbgf/rowSums(tbgf),Prokaryota=tbgb/rowSums(tbgb)),
        sprintf("%s/NoCLR_ExpVar_clrRA_target_taxa.rds",save.dir))

saveRDS(list(Fungi=tgf,Prokaryota=tgb),
        sprintf("%s/ExpVar_clrRA_target_taxa.rds",save.dir))

###########################################

fb <- "Fungi"

relf <- tb_m$Fungi/rowSums(tb_m$Fungi)

relf2 <- relf[,which(colnames(relf) != "Unidentified" & colSums(relf>0)/nrow(relf) > 0.03)]

grelf <- gather(cbind(ID=rownames(relf2),as.data.frame(relf2)),
       key=taxa,value="abundance",-1) 

grelf2 <- grelf[grelf$abundance!=0,]

grelf2$logRA <- log(grelf2$abundance)

grelf2$noise <- "Pres"
th_df <- NULL
for(i in 1:length(unique(grelf2$taxa))){#i <- 1
  show.progress(i,1:length(unique(grelf2$taxa)))
  grf <- grelf2[which(grelf2$taxa==unique(grelf2$taxa)[i]),] 
  ra <- grf$logRA
  
  grf[which(grf$logRA < mean(ra)-2*sd(ra)),"noise"] <- "low (< M-2SD)"
  grelf2[which(grelf2$taxa==unique(grelf2$taxa)[i]),"noise"] <- grf$noise
  
  th_df <- rbind(th_df,
                 data.frame(unique(grelf2$taxa)[i],mean(ra)-2*sd(ra)))
}

dimnames(th_df) <- list(th_df[,1],c("taxa","th_2SD"))
saveRDS(th_df,sprintf("%s/threshold_%s_%s.rds",save.dir,fb,taxa))

 grelf2$taxa2 <- factor(grelf2$taxa,levels=colnames(relf2)[order(colSums(relf2))])
 saveRDS(grelf2,sprintf("%s/graphics_basefile_abundance_threshold_%s_%s.rds",save.dir,fb,taxa))
 
 g <- ggplot(grelf2,aes(x=taxa2,y=logRA))+
   geom_point(aes(color=noise))+
   geom_boxplot(alpha=0.5,outlier.color = NA)+
   theme(axis.text.x = element_text(angle = 90, hjust = 1))+
   scale_color_manual(values=c("low (< M-2SD)"="gray","Pres"="red"))+
   labs(x=taxa,y="log(Relative abundance)",fill="Noise level")+
   ggtitle("Fungi Family")+
   theme(plot.title = element_text(hjust = 0.5))+
   theme_bw()+
   coord_flip()
 g
 
ggsave(sprintf("%s/outlier_Family_%s.pdf",save.dir,fb),g,width=6,height=9)

ocmat <- relf2
for(i in 1:ncol(relf2)){
  ocmat[,i] <- ifelse(ocmat[,i]<exp(th_df[colnames(ocmat)[i],2]),0,1)
}

saveRDS(ocmat,sprintf("%s/ocmat_remove_M2SD_%s_%s.rds",save.dir,fb,taxa))
saveRDS(relf2,sprintf("%s/ramat_for_%s_%s.rds",save.dir,fb,taxa))
###########################

fb <- "Prokaryota"

relf <- tb_m[[fb]]/rowSums(tb_m[[fb]])

relf2 <- relf[,which(colnames(relf) != "Unidentified" & colSums(relf>0) > 30)]

grelf <- gather(cbind(ID=rownames(relf2),as.data.frame(relf2)),
                key=taxa,value="abundance",-1) 

grelf2 <- grelf[grelf$abundance!=0,]

grelf2$logRA <- log(grelf2$abundance)

grelf2$noise <- "Pres"
th_df <- NULL
for(i in 1:length(unique(grelf2$taxa))){#i <- 1
  show.progress(i,1:length(unique(grelf2$taxa)))
  grf <- grelf2[which(grelf2$taxa==unique(grelf2$taxa)[i]),] 
  ra <- grf$logRA
  
  grf[which(grf$logRA < mean(ra)-2*sd(ra)),"noise"] <- "low (< M-2SD)"
  grelf2[which(grelf2$taxa==unique(grelf2$taxa)[i]),"noise"] <- grf$noise
  
  th_df <- rbind(th_df,
                 data.frame(unique(grelf2$taxa)[i],mean(ra)-2*sd(ra)))
}

dimnames(th_df) <- list(th_df[,1],c("taxa","th_2SD"))
saveRDS(th_df,sprintf("%s/threshold_%s_%s.rds",save.dir,fb,taxa))

grelf2$taxa2 <- factor(grelf2$taxa,levels=colnames(relf2)[order(colSums(relf2))])
saveRDS(grelf2,sprintf("%s/graphics_basefile_abundance_threshold_%s_%s.rds",save.dir,fb,taxa))

g <- ggplot(grelf2,aes(x=taxa2,y=logRA))+
  geom_point(aes(color=noise))+
  geom_boxplot(alpha=0.5,outlier.color = NA)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_color_manual(values=c("low (< M-2SD)"="gray","Pres"="red"))+
  labs(x=taxa,y="log(Relative abundance)",fill="Noise level")+
  ggtitle("Fungi Family")+
  theme(plot.title = element_text(hjust = 0.5))+
  theme_bw()+
  coord_flip()
g

ggsave(sprintf("%s/outlier_Family_%s.pdf",save.dir,fb),g,width=6,height=19)

ocmat <- relf2
for(i in 1:ncol(relf2)){
  ocmat[,i] <- ifelse(ocmat[,i]<exp(th_df[colnames(ocmat)[i],2]),0,1)
}

saveRDS(ocmat,sprintf("%s/ocmat_remove_M2SD_%s_%s.rds",save.dir,fb,taxa))
saveRDS(relf2,sprintf("%s/ramat_for_%s_%s.rds",save.dir,fb,taxa))