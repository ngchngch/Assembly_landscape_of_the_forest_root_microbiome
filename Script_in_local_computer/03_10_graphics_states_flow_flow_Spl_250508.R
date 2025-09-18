
##########################################################################
#function
set.seed(1234)

####################

#install.packages("mgcv")
library(ggplot2)
library(ggtext)
library(parallel)
library('RColorBrewer')
#library('scatterpie')
library('patchwork')
library("rELA")
library('tidyverse')
library("doParallel")
library(foreach)
#library(renv)
library("Rcpp")
library("RcppArmadillo")
library(vegan)
library(ggalluvial)
library(dplyr)
library(ComplexHeatmap)

library(circlize)

#######


dir_02_05 <- "Output_supercomputer/02_05_summary_taxa_select"

dir_03_01 <- 'Output_supercomputer/03_01_ELA_withRA_4step_250227'


dir_03_05 <- "Output_supercomputer/03_05_states_flow_diagram_select_250313"

dir_03_07 <- "Output/03_07_graphics_Zconv_landchanges_each_biplot_250312"
dir_02_06 <- "Output_supercomputer/02_06_ELA"

dir_02_07 <- "Output_supercomputer/02_07_assemblygraph_250311"

dir_02_06_02 <- "Output/02_06_02_Order_SSheatmap_250918"
#########################################################################
save.dir <- "Output/03_10_graphics_states_flow_flow_Spl_250508"
dir.create(save.dir)

#temporary
taxa_order <- list(Fungi=readRDS(sprintf("%s/Heatmap_taxa_order_fungi.rds",dir_02_06_02)),
                   Prokaryote=readRDS("%s/Heatmap_taxa_order_prokaryote.rds",dir_02_06_02))
##

full_sa <- list(Fungi=readRDS(sprintf("%s/Fungi/runSA_ST_full_Fungi_Family.rds",dir_02_06)),
                Prokaryote=readRDS(sprintf("%s/Prokaryota/runSA_ST_full_Prokaryota_Family.rds",dir_02_06)))


full_ss <-  list(Fungi=readRDS(sprintf("%s/Fungi/SScomposition_Fungi.rds",dir_02_06)),
                 Prokaryote=readRDS(sprintf("%s/Prokaryota/SScomposition_Prokaryota.rds",dir_02_06)))


coord_pcoa <- list(Fungi=readRDS(sprintf("%s/Fungi/coord_PCoA.rds",dir_02_07)),
                   Prokaryote=readRDS(sprintf("%s/Prokaryota/coord_PCoA.rds",dir_02_07)))

ocmat <- list(Fungi=readRDS(sprintf("%s/ELA_input_ocmat_Fungi.rds",dir_02_05)),
              Prokaryote=readRDS(sprintf("%s/ELA_input_ocmat_Prokaryote.rds",dir_02_05)))

mnmx_list <- list.files(dir_03_05,
                        "MaxMin_ab_",full.names = TRUE,recursive = TRUE)

l_pcoa <- list.files(pattern = "PCoA_new_SS_",
                      dir_03_05,full.names=TRUE,recursive = TRUE)

f_info <- t(sapply(strsplit(l_pcoa[grep("Fungi",l_pcoa)],"/"),function(x){
  c(x[12],gsub("PCoA_new_SS_Fungi_|.rds","",x[13]))
}))

df_pcoa_f <- foreach(nl=1:sum(grepl("Fungi",l_pcoa)),.combine=rbind) %do% {
  l <- readRDS(l_pcoa[grep("Fungi",l_pcoa)][nl])
  return(cbind(l,Taxa=f_info[nl,1],plant=f_info[nl,2]))
}

p_info <- t(sapply(strsplit(l_pcoa[grep("Prok",l_pcoa)],"/"),function(x){
  c(x[12],gsub("PCoA_new_SS_Prokaryota_|.rds","",x[13]))
}))

df_pcoa_p <- foreach(nl=1:sum(grepl("Prok",l_pcoa)),.combine=rbind) %do% {
  l <- readRDS(l_pcoa[grep("Prok",l_pcoa)][nl])
  return(cbind(l,Taxa=p_info[nl,1],plant=p_info[nl,2]))
}

#evenness
e_change <- list.files(dir_03_05,
                       "evenness_change_",full.names = TRUE,recursive = TRUE)

targ <- list.files(pattern = "Flow_data",dir_03_05,full.names=TRUE,recursive = TRUE)
ssims <- lapply(targ,readRDS)

names(ssims) <- gsub("Prokaryota","Prokaryote",sapply(strsplit(split = "/",targ),
                       function(x)paste(x[11],x[12],sep="_")))

sa_list <- list.files(dir_03_01,pattern = "runSA_",full.names = TRUE,recursive = TRUE)

topmat <- rbind(read.csv(sprintf("%s/delta_landscape_topmat.csv",dir_03_07)),
                 read.csv(sprintf("%s/delta_landscape_topmat_e.csv",dir_03_07)))

pick_taxa_f <- unique(topmat[which(topmat$z_land>0&topmat$fb=="Fungi"&topmat$ra=="median"),"Taxa"])

pf <- foreach(tx=pick_taxa_f,.combine=rbind) %do% {
  #tx <- pick_taxa_f[1]
  tm <- topmat[which(topmat$Taxa==tx&topmat$fb=="Fungi"&topmat$ra=="median"),]
  cbind(tx,tm[which.max(tm$z_land),"plant"])
}

pick_taxa_p <- unique(topmat[which(topmat$z_land>0&topmat$fb=="Prokaryote"&topmat$ra=="median"),"Taxa"])

pp <- foreach(tx=pick_taxa_p,.combine=rbind) %do% {
  #tx <- pick_taxa_p[1]
  tm <- topmat[which(topmat$Taxa==tx&topmat$fb=="Prokaryote"&topmat$ra=="median"),]
  cbind(tx,tm[which.max(tm$z_land),"plant"])
}

top2mat1 <-rbind(cbind(topmat[topmat$fb=="Fungi",],Taxa2=paste0("Fungi_",topmat[topmat$fb=="Fungi",]$Taxa),
                      Main=apply(topmat[topmat$fb=="Fungi",],1,function(x){
                        ifelse(x["Taxa"] %in% pf[,1],
                               x["plant"] == pf[pf[,1] == x["Taxa"],2] & x["ra"] == "median",
                               FALSE)
                        })),
                cbind(topmat[topmat$fb=="Prokaryote",],Taxa2=paste0("Prokaryote_",topmat[topmat$fb=="Prokaryote",]$Taxa),
                      Main=apply(topmat[topmat$fb=="Prokaryote",],1,function(x){
                        ifelse(x["Taxa"] %in% pp[,1],
                               x["plant"] == pp[pp[,1] == x["Taxa"],2] & x["ra"] == "median",
                               FALSE)
                      })))
                      
                      
top2mat <- top2mat1[which(top2mat1$Main==TRUE),]

sss <- NULL
sim1 <- list(NULL)
msmats <- NULL
rep_sts <- list(NULL)
for(i in 1:nrow(top2mat)){
  cat(i,"/",nrow(top2mat),"\r")
#i <- 23  
  ocm <- ocmat[[top2mat[i,"fb"]]]
  sa1 <- readRDS(sa_list[intersect(grep(ifelse(top2mat[i,"fb"]=="Fungi","Fungi","Prokaryota"),sa_list),
                                  grep(gsub("Fungi_|Prokaryote_","",top2mat[i,"Taxa2"]),sa_list))])
  
  fsa <- rownames(full_sa[[top2mat[i,"fb"]]][[1]])
  
  sim1[[i]] <- as.data.frame(ssims[[top2mat[i,"Taxa2"]]][[top2mat[i,"plant"]]])
  
  #a <- readRDS("/Users/ngch/Desktop/Now_using/data/sugaFB_ELA/analysis_series_SpCom/02_ELA_241216/tmp/Russula/Flow_data_Fungi_Russula.rds")[["Acer"]]
  #sim1[[i]] <- as.data.frame(a)
  maj_ss0 <-c()
 for(st in 1:length(unique(sim1[[i]]$time))){#st <- 1
   summ <- sim1[[i]][sim1[[i]]$time == unique(sim1[[i]]$time)[st],]%>%
     dplyr::group_by(group) %>% dplyr::summarise_each(dplyr::funs(sum),
                                                      vars = c(count_total = count))
   ssfreq <- summ$vars...count_total/sum(summ$vars...count_total)
   
   maj_ss0 <- c(maj_ss0,summ$group[ssfreq>0.05])
 } 
 
  maj_ss1 <- unique(maj_ss0)
 
  #pickup top5 dominant states
  th <- 2
  summ_1 <- sim1[[i]][sim1[[i]]$time == 1,]%>%
    dplyr::group_by(group) %>% dplyr::summarise_each(dplyr::funs(sum),
                                                     vars = c(count_total = count))
  summ_16 <- sim1[[i]][sim1[[i]]$time == 16,]%>%
    dplyr::group_by(group) %>% dplyr::summarise_each(dplyr::funs(sum),
                                                     vars = c(count_total = count))
  summ_32 <- sim1[[i]][sim1[[i]]$time == 32,]%>%
    dplyr::group_by(group) %>% dplyr::summarise_each(dplyr::funs(sum),
                                                     vars = c(count_total = count))
  
  if(nrow(summ_1)>th){
    dm_sts_1 <- gsub("NA",0,summ_1$group[order(summ_1$vars...count_total,decreasing = TRUE)[1:th]])
  }else{
    dm_sts_1 <- gsub("NA",0,summ_1$group)
  }
  
  if(nrow(summ_16)>th){
    dm_sts_16 <- gsub("NA",1,summ_16$group[order(summ_16$vars...count_total,decreasing = TRUE)[1:th]])
  }else{
    dm_sts_16 <- gsub("NA",1,summ_16$group)
  }
  
  if(nrow(summ_32)>th){
    dm_sts_32 <- gsub("NA",1,summ_32$group[order(summ_32$vars...count_total,decreasing = TRUE)[1:th]])
  }else{
    dm_sts_32 <- gsub("NA",1,summ_32$group)
  }
  
  dm_sts <- t(sapply(strsplit(unique(c(dm_sts_1,dm_sts_16,dm_sts_32)),""),as.numeric))
  dimnames(dm_sts) <- list(unique(c(dm_sts_1,dm_sts_16,dm_sts_32)),fsa)
  ####
  
  msmat1 <- unique(rbind(cbind(state_id=maj_ss1,
                  T(sapply(strsplit(gsub("A","",maj_ss1),""),function(x){
               if(all(x!="N")){return(as.numeric(x))}else{
                 x[x=="N"] <- 0
                 return(as.numeric(x))
               }}
                  ))),
               cbind(state_id=maj_ss1,
                     T(sapply(strsplit(gsub("A","",maj_ss1),""),function(x){
                       if(all(x!="N")){return(as.numeric(x))}else{
                         x[x=="N"] <- 1
                         return(as.numeric(x))
                       }}
                     )))))
  
  msmat0 <- unique(rbind(cbind(state_id=maj_ss1,
                               T(sapply(strsplit(gsub("A","",maj_ss1),""),function(x){
                                 if(all(x!="N")){return(as.numeric(x))}else{
                                   x[x=="N"] <- NA
                                   return(as.numeric(x))
                                 }}
                               )))))
  
  msmat0$state_id <- apply(msmat0[,-1],1,paste0,collapse="")
  dimnames(msmat0) <- list(msmat0$state_id,c("state_id",fsa))
  
  msmat1$state_id <- apply(msmat1[,-1],1,paste0,collapse="")
  dimnames(msmat1) <- list(msmat1$state_id,c("state_id",fsa))
  
  
  sss <- rbind(sss,cbind(top2mat[i,"fb"],c(msmat1$state_id)))
  
  ##outputs
  msmats[[top2mat[i,"fb"]]] <- rbind(msmats[[top2mat[i,"fb"]]],msmat1)
  rep_sts[[i]] <-  dm_sts
}


msmats2 <-  lapply(msmats,function(x){unique(x)})

## Fungi
set.seed(123)
fb <- "Fungi"
msm <- msmats2[["Fungi"]][,-1];rownames(msm) <- msmats2[["Fungi"]][,1]
mss <-  bind_rows(as.data.frame(full_ss[[fb]]),msm)[,colnames(msm)]
mss[is.na(mss)] <- 0

dmat_f <- vegdist(mss,method="jaccard")

hcf <- hclust(as.dist(dmat_f),method="ward.D2")

 
sf2 <- cbind(fb="Fungi",c(rownames(full_ss[[fb]]),msmats2[[fb]][,1])[hcf$order],
                 color[1:length(hcf$order)],
                 ifelse(grepl("F_B",c(rownames(full_ss[[fb]]),msmats2[[fb]][,1])[hcf$order]),
                        c(rownames(full_ss[[fb]]),msmats2[[fb]][,1])[hcf$order],
                        c(rep(NA,nrow(full_ss[[fb]])),
                        paste0("F_pB",formatC(1:length(msmats2[[fb]][,1]),width=2,flag=0)))[hcf$order]))

rownames(sf2) <- sf2[,2]

########################
fb <- "Prokaryote"
msm2 <- msmats2[["Prokaryote"]][,-1]#[,-1];rownames(msm2) <- msmats2[["Prokaryote"]][,1]
mss2 <- bind_rows(as.data.frame(full_ss[[fb]]),msm2)[,colnames(msm2)]
mss2[is.na(mss2)] <- 0

dmat_b <- vegdist(mss2,method="jaccard")

hcb <- hclust(as.dist(dmat_b),method="ward.D2")

sb2 <- cbind(fb="Prokaryote",
             c(rownames(full_ss[[fb]]),
               msmats2[["Prokaryote"]][,1])[hcb$order],
                 color[length(hcb$order)+1:length(hcb$order)],
             ifelse(grepl("P_B",c(rownames(full_ss[[fb]]),msmats2[[fb]][,1])[hcb$order]),
                    c(rownames(full_ss[[fb]]),msmats2[[fb]][,1])[hcb$order],
                    ifelse(grepl("P_B",c(rownames(full_ss[[fb]]),msmats2[[fb]][,1])[hcb$order]),
                           c(rownames(full_ss[[fb]]),msmats2[[fb]][,1])[hcb$order],
                           c(rep(NA,nrow(full_ss[[fb]])),
                             paste0("P_pB",formatC(1:length(msmats2[[fb]][,1]),width=2,flag=0)))[hcb$order])))
             
rownames(sb2) <- sb2[,2]

# #merge & set colors
ms <- rbind(sf2,sb2)#;rownames(ms) <- ms[,4]

#Fungi
fb <- "Fungi"
msa_f <- mss#[,which(colSums(mss,na.rm = TRUE)>0)]

#####

msa_f2 <- msa_f[grep("F_B",rownames(msa_f),invert=TRUE),]

####
alpha_ss_col <- 0.8

###skip#####
#Fungi


mds_result_f <- cbind(state_id=rownames(mss),
                    as.data.frame(cmdscale(dmat_f, k = 3)))  # 3次元空間にマッピング

colnames(mds_result_f) <- c("state_id","PCo1","PCo2","PCo3")
#mds_result <- cbind(state_id=rownames(df_pcoa2),as.data.frame(df_pcoa2))
#mds_result$ssid[is.na(mds_result$ssid)] <- paste0("F_pB",1:sum(is.na(mds_result$ssid))) 

# RGB色を生成
# 各次元を0〜1の範囲に正規化
mds_scaled_f <- apply(mds_result_f[,grep("PCo",colnames(mds_result_f))],
                      2, function(x) (x - min(x)) / (max(x) - min(x)))

# RGB色を生成
colors_f <- rgb(mds_scaled_f[, "PCo3"], mds_scaled_f[, "PCo1"], mds_scaled_f[, "PCo2"],alpha=alpha_ss_col)

names(colors_f) <- rownames(msa_f)

plot(mds_result_f[, "PCo1"], mds_result_f[, "PCo2"], col = colors_f, pch = 19, cex = 2,
     xlab = "Dimension 1", ylab = "Dimension 2", main = "Color Mapping Based on Dissimilarity")

###skip end####

#Prokaryote
fb <- "Prokaryote"
msa_p <- mss2#[,which(colSums(mss2,na.rm = TRUE)>0)]

#####

msa_p2 <- msa_p[grep("P_B",rownames(msa_p),invert=TRUE),]


##skip###
####
mds_result <- cbind(state_id=rownames(mss2),as.data.frame(cmdscale(dmat_b, k = 3)))  # 3次元空間にマッピング
colnames(mds_result) <- c("state_id","PCo1","PCo2","PCo3")
#mds_result$ssid[is.na(mds_result$ssid)] <- paste0("F_pB",1:sum(is.na(mds_result$ssid))) 

# RGB色を生成
# 各次元を0〜1の範囲に正規化
mds_scaled <- apply(mds_result[,grep("PCo",colnames(mds_result))], 2, function(x) (x - min(x)) / (max(x) - min(x)))

# RGB色を生成
colors_p <- rgb(mds_scaled[, "PCo3"], mds_scaled[, "PCo1"], mds_scaled[, "PCo2"],alpha=alpha_ss_col)

names(colors_p) <- rownames(msa_p)

plot(mds_result[, "PCo1"], mds_result[, "PCo2"], col = colors_p, pch = 19, cex = 2,
     xlab = "Dimension 1", ylab = "Dimension 2", main = "Color Mapping Based on Dissimilarity")

###
colvec <- c(colors_f,
            colors_p,
            "gray80") 
names(colvec) <- c(sf2[names(colors_f),4],sb2[names(colors_p),4],"Minor basins")

saveRDS(colvec,sprintf("%s/states_colvector.rds",save.dir))

##make pcoa data frame
df_pcoa <- df_pcoa_f
df_pcoa$state_id[grep("F_B",rownames(df_pcoa),invert=TRUE)] <- ifelse(!grepl("NA",rownames(df_pcoa[grep("F_B",rownames(df_pcoa),invert=TRUE),])),
                                                                      substr(rownames(df_pcoa[grep("F_B",rownames(df_pcoa),invert=TRUE),]),1,min(nchar(rownames(msa_f2)))),
                                                                      substr(rownames(df_pcoa[grep("F_B",rownames(df_pcoa),invert=TRUE),]),1,min(nchar(rownames(msa_f2)))+1)
)

df_pcoa$state_id <- ifelse(df_pcoa$step==1,gsub("NA","0",df_pcoa$state_id),gsub("NA","1",df_pcoa$state_id))
df_pcoa$state_id[grep("F_B",rownames(df_pcoa))] <-df_pcoa$ssid[grep("F_B",rownames(df_pcoa))] 

mds <- merge(df_pcoa[,grep("PCo",colnames(df_pcoa),invert = TRUE)],
             mds_result_f[,c("PCo1","PCo2","PCo3","state_id")],
             by="state_id")[,c("state_id","step","SS","ssid","Taxa","plant","PCo1","PCo2","PCo3")]
mds2 <- cbind(mds,
              group=ifelse(mds$state_id %in% c("F_B1","F_B2","F_B3"),"Group 1 (EcM)",
                           ifelse(mds$state_id %in% c("F_B4","F_B5","F_B7"),"Group 2 (AM)","Others")),
              ssid2=sf2[mds$state_id,4],
              full=factor(ifelse(grepl("F_B",mds$ssid),"In whole commnunity stability landscape",
                                 "In stability landscape without analyzed genus"),
                          levels=c("In whole commnunity stability landscape",
                                   "In stability landscape without analyzed genus")))

saveRDS(mds2,sprintf("%s/gradland_pcoa_all_Fungi.rds",save.dir))

g_full <- ggplot(mds2,aes(x=PCo1,y=PCo2))+
  geom_point(data=function(x){x[grep("F_B",x$ssid2),]},
             aes(fill=ssid2),shape=22,size=2)+
  geom_point(data=function(x){x[grep("F_pB",x$ssid2),]},
             aes(shape=SS),fill="white",size=4,alpha=0.5)+
  #geom_text_repel(aes(label=ssid),size=2)+
  labs(fill="Stable states",shape="",x="PCoA 1", y="PCoA 2")+
  theme_bw()+
  theme(legend.position = "bottom",
        legend.direction = "vertical")+
  scale_fill_manual(values=colvec)+
  scale_shape_manual(values=c(fullSS=22,partialSS=21))+
  guides(fill=guide_legend(nrow=1))

ggsave(sprintf("%s/legend_PCoA_Fungi.pdf",save.dir),
       g_legend(g_full),h=4,w=8)



##make pcoa data frame
df_pcoa <- df_pcoa_p
df_pcoa$state_id[grep("P_B",rownames(df_pcoa),invert=TRUE)] <- ifelse(!grepl("NA",rownames(df_pcoa[grep("P_B",rownames(df_pcoa),invert=TRUE),])),
                                                                      substr(rownames(df_pcoa[grep("P_B",rownames(df_pcoa),invert=TRUE),]),1,min(nchar(rownames(msa_p2)))),
                                                                      substr(rownames(df_pcoa[grep("P_B",rownames(df_pcoa),invert=TRUE),]),1,min(nchar(rownames(msa_p2)))+1)
)

df_pcoa$state_id <- ifelse(df_pcoa$step==1,gsub("NA","0",df_pcoa$state_id),gsub("NA","1",df_pcoa$state_id))

df_pcoa$state_id[grep("P_B",rownames(df_pcoa))] <-df_pcoa$ssid[grep("P_B",rownames(df_pcoa))] 

mds <- merge(df_pcoa[,grep("PCo",colnames(df_pcoa),invert = TRUE)],
             mds_result[,c("PCo1","PCo2","PCo3","state_id")],
             by="state_id")[,c("state_id","step","SS","ssid","Taxa","plant","PCo1","PCo2","PCo3")]
mds2_p <- cbind(mds,
                ssid2=sb2[mds$state_id,4],
                group=ifelse(mds$state_id %in% paste0("P_B",2:8),"Group 3",
                             ifelse(mds$state_id %in% paste0("P_B",9:15),"Group 2",
                                    ifelse(mds$state_id %in% paste0("P_B",16:24),"Group 1","Others"))),
              full=factor(ifelse(grepl("P_B",mds$ssid),"In whole commnunity stability landscape",
                                 "In stability landscape without analyzed genus"),
                          levels=c("In whole commnunity stability landscape",
                                   "In stability landscape without analyzed genus")))

saveRDS(mds2_p,sprintf("%s/gradland_pcoa_all_Prokaryote.rds",save.dir))

mds2_p$ssid3 <- factor(mds2_p$ssid2,
                       levels=c(paste0("P_B",1:sum(grepl("P_B",mds2_p$ssid2))),
                                paste0("P_pB",1:sum(grepl("P_pB",mds2_p$ssid2)))))
g_full <- ggplot(mds2_p,aes(x=PCo1,y=PCo2))+
  geom_point(data=function(x){x[grep("P_B",x$ssid2),]},
             aes(fill=ssid3),shape=22,size=2)+
  geom_point(data=function(x){x[grep("P_pB",x$ssid2),]},
             aes(shape=SS),fill="white",size=4)+
  #geom_text_repel(aes(label=ssid),size=2)+
  labs(fill="Stable states",shape="",x="PCoA 1",y="PCoA 2")+
  theme_bw()+
  theme(legend.position = "bottom",
        legend.direction = "vertical")+
  scale_fill_manual(values=colvec)+
  scale_shape_manual(values=c(fullSS=22,partialSS=21))+
  guides(fill=guide_legend(nrow=3))


ggsave(sprintf("%s/legend_PCoA_Prokaryotes.pdf",save.dir),
       g_legend(g_full),h=4,w=8)

#####################################
#fb <- "Fungi"
ga1 <- list(NULL)
mss_all_f <- NULL
mss_all_p <- NULL
dir.create(sprintf("%s/Top",save.dir),showWarnings = FALSE)
#dir.create(sprintf("%s/Top2",save.dir),showWarnings = FALSE)
dir.create(sprintf("%s/Top/Fungi",save.dir),showWarnings = FALSE)
dir.create(sprintf("%s/Top/Prokaryote",save.dir),showWarnings = FALSE)
#dir.create(sprintf("%s/Top2/Prokaryote",save.dir),showWarnings = FALSE)

for(i in 1:nrow(top2mat)){
 #i <- 209
  mnmx1 <- readRDS(mnmx_list[intersect(grep(ifelse(top2mat[i,"fb"]=="Fungi","Fungi","Prokaryota"),mnmx_list),
                                       grep(gsub("Fungi_|Prokaryote_","",top2mat[i,"Taxa2"]),mnmx_list))])
  
 
  
  dir.create(sprintf("%s/Top/%s/%s",
                     save.dir,top2mat[i,"fb"],
                     gsub("Fungi_|Prokaryote_","",top2mat[i,3])),
             showWarnings = FALSE)
  
  
  if(top2mat[i,"fb"]=="Fungi"){
    s2 <- sf2
    mds3 <- mds2
  }else{
      s2 <- sb2
      mds3 <- mds2_p
    }
  
  dat2 <- sim1[[i]]
  
  if(any(grepl("NA",dat2$group))){
    dat2$group[grep("NA",dat2$group)] <- ifelse(dat2$time==1,
                                                gsub("NA","0",dat2$group),
                                                gsub("NA","1",dat2$group))
  }
  
  clrab <- round(seq(from=mnmx1[[top2mat[i,"plant"]]][1],
      to=mnmx1[[top2mat[i,"plant"]]][2],length.out=32),3)
  
  dat2$alluvium2 <- sapply(strsplit(dat2$alluvium," "),
         function(x){paste(clrab[as.numeric(x[1])],
                           clrab[as.numeric(x[2])])})
  
  dat2$time2 <- ifelse(dat2$time==round(dat2$time,1),
                       clrab[dat2$time],clrab[dat2$time]+0.000001)
  dat2$ss <- ifelse(dat2$group %in% c(sf2[,2],sb2[,2]),
                    dat2$group, "Minor basins")
  
  dat2$ss <- ifelse(dat2$ss!= "Minor basins",
                   ms[,4][dat2$ss], "Minor basins")
  
  
  g_FLOW<- dat2 %>%
    filter(count>0) %>%
    ggplot(aes(x = time, y = count, stratum = group, alluvium = alluvium)) +
    geom_flow(aes(fill = ss),width =0.025, alpha = 1) +
    geom_stratum(data = dat2[dat2$time%%1==0,], aes(fill = ss), alpha = 1,linewidth=0.1) +
    labs(y="Prorortion of basins\n( 20,000 initial states )",
         x="Relative read count + 1 ( CLR-transformed )",
         subtitle = sprintf("*%s* in *%s*",gsub("Fungi_|Prokaryote_","",top2mat[i,"Taxa2"]),top2mat[i,"plant"]))+
    scale_fill_manual(values=colvec) +
    scale_x_continuous(breaks = seq(1, 32, 6),
                       labels = seq(clrab[1], clrab[32], length.out=6))+
    #geom_text(data = dat2[dat2$time%%1==0,], stat = "stratum", aes(label = after_stat(stratum))) +
    theme_bw()+
    theme(legend.position = "none",
          plot.subtitle = element_markdown(size = 15,hjust=0.5),
          plot.title = element_markdown(size = 16,hjust=0.5),
          axis.title = element_text(size=11),
          axis.text = element_text(size=8))
  
  
  ggsave(sprintf("%s/Top/%s/%s/Flow_diagram1_%s_%s.pdf",save.dir,
                 top2mat[i,"fb"],
                 gsub("Fungi_|Prokaryote_","",top2mat[i,"Taxa2"]),
                 gsub("Fungi_|Prokaryote_","",top2mat[i,"Taxa2"]),
                 top2mat[i,2]),g_FLOW,
         h=3,w=8)
  
  
  #PCoA of basins
  gobj <- mds3[mds3$Taxa==top2mat$Taxa[i] & mds3$plant==top2mat$plant[i],]
  
  sts <- rep_sts[[i]]
  
  if(top2mat[i,"fb"]=="Fungi"){
    gobj$rep_sts <- gobj$ssid2 %in% sf2[rownames(sts)[rownames(sts) %in% rownames(sf2)],4]
  }else{
    gobj$rep_sts <- gobj$ssid2 %in% sb2[rownames(sts)[rownames(sts) %in% rownames(sb2)],4]
  }
  
  g_PCOA <- ggplot(gobj,aes(x=PCo1,y=PCo2))+
    geom_point(data=function(x){x[!is.na(x$ssid),]},
               aes(fill=ssid2,shape=SS,size=SS))+
    geom_point(data=function(x){x[is.na(x$ssid),]},
               aes(fill=ssid2,shape=SS,size=SS))+
    #geom_text_repel(aes(label=ssid),size=2)+
    labs(subtitle =  sprintf("*%s* in *%s*",gsub("Fungi_|Prokaryote_","",top2mat[i,"Taxa2"]),top2mat[i,"plant"]),
         x="PCoA 1",y="PCoA 2")+
    theme_bw()+
    theme(#aspect.ratio=1,
          strip.text = element_blank(),
          plot.subtitle = element_markdown(size = 15,hjust=0.5),
          strip.background= element_blank(),
          plot.title = element_text(size=20),
          legend.position = "none",
          axis.title = element_text(size=11),
          axis.text = element_text(size=8))+
    scale_fill_manual(values=colvec)+
    scale_shape_manual(values=c(fullSS=22,partialSS=21))+
    scale_size_manual(values=c(fullSS=2,partialSS=4))+
    facet_wrap(~step)
  
  ggsave(sprintf("%s/Top/%s/%s/PCoA_%s_%s_nolab.pdf",save.dir,
                 top2mat[i,"fb"],
                 gsub("Fungi_|Prokaryote_","",top2mat[i,"Taxa2"]),
                 gsub("Fungi_|Prokaryote_","",top2mat[i,"Taxa2"]),
                 top2mat[i,"plant"]),g_PCOA+theme(aspect.ratio = 1),
         h=3,w=8)
  ggsave(sprintf("%s/Top/%s/%s/PCoA_%s_%s.pdf",save.dir,
                 top2mat[i,"fb"],
                 gsub("Fungi_|Prokaryote_","",top2mat[i,"Taxa2"]),
                 gsub("Fungi_|Prokaryote_","",top2mat[i,"Taxa2"]),
                 top2mat[i,"plant"]),
         g_PCOA+geom_text_repel(aes(label=ssid),size=1)+
           geom_text_repel(data = function(x){x[x$rep_sts,]},aes(label=ssid2),size=2)+
           theme(aspect.ratio = 1),
         h=3,w=8)
  
  #evenness
  
  ec <- readRDS(e_change[grepl(ifelse(top2mat[i,"fb"]=="Fungi","Fungi","Prokaryota"),
                               e_change)&grepl(top2mat$Taxa[i],e_change,fixed = TRUE)])
  
  mine <- min(ec[ec$plant==top2mat$plant[i],"evenness"][1:32])
  maxe <- max(ec[ec$plant==top2mat$plant[i],"evenness"][1:32])
  
  g_even <- ggplot(ec[ec$plant==top2mat$plant[i],][1:32,],
                   aes(x=ra_clr,y=evenness))+
    geom_line(color="darkorange",linewidth=3)+
     geom_point(data=function(x){x[which(!order(x$ra_clr) %in% c(1,16,32)),]},
                fill="transparent",size=3,shape=21)+
    geom_point(data=function(x){x[which(order(x$ra_clr) %in% c(1,16,32)),]},
               fill="red",size=4,shape=23)+
    labs(y="Stable state\nevenness",x="Relative read count + 1 ( CLR-transformed )",
         subtitle =  sprintf("*%s* in *%s*",gsub("Fungi_|Prokaryote_","",top2mat[i,"Taxa2"]),top2mat[i,"plant"]))+
    theme_bw()+
    theme(#aspect.ratio=0.3,
      axis.title.x = element_blank(),
      plot.subtitle = element_markdown(size = 15,hjust=0.5),
      axis.text.x =  element_blank(),
          axis.title.y = element_text(size=11),
          axis.text = element_text(size=9))+
    coord_cartesian(ylim=c(mine-0.06*(maxe-mine),
                           maxe+0.06*(maxe-mine)))
 
  ggsave(sprintf("%s/Top/%s/%s/evenness_change_%s_%s.pdf",save.dir,
                 top2mat[i,"fb"],
                 gsub("Fungi_|Prokaryote_","",top2mat[i,"Taxa2"]),
                 gsub("Fungi_|Prokaryote_","",top2mat[i,"Taxa2"]),
                 top2mat[i,"plant"]),g_even+theme(aspect.ratio = 0.3),
         h=3,w=8)
  
  ##merge
  ga1[[i]]<- g_PCOA +# theme(panel.spacing = unit(1, "cm")) +
    g_even+theme(plot.subtitle = element_blank())+
    g_FLOW+theme(plot.subtitle = element_blank())+
    plot_layout(ncol=1,heights = c(2,1.3,2))#g_FLOW
  
  ggsave(sprintf("%s/Top/%s/%s/Flow_diagram_withPCoA_%s_%s_nolab.pdf",save.dir,
                 top2mat[i,"fb"],
                 gsub("Fungi_|Prokaryote_","",top2mat[i,"Taxa2"]),
                 gsub("Fungi_|Prokaryote_","",top2mat[i,"Taxa2"]),
                 top2mat[i,"plant"]),ga1[[i]],
         h=6,w=6)
  ggsave(sprintf("%s/Top/%s/%s/Flow_diagram_withPCoA_%s_%s.pdf",save.dir,
                 top2mat[i,"fb"],
                 gsub("Fungi_|Prokaryote_","",top2mat[i,"Taxa2"]),
                 gsub("Fungi_|Prokaryote_","",top2mat[i,"Taxa2"]),
                 top2mat[i,"plant"]),
         g_PCOA +geom_text_repel(aes(label=ssid),size=2)+# theme(panel.spacing = unit(1, "cm")) +
           g_even+theme(plot.subtitle = element_blank())+
           g_FLOW+theme(plot.subtitle = element_blank())+
           plot_layout(ncol=1,heights = c(2,1.3,2)),#g_FLOW,
         h=6,w=6)
  
}

##make dual mycorrhizal stable state composition heatmap


jac_dist <- function(x) {
    vegdist(x, method = "jaccard")

}

mss_dual <- mss[!grepl("F_B",rownames(mss)) & mss[,"Russulaceae"]==1 & mss[,"Glomeraceae"]==1,]
mss_ecm <- mss[!grepl("F_B",rownames(mss)) & mss[,"Russulaceae"]==1 & mss[,"Glomeraceae"]==0,]
mss_am <- mss[!grepl("F_B",rownames(mss)) & mss[,"Russulaceae"]==0 & mss[,"Glomeraceae"]==1,]

mss_dual1 <- mss_dual[,colSums(mss_dual)>0]
ht_ord <- draw(Heatmap(as.matrix(mss_dual1),
                       clustering_distance_rows = jac_dist,
                       clustering_distance_columns = jac_dist, 
                       clustering_method_rows = "ward.D2",
                       clustering_method_columns = "ward.D2"))

rord <- row_order(ht_ord)
cord <- column_order(ht_ord)

mss_dual2 <- mss_dual1[rord,cord]
cust_col <- foreach(i = 1:nrow(mss_dual2), .combine = rbind) %do% {
  ifelse(mss_dual2[i, ]==1,colors_f[rownames(mss_dual2)[i]],"white")
} 
rownames(cust_col) <- rownames(mss_dual2)


col_fun = colorRamp2(c(0, 1), c("white", "white"))


sscomp <- Heatmap(mss_dual2,
                  row_names_side = "left", column_names_side = "bottom", 
                  col = col_fun,
                  show_column_dend = FALSE,
                  show_row_names = FALSE,
                  show_heatmap_legend = FALSE,
                  row_dend_side = "left", 
                  rect_gp = gpar(fill = as.vector(cust_col),
                                 col = "black"), 
                  row_names_gp = gpar(cex=0.9), 
                  column_names_gp = gpar(cex=0.9), 
                  row_dend_width = unit(2, "cm"),
                  column_dend_height = unit(2, "cm"), 
                  column_title = "Family composition of dual mycorrhizal stable states",
                  column_names_rot = 45, row_title = "Stable states", 
                  clustering_distance_rows = jac_dist,
                  clustering_distance_columns = jac_dist, 
                  clustering_method_rows = "ward.D2",
                  clustering_method_columns = "ward.D2")

pdf(sprintf("%s/SScomposition_dual_mycorrhizal.pdf",save.dir),
    width = 8, height = 5)
sscomp
dev.off()


g_dual <- ggplot(mds2,aes(x=PCo1,y=PCo2))+
  geom_point(data=function(x){x[grep("F_B",x$ssid2),]},
             aes(fill=ssid2),shape=22,size=2,
             show.legend = FALSE)+
  geom_point(data=function(x){x[grep("F_pB",x$ssid2),]},
             aes(shape=SS),fill="white",size=4,alpha=0.5,
             show.legend = FALSE)+
  geom_point(data=function(x){x[x$state_id %in% rownames(mss_dual2),]},
             aes(shape=SS,fill=ssid2),size=4,alpha=0.5,
             show.legend = FALSE)+
  # geom_text_repel(data=function(x){x[x$state_id %in% rownames(mss_dual2),]},
  #                 aes(label=ssid2),size=2)+
  labs(fill="Stable states",shape="",x="PCoA 1", y="PCoA 2",title="Dual muycorrhizal states")+
  theme_bw()+
  theme(legend.position = "bottom",
        aspect.ratio = 1,
        axis.text = element_text(size=12),
        axis.title = element_text(size=14),
        legend.direction = "vertical")+
  scale_fill_manual(values=colvec)+
  scale_shape_manual(values=c(fullSS=22,partialSS=21))+
  guides(fill=guide_legend(nrow=1))


g_ecm <- ggplot(mds2,aes(x=PCo1,y=PCo2))+
  geom_point(data=function(x){x[grep("F_B",x$ssid2),]},
             aes(fill=ssid2),shape=22,size=2,
             show.legend = FALSE)+
  geom_point(data=function(x){x[grep("F_pB",x$ssid2),]},
             aes(shape=SS),fill="white",size=4,alpha=0.5,
             show.legend = FALSE)+
  geom_point(data=function(x){x[x$state_id %in% rownames(mss_ecm),]},
             aes(shape=SS,fill=ssid2),size=4,alpha=0.5,
             show.legend = FALSE)+
  # geom_text_repel(data=function(x){x[x$state_id %in% rownames(mss_dual2),]},
  #                 aes(label=ssid2),size=2)+
  labs(fill="Stable states",shape="",x="PCoA 1", y="PCoA 2", title="EcM states")+
  theme_bw()+
  theme(legend.position = "bottom",
        aspect.ratio = 1,
        axis.text = element_text(size=12),
        axis.title = element_text(size=14),
        legend.direction = "vertical")+
  scale_fill_manual(values=colvec)+
  scale_shape_manual(values=c(fullSS=22,partialSS=21))+
  guides(fill=guide_legend(nrow=1))


g_am <- ggplot(mds2,aes(x=PCo1,y=PCo2))+
  geom_point(data=function(x){x[grep("F_B",x$ssid2),]},
             aes(fill=ssid2),shape=22,size=2,
             show.legend = FALSE)+
  geom_point(data=function(x){x[grep("F_pB",x$ssid2),]},
             aes(shape=SS),fill="white",size=4,alpha=0.5,
             show.legend = FALSE)+
  geom_point(data=function(x){x[x$state_id %in% rownames(mss_am),]},
             aes(shape=SS,fill=ssid2),size=4,alpha=0.5,
             show.legend = FALSE)+
  # geom_text_repel(data=function(x){x[x$state_id %in% rownames(mss_dual2),]},
  #                 aes(label=ssid2),size=2)+
  labs(fill="Stable states",shape="",x="PCoA 1", y="PCoA 2", title="AM states")+
  theme_bw()+
  theme(legend.position = "bottom",
        aspect.ratio = 1,
        axis.text = element_text(size=12),
        axis.title = element_text(size=14),
        legend.direction = "vertical")+
  scale_fill_manual(values=colvec)+
  scale_shape_manual(values=c(fullSS=22,partialSS=21))+
  guides(fill=guide_legend(nrow=1))

g_am

g_merge <- g_am+g_ecm+g_dual

ggsave(sprintf("%s/SS_PCoA_Fungi.pdf",save.dir),
       g_merge,h=4,w=10)


###make dominat states heatmap

jac_dist <- function(x) {
  vegdist(x, method = "jaccard")
}

i <- 4
  cat(paste0("Processing ",i," of ",length(rep_sts),"\n"))
  sts1 <- rep_sts[[i]]
  if(top2mat[i,"fb"]=="Fungi"){
    rownames(sts1) <- sf2[rownames(sts1)[rownames(sts1) %in% rownames(sf2)],4]
  }else{
    rownames(sts1) <- sb2[rownames(sts1)[rownames(sts1) %in% rownames(sb2)],4]
  }
  
  ht_ord <- draw(Heatmap(sts1,clustering_distance_rows = jac_dist,
                         clustering_method_rows = "ward.D2",
                         clustering_method_columns = "ward.D2"))
  
  rord <- row_order(ht_ord)
  
  sts4 <- sts1[rord,colSums(sts1)>0]
  
  i <- 5
  cat(paste0("Processing ",i," of ",length(rep_sts),"\n"))
  sts2 <- rep_sts[[i]]
  if(top2mat[i,"fb"]=="Fungi"){
    rownames(sts2) <- sf2[rownames(sts2)[rownames(sts2) %in% rownames(sf2)],4]
  }else{
    rownames(sts2) <- sb2[rownames(sts2)[rownames(sts2) %in% rownames(sb2)],4]
  }
  
  ht_ord <- draw(Heatmap(sts2,clustering_distance_rows = jac_dist,
                         clustering_method_rows = "ward.D2",
                         clustering_method_columns = "ward.D2"))
  
  rord <- row_order(ht_ord)
  
  sts5 <- sts2[rord,colSums(sts2)>0]
  

  #i = 4
  i <- 4
  
  sts <- sts1[,c(taxa_order[[top2mat[i,"fb"]]],
                 setdiff(union(colnames(sts4),colnames(sts5)),
                         taxa_order[[top2mat[i,"fb"]]]))]
  
  cust_col <- foreach(k = 1:nrow(sts), .combine = rbind) %do% {
    ifelse(sts[k, ]==1,colvec[rownames(sts)[k]],"white")
  } 
  rownames(cust_col) <- rownames(sts)
  
  
  col_fun = colorRamp2(c(0, 1), c("white", "white"))
  
  saveRDS(list(mat=sts,color=cust_col),sprintf("%s/Top/%s/%s/SSheatmap_%s_%s.rds",save.dir,top2mat[i,"fb"],
                  gsub("Fungi_|Prokaryote_","",top2mat[i,"Taxa2"]),
                  gsub("Fungi_|Prokaryote_","",top2mat[i,"Taxa2"]),
                  top2mat[i,"plant"]))
  
  sscomp <- Heatmap(sts, row_names_side = "left", column_names_side = "bottom", 
                    col = col_fun,
                    show_column_dend = FALSE,
                    show_row_dend = FALSE,
                    cluster_columns = FALSE,
                    
                    cluster_rows = FALSE,
                    show_heatmap_legend = FALSE,
                    row_dend_side = "left", 
                    rect_gp = gpar(fill = as.vector(cust_col),
                                   col = "black"), 
                    row_names_gp = gpar(cex=0.9), 
                    column_names_gp = gpar(cex=0.9), 
                    row_dend_width = unit(2, "cm"),
                    column_dend_height = unit(2, "cm"), 
                    column_names_rot = 45, row_title = "Stable states")
  
  
  pdf(sprintf("%s/Top/%s/%s/SS_heatmap_%s_%s.pdf",save.dir,top2mat[i,"fb"],
              gsub("Fungi_|Prokaryote_","",top2mat[i,"Taxa2"]),
              gsub("Fungi_|Prokaryote_","",top2mat[i,"Taxa2"]),
              top2mat[i,"plant"]),width = 5.5,height=3.7)
  
  sscomp
  
  dev.off() 
  
  i <- 5
  
  sts <- sts2[,c(taxa_order[[top2mat[i,"fb"]]],
                 setdiff(union(colnames(sts4),colnames(sts5)),
                         taxa_order[[top2mat[i,"fb"]]]))]
  
  cust_col <- foreach(k = 1:nrow(sts), .combine = rbind) %do% {
    ifelse(sts[k, ]==1,colvec[rownames(sts)[k]],"white")
  } 
  rownames(cust_col) <- rownames(sts)
  
  
  col_fun = colorRamp2(c(0, 1), c("white", "white"))
  
  saveRDS(list(mat=sts,color=cust_col),sprintf("%s/Top/%s/%s/SSheatmap_%s_%s.rds",save.dir,top2mat[i,"fb"],
                                               gsub("Fungi_|Prokaryote_","",top2mat[i,"Taxa2"]),
                                               gsub("Fungi_|Prokaryote_","",top2mat[i,"Taxa2"]),
                                               top2mat[i,"plant"]))
  
  sscomp <- Heatmap(sts, row_names_side = "left", column_names_side = "bottom", 
                    col = col_fun,
                    show_column_dend = FALSE,
                    show_row_dend = FALSE,
                    cluster_columns = FALSE,
                    
                    cluster_rows = FALSE,
                    show_heatmap_legend = FALSE,
                    row_dend_side = "left", 
                    rect_gp = gpar(fill = as.vector(cust_col),
                                   col = "black"), 
                    row_names_gp = gpar(cex=0.9), 
                    column_names_gp = gpar(cex=0.9), 
                    row_dend_width = unit(2, "cm"),
                    column_dend_height = unit(2, "cm"), 
                    column_names_rot = 45, row_title = "Stable states")
  
  
  pdf(sprintf("%s/Top/%s/%s/SS_heatmap_%s_%s.pdf",save.dir,top2mat[i,"fb"],
              gsub("Fungi_|Prokaryote_","",top2mat[i,"Taxa2"]),
              gsub("Fungi_|Prokaryote_","",top2mat[i,"Taxa2"]),
              top2mat[i,"plant"]),width = 5.5,height=3.5)
  
  sscomp
  
  dev.off() 
#}

  i <- 23
  cat(paste0("Processing ",i," of ",length(rep_sts),"\n"))
  stsp1 <- rep_sts[[i]]
  if(top2mat[i,"fb"]=="Fungi"){
    rownames(stsp1) <- sf2[rownames(stsp1)[rownames(stsp1) %in% rownames(sf2)],4]
  }else{
    rownames(stsp1) <- sb2[rownames(stsp1)[rownames(stsp1) %in% rownames(sb2)],4]
  }
  
  ht_ord <- draw(Heatmap(stsp1,clustering_distance_rows = jac_dist,
                         clustering_method_rows = "ward.D2",
                         clustering_method_columns = "ward.D2"))
  
  rord <- row_order(ht_ord)
  
  sts4 <- stsp1[rord,colSums(stsp1)>0]
  
  i <- 29
  cat(paste0("Processing ",i," of ",length(rep_sts),"\n"))
  stsp2 <- rep_sts[[i]]
  if(top2mat[i,"fb"]=="Fungi"){
    rownames(stsp2) <- sf2[rownames(stsp2)[rownames(stsp2) %in% rownames(sf2)],4]
  }else{
    rownames(stsp2) <- sb2[rownames(stsp2)[rownames(stsp2) %in% rownames(sb2)],4]
  }
  
  ht_ord <- draw(Heatmap(stsp2,clustering_distance_rows = jac_dist,
                         clustering_method_rows = "ward.D2",
                         clustering_method_columns = "ward.D2"))
  
  rord <- row_order(ht_ord)
  
  sts5 <- stsp2[rord,colSums(stsp2)>0]
  
  
  #i = 4
  i <- 23
  
  sts <- stsp1[,c(taxa_order[[top2mat[i,"fb"]]],
                 setdiff(union(colnames(sts4),colnames(sts5)),
                         taxa_order[[top2mat[i,"fb"]]]))]
  
  cust_col <- foreach(k = 1:nrow(sts), .combine = rbind) %do% {
    ifelse(sts[k, ]==1,colvec[rownames(sts)[k]],"white")
  } 
  rownames(cust_col) <- rownames(sts)
  
  
  col_fun = colorRamp2(c(0, 1), c("white", "white"))
  
  saveRDS(list(mat=sts,color=cust_col),sprintf("%s/Top/%s/%s/SSheatmap_%s_%s.rds",save.dir,top2mat[i,"fb"],
                                               gsub("Fungi_|Prokaryote_","",top2mat[i,"Taxa2"]),
                                               gsub("Fungi_|Prokaryote_","",top2mat[i,"Taxa2"]),
                                               top2mat[i,"plant"]))
  
  sscomp <- Heatmap(sts, row_names_side = "left", column_names_side = "bottom", 
                    col = col_fun,
                    show_column_dend = FALSE,
                    show_row_dend = FALSE,
                    cluster_columns = FALSE,
                    
                    cluster_rows = FALSE,
                    show_heatmap_legend = FALSE,
                    row_dend_side = "left", 
                    rect_gp = gpar(fill = as.vector(cust_col),
                                   col = "black"), 
                    row_names_gp = gpar(cex=0.9), 
                    column_names_gp = gpar(cex=0.9), 
                    row_dend_width = unit(2, "cm"),
                    column_dend_height = unit(2, "cm"), 
                    column_names_rot = 45, row_title = "Stable states")
  
  
  pdf(sprintf("%s/Top/%s/%s/SS_heatmap_%s_%s.pdf",save.dir,top2mat[i,"fb"],
              gsub("Fungi_|Prokaryote_","",top2mat[i,"Taxa2"]),
              gsub("Fungi_|Prokaryote_","",top2mat[i,"Taxa2"]),
              top2mat[i,"plant"]),width = 9,height=3.7)
  
  sscomp
  
  dev.off() 
  
  i <- 29
  
  sts <- stsp2[,c(taxa_order[[top2mat[i,"fb"]]],
                 setdiff(union(colnames(sts4),colnames(sts5)),
                         taxa_order[[top2mat[i,"fb"]]]))]
  
  cust_col <- foreach(k = 1:nrow(sts), .combine = rbind) %do% {
    ifelse(sts[k, ]==1,colvec[rownames(sts)[k]],"white")
  } 
  rownames(cust_col) <- rownames(sts)
  
  
  col_fun = colorRamp2(c(0, 1), c("white", "white"))
  
  saveRDS(list(mat=sts,color=cust_col),sprintf("%s/Top/%s/%s/SSheatmap_%s_%s.rds",save.dir,top2mat[i,"fb"],
                                               gsub("Fungi_|Prokaryote_","",top2mat[i,"Taxa2"]),
                                               gsub("Fungi_|Prokaryote_","",top2mat[i,"Taxa2"]),
                                               top2mat[i,"plant"]))
  
  sscomp <- Heatmap(sts, row_names_side = "left", column_names_side = "bottom", 
                    col = col_fun,
                    show_column_dend = FALSE,
                    show_row_dend = FALSE,
                    cluster_columns = FALSE,
                    
                    cluster_rows = FALSE,
                    show_heatmap_legend = FALSE,
                    row_dend_side = "left", 
                    rect_gp = gpar(fill = as.vector(cust_col),
                                   col = "black"), 
                    row_names_gp = gpar(cex=0.9), 
                    column_names_gp = gpar(cex=0.9), 
                    row_dend_width = unit(2, "cm"),
                    column_dend_height = unit(2, "cm"), 
                    column_names_rot = 45, row_title = "Stable states")
  
  
  pdf(sprintf("%s/Top/%s/%s/SS_heatmap_%s_%s.pdf",save.dir,top2mat[i,"fb"],
              gsub("Fungi_|Prokaryote_","",top2mat[i,"Taxa2"]),
              gsub("Fungi_|Prokaryote_","",top2mat[i,"Taxa2"]),
              top2mat[i,"plant"]),width = 9,height=3)
  
  sscomp
  
  dev.off() 
  #}
  
# #merge heatmap
# 
# #Fungi
rownames(msa_f) <- sf2[rownames(msa_f),4]
# 
jac_dist <- function(x) {
  vegdist(x, method = "jaccard",na.rm=TRUE)
}

row_ha = rowAnnotation(bar = rownames(msa_f),
                       show_legend = FALSE,
                       show_annotation_name = FALSE,
                       annotation_label="Basin ID ",
                       annotation_name_rot =45,
                       annotation_name_gp = gpar(size=3),
                       gp = gpar(col = "black"),
                       col=list(bar=colvec[rownames(msa_f)]))

col_fun = circlize::colorRamp2(c(0, 1), c("white", "darkorange"))

sscomp <- Heatmap(na.omit(msa_f)[,taxa_order$Fungi], row_names_side = "left", column_names_side = "bottom", 
                  col = col_fun,
                  #row_split=2,
                  show_row_names = FALSE,
                  left_annotation = row_ha,
                  show_heatmap_legend = FALSE,
                  row_dend_side = "left", rect_gp = gpar(col = "black"), 
                  row_names_gp = gpar(cex=0.6, fontface = "bold"), 
                  column_names_gp = gpar(cex=0.6, fontface = "bold"), 
                  row_dend_width = unit(2, "cm"),
                  column_dend_height = unit(2, "cm"), 
                  column_title = "Fungi ( Family )",
                  column_names_rot = 45, row_title = "Stable states", 
                  clustering_distance_rows = jac_dist,
                  clustering_method_rows = "ward.D2",
                  cluster_columns = FALSE,)





pdf(sprintf("%s/Heatmap_basins_Fungi.pdf",save.dir),width=4,height=8)
 sscomp
# 
 dev.off()
# 
# #Prok
rownames(msa_p) <- sb2[rownames(msa_p),4]

jac_dist <- function(x) {
  vegdist(x, method = "jaccard",na.rm=TRUE)
}

row_ha = rowAnnotation(bar = rownames(msa_p),
                       show_legend = FALSE,
                       show_annotation_name = FALSE,
                       annotation_label="Basin ID ",
                       annotation_name_rot =45,
                       annotation_name_gp = gpar(angle = 90),
                       gp = gpar(col = "black"),
                       col=list(bar=colvec[rownames(msa_p)]))

col_fun = colorRamp2(c(0, 1), c("white", "darkgreen"))

sscomp <- Heatmap(msa_p[,taxa_order$Prokaryote], row_names_side = "left", column_names_side = "bottom", 
                  col = col_fun,
                  #row_split = 10,
                  show_row_names = FALSE,
                  left_annotation = row_ha,
                  show_heatmap_legend = FALSE,
                  row_dend_side = "left", rect_gp = gpar(col = "black"), 
                  row_names_gp = gpar(cex=0.6, fontface = "bold"), 
                  column_names_gp = gpar(cex=0.6, fontface = "bold"), 
                  row_dend_width = unit(2, "cm"),
                  column_dend_height = unit(2, "cm"), 
                  column_title = "Prokaryotes ( Family )",
                  column_names_rot = 45, row_title = "Stable states", 
                  cluster_columns = FALSE,
                  clustering_distance_rows = jac_dist,
                  clustering_method_rows = "ward.D2")


pdf(sprintf("%s/Heatmap_basins_Prokaryote.pdf",save.dir),width=8,height=15)
sscomp
dev.off()

