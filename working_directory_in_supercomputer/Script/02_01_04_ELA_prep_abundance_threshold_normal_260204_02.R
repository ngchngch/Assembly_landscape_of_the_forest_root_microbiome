











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

save.dir <- "02_01_04_ELA_prep_abundance_threshold_normal_260204_02"
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

inputdir <- "threadlipper_analysis_series/240430_caverage_rarefaction"
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


#############
##1270 samples in tb_all
tb_all <- list(Fungi=tb_al$Fungi[which(rownames(tb_al$Fungi) %in% rownames(sp_info)),],
               Prokaryota=tb_al$Prokaryota[which(rownames(tb_al$Prokaryota) %in% rownames(sp_info)),])
tb_m <- list(Fungi=tb_m1$Fungi[which(rownames(tb_m1$Fungi) %in% rownames(sp_info)),],
             Prokaryota=tb_m1$Prokaryota[which(rownames(tb_m1$Prokaryota) %in% rownames(sp_info)),])
saveRDS(tb_m,sprintf("%s/Df_RA_%s.rds",save.dir,taxa))
##########################################

###########################################

fb <- "Fungi"

relf <- tb_m$Fungi/rowSums(tb_m$Fungi)

relf2 <- relf[,which(colnames(relf) != "Unidentified" & colSums(relf>0)/nrow(relf) > 0.03)]


for(rel_th in seq(0,0.1,by=0.01)){
ocmat <- ifelse(relf2<=rel_th,0,1)
saveRDS(ocmat,sprintf("%s/ocmat_remove_M2SD_%s_%s_relth%s.rds",save.dir,fb,taxa,rel_th))
}

saveRDS(relf2,sprintf("%s/ramat_for_%s_%s.rds",save.dir,fb,taxa))
###########################

fb <- "Prokaryota"

relf <- tb_m[[fb]]/rowSums(tb_m[[fb]])

relf2 <- relf[,which(colnames(relf) != "Unidentified" & colSums(relf>0) > 30)]

for(rel_th in seq(0,0.1,by=0.01)){
ocmat <- ifelse(relf2<=rel_th,0,1)
saveRDS(ocmat,sprintf("%s/ocmat_remove_M2SD_%s_%s_relth%s.rds",save.dir,fb,taxa,rel_th))
}
saveRDS(relf2,sprintf("%s/ramat_for_%s_%s.rds",save.dir,fb,taxa))