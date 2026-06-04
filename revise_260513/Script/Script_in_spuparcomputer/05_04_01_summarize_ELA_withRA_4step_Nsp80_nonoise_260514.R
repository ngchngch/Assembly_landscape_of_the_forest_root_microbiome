











##########################################################################
#install.packages("/Users/ngch/Desktop/Now_using/data/sugaFB_ELA/packages/rELAv0.51_fujita240702.tar.gz")
#setwd("/Users/ngch/Desktop/Now_using/data/sugaFB_ELA")


library(ggplot2)
library(ggstar)
library(parallel)
library(foreach)
library(vegan)
library("Rcpp")
library("RcppArmadillo")
library("doParallel")
library('tidyverse')
library('gtools')
library('igraph')
library('RColorBrewer')
library("stringdist")
#library('scatterpie')
library("rELA")
library("ggstar")



#ELA_prep_dir='02_01_ELA_prep_abundance_threshold'
#n.core=8
#dir_02_03='02_03_summarize_occ_taxa_th'




#dir_05_01 <- "05_01_gLV_simulation"
#dir_05_02 <- "05_02_ELA_withRA_multistep"
#########################################################################
save.dir <- "05_04_01_summarize_ELA_withRA_4step_Nsp80_nonoise_260514"
dir.create(save.dir)

########################################################################
#read original functions
source("packages/01_1_function.R")

set.seed(12345)
dl <- NULL

sc <- readRDS(sprintf("%s/top3_seed_across.rds",dir_05_01_00))

for(r in 1:15){
seed <- sc[r,2]
conn <- sc[r,1]

  connect <- paste0("conn",conn)
    
    delta_list <- list.files(sprintf("%s/seed%s/%s",dir_05_02,seed,connect),
                             pattern="ELA_SSprob",
                             recursive=TRUE,full.names = TRUE)
    kness <- readRDS(sprintf("%s/seed%s/keystonness_vs_index_%s.rds",dir_05_01,seed,connect))
    
    k_summary <- readRDS(sprintf("%s/seed%s/keystonness_summary_%s.rds",dir_05_01,seed,connect))
    
    key_parm <- unique(kness[,-c(3,8:9)])
    
    dd <- do.call(rbind,lapply(delta_list,function(x){x2 <- readRDS(x)
    x3<-x2[,c(1,(ncol(x2)-2:0))];colnames(x3) <- c("id","ra","d_land","d_even")
    return(x3)
    }))
    
    head(dd)
    dl <- rbind(dl,cbind(seed=seed,
                         connectance=as.numeric(gsub("conn","",connect)),
                         dd,
                         k_summary[match(dd$id,1:80),],
                         key_parm[match(dd$id,key_parm[,"id"]),]
    ))
}

saveRDS(dl,sprintf("%s/ELA_withRA_summary_allseed.rds",save.dir))