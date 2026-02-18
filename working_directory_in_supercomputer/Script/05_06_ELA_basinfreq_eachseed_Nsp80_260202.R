
##########################################################################
#install.packages("/Users/ngch/Desktop/Now_using/data/sugaFB_ELA/packages/rELAv0.51_fujita240702.tar.gz")
#setwd("/Users/ngch/Desktop/Now_using/data/sugaFB_ELA")

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





#########################################################################
save.dir <- "05_06_ELA_basinfreq_eachseed_Nsp80_260202"
dir.create(save.dir)

########################################################################
#read original functions
source("packages/01_1_function.R")

#from fujita
Rcpp::sourceCpp("packages/ELA_functions_v060.cpp")

########################################################################
sc <- readRDS(sprintf("%s/top3_seed_across.rds",dir_05_01_00))
seed <- sc[r,2]
conn0 <- sc[r,1]

set.seed(seed)

set.seed(seed)
cpath0 <- list.files(sprintf("%s/seed%s",dir_05_01,seed),pattern = "community_matrix",full.names=TRUE)
cpath0_wi <- list.files(sprintf("%s/seed%s",dir_05_01,seed),pattern = "dyn_wi",full.names=TRUE)
cpath0_woi <- list.files(sprintf("%s/seed%s",dir_05_01,seed),pattern = "dyn_woi",full.names=TRUE)
cpath <- cpath0[grep(paste0("conn",conn0),cpath0)]
cpath_wi <- cpath0_wi[grep(paste0("conn",conn0),cpath0_wi)]
cpath_woi <- cpath0_woi[grep(paste0("conn",conn0),cpath0_woi)]

cnam <- sapply(strsplit(cpath,"_"),function(x){gsub(".rds","",x[length(x)])})
dir.create(sprintf("%s/seed%s",save.dir,seed))

#ELA parameters
qth <- 10^-5 # do not change!!
SS.itr <- 150000
even <- NULL

for(conn in 1:length(cpath)){
  comm_list <- readRDS(cpath[conn])
  comm_list_wi <- readRDS(cpath_wi[conn])
  comm_list_woi <- readRDS(cpath_woi[conn])
  cnam1 <- cnam[conn]
  dir.create(sprintf("%s/seed%s/%s",save.dir,seed,cnam1))
  
  bin_th <- comm_list$bin_th
  mat_wi0 <- t(sapply(comm_list_wi,function(x){
  ifelse(x[[sp]][sample(1:80,1),-sp]>bin_th[1],1,0)}))
  mat_woi0 <- t(sapply(comm_list_woi,function(x){
    ifelse(x[[sp]][sample(1:80,1),-sp]>bin_th[1],1,0)}))
  
  dimnames(mat_wi0) <- list(paste0("S",1:nrow(mat_wi0)),paste0("Sp",1:ncol(mat_wi0)))
  dimnames(mat_woi0) <- list(paste0("S",1:nrow(mat_woi0)),paste0("Sp",1:ncol(mat_woi0)))
  #Rcpp::sourceCpp("for_Github/packages/ELA_functions_v060.cpp")
  #sp <- 1
  
  mat_wi <- mat_wi0[,colMeans(mat_wi0)>0.1&colMeans(mat_wi0)<0.9]
  mat_woi <- mat_woi0[,colnames(mat_wi)]
  
  ###################################################
  
  part <- sprintf("%s/seed%s/%s/%s",save.dir,seed,cnam1,paste0("ID",sp))
  dir.create(part)
  
  cat(sprintf("%s\n",sp))
  
  #list[ocmatf, abmatf, enmatf, samplelabelf, specieslabelf, factorlabelf] <-
  
  sa_wi <- runSA(ocmat=as.matrix(mat_wi),
              qth=qth, rep=1280, threads=n.core)
  saveRDS(sa_wi, file=sprintf("%s/runSA_wi_wo%s_%s.rds",part,sp,cnam1))
  
  sa_woi <- runSA(ocmat=as.matrix(mat_woi),
                 qth=qth, rep=1280, threads=n.core)
  saveRDS(sa_woi, file=sprintf("%s/runSA_woi_wo%s_%s.rds",part,sp,cnam1))
  
  #sa <- readRDS("/Users/ngch/Desktop/Now_using/data/sugaFB_ELA/analysis_series/03_03_same_samp_ELA_withRA_Genus_240703/Fungi/Acidibacter/runSA_Fungi_woAcidibacter.rds")
  #make input start sets
  hg <- sa2params(sa_wi)[[4]]
  state <- foreach(i=1:SS.itr,.combine="rbind")%do%{
    st <- runif(length(hg), 0, 2) |> as.integer()
  }
  rownames(state) <- sprintf("Start_%05d",1:SS.itr)
  saveRDS(state, file=sprintf("%s/start_%s.rds",part,cnam1))
  
  param_wi <- sa2params(sa_wi)
  param_woi <- sa2params(sa_woi)
  
  hge_wi <- param_wi[[4]]
  je_wi <- param_wi[[2]]
  hge_woi <- param_woi[[4]]
  je_woi <- param_woi[[2]]
  
  ss_wi <- SSestimate_given(hge_wi, je_wi, state)
  ss_woi <- SSestimate_given(hge_woi, je_woi, state)
  
  mst_wi <- ss_wi[,-ncol(ss_wi)]
  mst_woi <- ss_woi[,-ncol(ss_woi)]

  id_wi <- apply(mst_wi, 1, paste, collapse='')
  id_woi <- apply(mst_woi, 1, paste, collapse='')
  
  prop_wi <- as.numeric(table(id_wi))/SS.itr
  prop_woi <- as.numeric(table(id_woi))/SS.itr
  
  even <- rbind(even,data.frame(id=sp,
                                seed=seed,
             conn=cnam1,
             nbasin_wi=length(table(id_wi)),
             nbasin_woi=length(table(id_woi)),
             d_nbasin=length(table(id_wi))-length(table(id_woi)),
             shannon_wi=-sum(prop_wi*log(prop_wi)),
             shannon_woi=-sum(prop_woi*log(prop_woi)),
             d_shannon=-sum(prop_wi*log(prop_wi))+sum(prop_woi*log(prop_woi)))
             )
  
  cat("|\n")  
}

saveRDS(even, file=sprintf("%s/seed%s/ELA_evenness_wiwoi_sp%s_conn%s.rds",save.dir,seed,sp,conn0))