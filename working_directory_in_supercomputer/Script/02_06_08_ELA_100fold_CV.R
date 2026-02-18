
##########################################################################
set.seed(1234)

library(ggplot2)
library(ggstar)
library(parallel)
library(foreach)
library(vegan)
library("Rcpp")
library("RcppArmadillo")
library("stringdist")
library("doParallel")
library('tidyverse')
library('gtools')
library('igraph')
library('RColorBrewer')
library("stringdist")
#library('scatterpie')
library("rELA")
library("ggtext")


library("ggforce")

#install.packages("/Users/ngch/Desktop/Now_using/data/sugaFB_ELA/packages/rELA.v0.51.tar.gz")





#########################################################################
save.dir <- "02_06_08_ELA_100fold_CV"
dir.create(save.dir)


########################################################################
#read original functions
source("packages/01_1_function.R")

########################################################################
tx_f <- readRDS("Base_data/Fungi/taxa_list_mod.rds")
tx_b <- as.data.frame(readRDS("Base_data/OTU_basedata_set/taxonomy_list_dada2.rds"))
####################################################
taxa <- "Family"





#read data
rel_df <- readRDS(sprintf("%s/Df_RA_%s.rds",ELA_prep_dir,taxa))
ocmat <- list(Fungi=readRDS(sprintf("%s/ELA_input_ocmat_Fungi.rds",dir_02_05)),
              Prokaryota=readRDS(sprintf("%s/ELA_input_ocmat_Prokaryote.rds",dir_02_05)))

info <- readRDS(sprintf("%s/comp_sample_info_plant2.rds",ELA_prep_dir))
sp_info <- readRDS(sprintf("%s/ELA_input_plant.rds",ELA_prep_dir))
SA_full <- list(Fungi=readRDS(sprintf("%s/Fungi/runSA_ST_full_Fungi_%s.rds",dir_02_06,taxa)),
                Prokaryota=readRDS(sprintf("%s/Prokaryota/runSA_ST_full_Prokaryota_%s.rds",dir_02_06,taxa)))
###################################################
#ELA parameters
qth <- 10^-5 # do not change!!
SS.itr <- 150000
fb <-"Fungi"
dir <- sprintf("%s/%s",save.dir,fb)
dir.create(dir)

  ocmatf0 <- ocmat[[fb]]

  if(nfold==1){
    fold_id <- sample(rep(1:100,length.out=nrow(ocmatf0)))
    saveRDS(file=sprintf("%s/fold_id_%s.rds",dir,fb),fold_id)
  }else{
    fold_id <- readRDS(sprintf("%s/fold_id_%s.rds",dir,fb))
  }
  
  ocmatf <- ocmatf0[fold_id!=nfold,]
  ans_ocmatf <- ocmatf0[fold_id==nfold,]
  
  sa_full <- SA_full[[fb]]
  sa <- runSA(ocmat=as.matrix(ocmatf),
              enmat = sp_info[rownames(ocmatf),-c(1,ncol(sp_info))],
              qth=qth, rep=1280, threads=n.core)
  
  saveRDS(sa, file=sprintf("%s/runSA_ST_full_%s_%s_fold%s.rds",dir,fb,taxa,nfold))
  
  pltab <- sp_info[rownames(ans_ocmatf),-c(1,ncol(sp_info))]
  
  cv_res <- NULL
  for(pl in 1:nrow(pltab)){#pl <- 4
    if(sum(pltab[pl,]==1)){
      plnam <- colnames(pltab)[which(pltab[pl,]==1)]
    }else{
      plnam <- colnames(sp_info)[ncol(sp_info)]
    }
    
    sa2p <- sa2params(sa,as.numeric(pltab[pl,]))
    sa2p_f <- sa2params(sa_full,as.numeric(pltab[pl,]))
    
    hgestp <- sa2p[[4]]
    jestp <- sa2p[[2]]
    hgestp_f <- sa2p_f[[4]]
    jestp_f <- sa2p_f[[2]]
    
    b_cv <- Bi(ans_ocmatf[pl,],hgestp,jestp)[[1]]
    b_f <- Bi(ans_ocmatf[pl,],hgestp_f,jestp_f)[[1]]
    
    d <- vegdist(rbind(id2bin(b_cv,ncol(ans_ocmatf)),
                  id2bin(b_f,ncol(ans_ocmatf))),method="jaccard")
    
    cv_res <- rbind(cv_res,data.frame(plant=plnam,
               fold_id=nfold,
               Basin_CV=b_cv,
               Basin_full=b_f,
               dist_jaccard=as.numeric(d[1])))
  }
    
  saveRDS(cv_res,file=sprintf("%s/CV_results_%s_%s_fold%s.rds",dir,fb,taxa,nfold))
  
  
#######################################################################
fb <-"Prokaryota"
  dir <- sprintf("%s/%s",save.dir,fb)
  dir.create(dir)
  
  ocmatf0 <- ocmat[[fb]]
  
  if(nfold==1){
    fold_id <- sample(rep(1:100,length.out=nrow(ocmatf0)))
    saveRDS(file=sprintf("%s/fold_id_%s.rds",dir,fb),fold_id)
  }else{
    fold_id <- readRDS(sprintf("%s/fold_id_%s.rds",dir,fb))
  }
  
  ocmatf <- ocmatf0[fold_id!=nfold,]
  ans_ocmatf <- ocmatf0[fold_id==nfold,]
  
  sa_full <- SA_full[[fb]]
  sa <- runSA(ocmat=as.matrix(ocmatf),
              enmat = sp_info[rownames(ocmatf),-c(1,ncol(sp_info))],
              qth=qth, rep=1280, threads=n.core)
  
  saveRDS(sa, file=sprintf("%s/runSA_ST_full_%s_%s_fold%s.rds",dir,fb,taxa,nfold))
  
  pltab <- sp_info[rownames(ans_ocmatf),-c(1,ncol(sp_info))]
  
  cv_res <- NULL
  for(pl in 1:nrow(pltab)){#pl <- 4
    if(sum(pltab[pl,]==1)){
      plnam <- colnames(pltab)[which(pltab[pl,]==1)]
    }else{
      plnam <- colnames(sp_info)[ncol(sp_info)]
    }
    
    sa2p <- sa2params(sa,as.numeric(pltab[pl,]))
    sa2p_f <- sa2params(sa_full,as.numeric(pltab[pl,]))
    
    hgestp <- sa2p[[4]]
    jestp <- sa2p[[2]]
    hgestp_f <- sa2p_f[[4]]
    jestp_f <- sa2p_f[[2]]
    
    b_cv <- Bi(ans_ocmatf[pl,],hgestp,jestp)[[1]]
    b_f <- Bi(ans_ocmatf[pl,],hgestp_f,jestp_f)[[1]]
    
    d <- vegdist(rbind(id2bin(b_cv,ncol(ans_ocmatf)),
                       id2bin(b_f,ncol(ans_ocmatf))),method="jaccard")
    
    cv_res <- rbind(cv_res,data.frame(plant=plnam,
                                      fold_id=nfold,
                                      Basin_CV=b_cv,
                                      Basin_full=b_f,
                                      dist_jaccard=as.numeric(d[1])))
  }
  
  saveRDS(cv_res,file=sprintf("%s/CV_results_%s_%s_fold%s.rds",dir,fb,taxa,nfold))
  