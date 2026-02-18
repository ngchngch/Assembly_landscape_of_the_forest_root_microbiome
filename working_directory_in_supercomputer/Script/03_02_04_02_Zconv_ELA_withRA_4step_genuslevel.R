
##########################################################################
#install.packages("/Volumes/8TBHDD_NGCH/sugadaira_bacteria_2023/240801_SSchange_randamize/packages/rELA.v0.51.tar")

#setwd("/Users/ngch/Desktop/Now_using/data/sugaFB_ELA/analysis_series_SpCom/02_ELA_241002")
set.seed(1234)

library(parallel)
library(foreach)
library(vegan)
library("Rcpp")
library("RcppArmadillo")
library('gtools')
library('RColorBrewer')
#library('scatterpie')
library("rELA")
#library('tidygraph')
















#########################################################################
save.dir <- "03_02_04_02_Zconv_ELA_withRA_4step_genuslevel"
dir.create(save.dir)

########################################################################
#read original functions
source("packages/01_1_function.R")

list_or <- list.files(path=dir_03_01_03,"ELA_SSprob_diff", recursive = TRUE,full.names=TRUE)
rand_dir <- list.dirs(path=dir_03_02_04, recursive = FALSE)

tagmat <- readRDS(sprintf("%s/ExpVar_clrRA_target_taxa.rds",ELA_prep_dir))

tagtaxs <- setdiff(unique(c(colnames(tagmat$Fungi),colnames(tagmat$Prokaryota))),"Unidentified")

#t2 <- unique(sapply(strsplit(orl,"/"),function(x)x[11]))
#tagtaxs <- intersect(tagtaxs,t2)
#5;29;35;40;44;109;166
#tagtaxs <- c("Amanita","Pseudonocardia")
#setdiff(tagtaxs,zval_l[,"target"])

  tx  <- tagtaxs[sp]
  
  fb <- "Fungi"
  dir.create(sprintf("%s/%s",save.dir,fb))
  
  dir.create(sprintf("%s/%s/%s",save.dir,fb,tx))
  
  or_dle <- readRDS(list_or[grepl(paste(fb,tx,sep="/"),list_or,fixed = TRUE)])
  
  rand_lists1 <- list.dirs(rand_dir[grep(fb,rand_dir)],
                           recursive = FALSE)
  
  rand_lists <- list.files(rand_lists1[grepl(tx,rand_lists1,fixed = TRUE)],
                           recursive = TRUE,full.names = TRUE)
  
  rf <- mclapply(rand_lists,
               readRDS,mc.cores=n.core)
  
  dle <- sapply(rf,function(x){
    c(x[,"d_land"],x[,"d_even"])
  })
  
  dland <- dle[c(1:nrow(or_dle)),]
  deven <- dle[c(nrow(or_dle)+1:nrow(or_dle)),]
    
  saveRDS(dland,sprintf("%s/%s/%s/dland_1_3000.rds",save.dir,fb,tx))
  saveRDS(deven,sprintf("%s/%s/%s/deven_1_3000.rds",save.dir,fb,tx))
  
  
  z_dland <- (or_dle[,"d_land"]-rowMeans(dland))/apply(dland,1,sd)
  z_deven <- (or_dle[,"d_even"]-rowMeans(deven))/apply(deven,1,sd)
  
    p_dland <- rowSums(or_dle[,"d_land"] <= dland)/ncol(dland)
p_deven <- ifelse(or_dle[,"d_even"]>0,
                  rowSums(or_dle[,"d_even"] <= deven)/ncol(deven),
                  rowSums(or_dle[,"d_even"] >= deven)/ncol(deven))
  
  
  res_df <- data.frame(ID=sp,fb=fb,taxa=tx,
                                    plant=or_dle[,"plant"],
                                    ra=or_dle[,"ra"],
                                    ab=or_dle[,"ab"],
                                    scale_ab=or_dle[,"std_ab_clr"],
                                    d_land=or_dle[,"d_land"],
                                    z_dland=z_dland,
                                    p_dland=p_dland,
                                    d_even=or_dle[,"d_even"],
                                    z_deven=z_deven,
                                    p_deven=p_deven)
  saveRDS(res_df,
          sprintf("%s/%s/%s/Zconv_ELA_results.rds",save.dir,fb,tx))
  ##########
  fb <- "Prokaryota"
  dir.create(sprintf("%s/%s",save.dir,fb))
  
  dir.create(sprintf("%s/%s/%s",save.dir,fb,tx))
  
  or_dle <- readRDS(list_or[grepl(paste(fb,tx,sep="/"),list_or,fixed = TRUE)])
  
  rand_lists1 <- list.dirs(rand_dir[grep(fb,rand_dir)],
                           recursive = FALSE)
  
  rand_lists <- list.files(rand_lists1[grepl(tx,rand_lists1,fixed = TRUE)],
                           recursive = TRUE,full.names = TRUE)
  
  rf <- mclapply(rand_lists,
                 readRDS,mc.cores=n.core)
  
  dle <- sapply(rf,function(x){
    c(x[,"d_land"],x[,"d_even"])
  })
  
  
  dland <- dle[c(1:nrow(or_dle)),]
  deven <- dle[c(nrow(or_dle)+1:nrow(or_dle)),]
    
  saveRDS(dland,sprintf("%s/%s/%s/dland_1_3000.rds",save.dir,fb,tx))
  saveRDS(deven,sprintf("%s/%s/%s/deven_1_3000.rds",save.dir,fb,tx))
  
  
  z_dland <- (or_dle[,"d_land"]-rowMeans(dland))/apply(dland,1,sd)
  z_deven <- (or_dle[,"d_even"]-rowMeans(deven))/apply(deven,1,sd)
  
  p_dland <- rowSums(or_dle[,"d_land"] <= dland)/ncol(dland)
p_deven <- ifelse(or_dle[,"d_even"]>0,
                  rowSums(or_dle[,"d_even"] <= deven)/ncol(deven),
                  rowSums(or_dle[,"d_even"] >= deven)/ncol(deven))
  
  res_df <- data.frame(ID=sp,fb=fb,taxa=tx,
                                    plant=or_dle[,"plant"],
                       ra=or_dle[,"ra"],
                       ab=or_dle[,"ab"],
                       scale_ab=or_dle[,"std_ab_clr"],
                       d_land=or_dle[,"d_land"],
                       z_dland=z_dland,
                       p_dland=p_dland,
                       d_even=or_dle[,"d_even"],
                       z_deven=z_deven,
                       p_deven=p_deven)
  saveRDS(res_df,
          sprintf("%s/%s/%s/Zconv_ELA_results.rds",save.dir,fb,tx))