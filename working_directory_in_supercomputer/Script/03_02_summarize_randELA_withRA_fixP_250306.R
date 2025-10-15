##00_install_packages
##01_01_covarage_based_rarefaction_all
##01_02_Barplot 
##02_01_ELA_threshold
##02_02_ELA 
##02_03_graphics_ELA
##03_01_covarage_based_rarefaction_each
##03_02_ELA_withRAgradient #now
##03_03_graphics_ELA_withRAgradient
##04_01_SSclustering
##04_02_Change_eachSScluster_withRA
##04_03_graphics_Change_eachSScluster
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
save.dir <- "03_02_summarize_randELA_withRA_fixP"
dir.create(save.dir)

########################################################################
#read original functions
source("packages/01_1_function.R")

rand_dir <- c(list.dirs(path="03_02_randELA_withRA_4s_fixPS_1_3000_250227", recursive = FALSE),
              list.dirs(path="03_02_randELA_withRA_4s_fixPS_3001_6000_250227", recursive = FALSE),
              list.dirs(path="03_02_randELA_withRA_4s_fixPS_9001_10500_250227", recursive = FALSE),
              list.dirs(path="03_02_randELA_withRA_4s_fixPS_6001_9000_250227", recursive = FALSE))

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
  
  rand_lists1 <- list.dirs(rand_dir[grep(fb,rand_dir)],
                           recursive = FALSE)
  
  rand_lists <- list.files(rand_lists1[grepl(tx,rand_lists1,fixed = TRUE)],
                           recursive = TRUE,full.names = TRUE)
  
  rf <- mclapply(rand_lists,
               readRDS,mc.cores=n.core)
  
  dle <- sapply(rf,function(x){
    c(x[,"d_land"],x[,"d_even"])
  })
  
  dland <- dle[c(1:18),]
  deven <- dle[c(19:36),]
    
  saveRDS(dland,sprintf("%s/%s/%s/dland_1_3000_6001_10500.rds",save.dir,fb,tx))
  saveRDS(deven,sprintf("%s/%s/%s/deven_1_3000_6001_10500.rds",save.dir,fb,tx))
  
  fb <- "Prokaryota"
  dir.create(sprintf("%s/%s",save.dir,fb))
  
  dir.create(sprintf("%s/%s/%s",save.dir,fb,tx))
  
  rand_lists1 <- list.dirs(rand_dir[grep(fb,rand_dir)],
                           recursive = FALSE)
  
  rand_lists <- list.files(rand_lists1[grepl(tx,rand_lists1,fixed = TRUE)],
                           recursive = TRUE,full.names = TRUE)
  
  rf <- mclapply(rand_lists,
                 readRDS,mc.cores=n.core)
  
  dle <- sapply(rf,function(x){
    c(x[,"d_land"],x[,"d_even"])
  })
  
  dland <- dle[c(1:18),]
  deven <- dle[c(19:36),]
  
  saveRDS(dland,sprintf("%s/%s/%s/dland_1_10500.rds",save.dir,fb,tx))
  saveRDS(deven,sprintf("%s/%s/%s/deven_1_10500.rds",save.dir,fb,tx))
  