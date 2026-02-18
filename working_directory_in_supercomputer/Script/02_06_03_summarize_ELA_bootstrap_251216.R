











##########################################################################
set.seed(1234)


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








#########################################################################
save.dir <- "02_06_03_summarize_ELA_bootstrap_251216"
dir.create(save.dir)


########################################################################
#read original functions
source("packages/01_1_function.R")

fb <- "Fungi"
boot_ela <- do.call(rbind,lapply(list.files(sprintf("%s/%s",dir_02_06_02,fb),
                                            pattern="energies",full.names = TRUE)[1:1000],readRDS))

be <- unique(boot_ela[,-ncol(boot_ela)])

be$lCI <- NA
be$uCI <- NA

for(i in 1:nrow(be)){
   be[i,c("lCI","uCI")] <- quantile(boot_ela[boot_ela[,"plant"]==be[i,"plant"] &
                                         boot_ela[,"state_id"]==be[i,"state_id"],
                    "energy"],c(0.025,0.975))
}

saveRDS(be,
        sprintf("%s/energy_CI_bootstrap_%s.rds",save.dir,fb))

fb <- "Prokaryota"
boot_ela <- do.call(rbind,lapply(list.files(sprintf("%s/%s",dir_02_06_02,fb),
                                            pattern="energies",full.names = TRUE)[1:1000],readRDS))

be <- unique(boot_ela[,-ncol(boot_ela)])

be$lCI <- NA
be$uCI <- NA

for(i in 1:nrow(be)){
  be[i,c("lCI","uCI")] <- quantile(boot_ela[boot_ela[,"plant"]==be[i,"plant"] &
                                              boot_ela[,"state_id"]==be[i,"state_id"],
                                            "energy"],c(0.025,0.975))
}

saveRDS(be,
        sprintf("%s/energy_CI_bootstrap_%s.rds",save.dir,fb))
