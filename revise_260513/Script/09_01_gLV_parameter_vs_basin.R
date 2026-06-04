setwd("/Users/ngch/Desktop/Now_using/data/sugaFB_ELA/revise_260513")
source("packages/fun_gLV_simulation_revise_260514.R")

# %% loading library and functions
library(tidyverse)
library(ggplot2)
library(doParallel)
library(igraph)
library(deSolve)
library(vegan)
library(cluster)
library(rELA)

save.dir <- "Analysis/09_01_gLV_parameter_vs_basin"
dir.create(save.dir)

dir_05_01 <- "Basedata/05_01_parameters_and_matrices" 

Alist <- list.files(dir_05_01,pattern="Amatrix",full.names=TRUE,recursive=TRUE)

comm_list <- list.files(dir_05_01,pattern="community_matrix",full.names=TRUE,recursive=TRUE)

n.core <- 8
#ELA parameters
qth <- 10^-5 # do not change!!
SS.itr <- 20000

seeds <- sapply(strsplit(Alist,"/"),"[",3)
conn <- sapply(strsplit(Alist,"_"),function(x){as.numeric(gsub("seed|.rds","",x[[length(x)]]))})

cluster_summary <- NULL

for(n in 1:length(Alist)){
  
  cat("\r", n)
  
  A <- readRDS(Alist[n])[[1]]
  
  comm <- readRDS(
    comm_list[
      grepl(paste0(seeds[n],"/"), comm_list) &
        grepl(conn[n], comm_list)
    ]
  )
  
  bin0 <- as.matrix(comm$bin_mat);dimnames(bin0) <- list(paste0("S",1:nrow(bin0)),paste0("Sp",1:ncol(bin0)))
  
  bin <- bin0[
    ,
    colMeans(bin0) > 0.1 &
      colMeans(bin0) < 0.9
  ]
  
  sa <- runSA(
    ocmat  = bin,
    qth    = qth,
    rep    = 1280,
    threads= n.core
  )
  
  saveRDS(
    sa,
    file = sprintf(
      "%s/runSA_%s_%s.rds",
      save.dir,
      seeds[n],
      conn[n]
    )
  )
  
  para <- sa2params(sa)
  
  states <- SSestimate(
    para[[1]],
    para[[2]],
    itr = 150000
  )
  
  n_basin <- nrow(unique(states[,-ncol(states)]))
  
  # =========================================================
  # network properties
  # =========================================================
  
  ## connectance
  connectance <- conn[n]
  
  ## mean interaction strength
  ## （非ゼロ要素のみ）
  nonzero_A <- A[A != 0]
  
  mean_interaction_strength <- mean(abs(nonzero_A))
  
  
  # =========================================================
  # summary table
  # =========================================================
  
  cluster_summary <- rbind(
    cluster_summary,
    data.frame(
      seed = seeds[n],
      connectance = connectance,
      mean_interaction_strength =
        mean_interaction_strength,
      n_basin = n_basin
    )
  )
  
  write.csv(
    cluster_summary,
    sprintf("%s/cluster_summary.csv", save.dir),
    row.names = FALSE
  )
}
