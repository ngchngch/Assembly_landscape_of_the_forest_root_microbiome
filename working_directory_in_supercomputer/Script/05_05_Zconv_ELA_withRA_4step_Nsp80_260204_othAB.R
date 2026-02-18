
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
save.dir <- "05_05_Zconv_ELA_withRA_4step_Nsp80_260204_othAB"
dir.create(save.dir)

########################################################################
#read original functions
source("packages/01_1_function.R")

p_seq <- c("perc15","perc35","perc65","perc85")

res_df <- NULL

sc <- readRDS(sprintf("%s/top3_seed_across.rds",dir_05_01_00))

for(r in 1:15){
seed <- sc[r,2]
connect <- sc[r,1]
    
    set.seed(seed)
    dir.create(sprintf("%s/seed%s",save.dir,seed))
    dir.create(sprintf("%s/seed%s/conn%s",save.dir,seed,connect))
    
      or_dir <- list.dirs(path=paste0(dir_05_02,"/seed",seed,"/conn",connect), recursive = FALSE)
      rand_dir <- list.dirs(path=paste0(dir_05_03,"/seed",seed,"/conn",connect), recursive = FALSE)
      
      for(sp in 1:80){
        cat("Processing seed:",seed," connect:",connect," Sp:",sp,"\n")
      tx  <- paste0("ID",sp)
      
      dir.create(sprintf("%s/seed%s/conn%s/%s",save.dir,seed,connect,tx))
      
      rand_lists <- list.files(rand_dir[grep(paste0(tx,"$"),rand_dir)],pattern="ELA_SSprob_diff",
                               recursive = TRUE,full.names = TRUE)
      
      or_dle <- readRDS(list.files(pattern="ELA_SSprob_diff",
                                   or_dir[grep(paste0(tx,"$"),or_dir)],
                                   recursive = TRUE,full.names = TRUE))
      
      if(any(!is.na(or_dle[,"d_land"]))){
        rf <- mclapply(rand_lists,
                       readRDS,mc.cores=n.core)
        
        dle <- sapply(rf,function(x){
          c(x[,"d_land"],x[,"d_even"])
        })
        
        dland <- dle[c(1:4),]
        deven <- dle[c(5:8),]
        
        saveRDS(dland,sprintf("%s/seed%s/conn%s/%s/dland_1_3000.rds",save.dir,seed,connect,tx))
        saveRDS(deven,sprintf("%s/seed%s/conn%s/%s/deven_1_3000.rds",save.dir,seed,connect,tx))
        
        z_dland <- (or_dle[match(p_seq,or_dle[,"ra"]),"d_land"]-rowMeans(dland))/apply(dland,1,sd)
        z_deven <- (or_dle[match(p_seq,or_dle[,"ra"]),"d_even"]-rowMeans(deven))/apply(deven,1,sd)
        
        p_dland <- c()
        p_deven <- c()
        for(perc in 1:length(p_seq)){
          p_dland[perc] <- sum(or_dle[match(c(p_seq[perc]),or_dle[,"ra"]),"d_land"] <= dland[perc,])/ncol(dland)
          if(or_dle[match(c(p_seq[perc]),or_dle[,"ra"]),"d_even"]>0){
            p_deven[perc] <- sum(or_dle[match(c(p_seq[perc]),or_dle[,"ra"]),"d_even"] <= deven[perc,])/ncol(deven)
          }else{
            p_deven[perc] <- sum(or_dle[match(c(p_seq[perc]),or_dle[,"ra"]),"d_even"] >= deven[perc,])/ncol(deven)
          }
        }
        
        res_df <- rbind(res_df,data.frame(ID=sp,
                                          seed=seed,
                                          connectance=connect,
                                          ra=p_seq,
                                          ab=or_dle[match(p_seq,or_dle[,"ra"]),"ab"],
                                          scale_ab=or_dle[match(p_seq,or_dle[,"ra"]),"scale_ab"],
                                          d_land=or_dle[match(p_seq,or_dle[,"ra"]),"d_land"],
                                          z_dland=z_dland,
                                          p_dland=p_dland,
                                          d_even=or_dle[match(p_seq,or_dle[,"ra"]),"d_even"],
                                          z_deven=z_deven,
                                          p_deven=p_deven))
        
      }else{
        res_df <- rbind(res_df,data.frame(ID=sp,
                                          seed=seed,
                                          connectance=connect,
                                          ra=p_seq,
                                          ab=rep(NA,length(p_seq)),
                                          scale_ab=rep(NA,length(p_seq)),
                                          d_land=rep(NA,length(p_seq)),
                                          z_dland=rep(NA,length(p_seq)),
                                          p_dland=rep(NA,length(p_seq)),
                                          d_even=rep(NA,length(p_seq)),
                                          z_deven=rep(NA,length(p_seq)),
                                          p_deven=rep(NA,length(p_seq))))
      }
      
      
    }
}

saveRDS(res_df,sprintf("%s/summarized_ELA_withRA_4step_Nsp80.rds",save.dir))
