









##########################################################################
set.seed(1234)
library(scales)
library(ggplot2)
library(ggtext)
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
save.dir <- "02_06_04_ELA_initial_state_numebr_260130"
dir.create(save.dir)


########################################################################
#read original functions
source("packages/01_1_function.R")
ocmat <- list(Fungi=readRDS(sprintf("%s/ELA_input_ocmat_Fungi.rds",dir_02_05)),
              Prokaryota=readRDS(sprintf("%s/ELA_input_ocmat_Prokaryote.rds",dir_02_05)))

sp_info <- readRDS(sprintf("%s/ELA_input_plant.rds",ELA_prep_dir))

fb <-"Fungi"

ocmatf <- ocmat[[fb]]

sa <- readRDS(sprintf("%s/%s/runSA_ST_full_%s_Family.rds",
                      dir_02_06,fb,fb))

pltab <- unique(sp_info[rownames(ocmatf),-c(1,ncol(sp_info))])
plnam <- c()

res <- NULL
rep <- c(100,1000,5000,10000,20000,40000,70000,100000,150000,200000)
#rep <- c(1000,100000,150000,300000)

for(pl in 1:nrow(pltab)){#pl <- 4
  cat("plant ",pl,"\n")
  if(sum(pltab[pl,]==1)){
    plnam[pl] <- colnames(pltab)[which(pltab[pl,]==1)]
  }else{
    plnam[pl] <- colnames(sp_info)[ncol(sp_info)]
  }
  
sa2p <- sa2params(sa,as.numeric(pltab[pl,]))

hgestp <- sa2p[[4]]
jestp <- sa2p[[2]]

uqm <- c()
for(i in 1:length(rep)){
  SS.itr <- rep[i]
    minsets <- as.data.frame(SSestimate(hgestp, jestp, itr = SS.itr))
  uqm[i] <- nrow(unique(minsets[, -length(minsets[1, ])]))
}

res <- rbind(res,data.frame(plant=plnam[pl],init_number=rep,basin_number=uqm))

}

saveRDS(res,sprintf("%s/ELA_initial_state_for_basin_%s.rds",
                    save.dir,fb))

res$plant <- sprintf("*%s*",res$plant)
g <- ggplot(res,aes(x=init_number,y=basin_number))+
  geom_point()+
  geom_line()+
  scale_y_continuous(breaks = scales::breaks_width(width = 1)) +
  xlab("Number of initial states")+
  ylab("Number of detected basins")+
  theme_bw()+
  theme(
    aspect.ratio=1,
    strip.text = element_markdown(size=14),
    axis.text=element_text(size=10),
    axis.title=element_text(size=14)
  )+
  facet_wrap(~plant,scales="free")

ggsave(sprintf("%s/ELA_initial_state_for_basin_%s.pdf",
                 save.dir,fb),plot=g,
       width=8,height=6)
################
fb <-"Prokaryota"

sa <- readRDS(sprintf("%s/%s/runSA_ST_full_%s_Family.rds",
                      dir_02_06,fb,fb))
pltab <- unique(sp_info[rownames(ocmatf),-c(1,ncol(sp_info))])
plnam <- c()


res <- NULL
rep <- c(100,1000,5000,10000,20000,40000,70000,100000,150000,200000)
#rep <- c(1000,100000,150000,300000)

for(pl in 1:nrow(pltab)){#pl <- 4
  cat("plant ",pl,"\n")
  if(sum(pltab[pl,]==1)){
    plnam[pl] <- colnames(pltab)[which(pltab[pl,]==1)]
  }else{
    plnam[pl] <- colnames(sp_info)[ncol(sp_info)]
  }
  
  sa2p <- sa2params(sa,as.numeric(pltab[pl,]))
  
  hgestp <- sa2p[[4]]
  jestp <- sa2p[[2]]
  
  uqm <- c()
  for(i in 1:length(rep)){
    SS.itr <- rep[i]
    minsets <- as.data.frame(SSestimate(hgestp, jestp, itr = SS.itr))
    uqm[i] <- nrow(unique(minsets[, -length(minsets[1, ])]))
  }
  
  res <- rbind(res,data.frame(plant=plnam[pl],init_number=rep,basin_number=uqm))
  
}

saveRDS(res,sprintf("%s/ELA_initial_state_for_basin_%s.rds",
                    save.dir,fb))

res$plant <- sprintf("*%s*",res$plant)
g <- ggplot(res,aes(x=init_number,y=basin_number))+
  geom_point()+
  geom_line()+
  scale_y_continuous(breaks = scales::breaks_width(width = 1)) +
  xlab("Number of initial states")+
  ylab("Number of detected basins")+
  theme_bw()+
  theme(
    aspect.ratio=1,
    strip.text = element_markdown(size=14),
    axis.text=element_text(size=10),
    axis.title=element_text(size=14)
  )+
  facet_wrap(~plant,scales="free")

ggsave(sprintf("%s/ELA_initial_state_for_basin_%s.pdf",
               save.dir,fb),plot=g,
       width=8,height=6)