











##########################################################################
set.seed(1234)

library(ggplot2)
library(ggstar)
library(parallel)
library(foreach)
library(vegan)
library("Rcpp")
library("RcppArmadillo")
library("patchwork")
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







#########################################################################
#setwd("/Users/ngch/Desktop/Now_using/data/sugaFB_ELA/revise_251212/output_spacom/result_for_revise/analysis_in_local_computer/Output")
save.dir <- "Output/02_06_01_02_PCoA_basins_samples"
dir.create(save.dir)

dir_02_01 <- "Output_supercomputer/02_01_ELA_prep_abundance_threshold" 
dir_02_06 <- "Output_supercomputer/02_06_ELA"
dir_03_10 <- "Output_supercomputer/03_10_graphics_states_flow_flow_Spl_250508"

ra_l <- readRDS(sprintf("%s/Df_RA_Family.rds",dir_02_01))

bcolor <- readRDS(sprintf("%s/states_colvector.rds",dir_03_10))
##################
fb <- "Fungi"

ocmat <- readRDS(sprintf("%s/%s/input_ocmat_%s.rds",dir_02_06,fb,fb))
ramat <- ra_l[[fb]]

info <- readRDS(sprintf("%s/comp_sample_info_plant2.rds",dir_02_01))
sp_info <- readRDS(sprintf("%s/ELA_input_plant.rds",dir_02_01))
sa <- readRDS(sprintf("%s/%s/runSA_ST_full_%s_Family.rds",dir_02_06,fb,fb))
pltab <- unique(sp_info[rownames(ocmat),-c(1,ncol(sp_info))])

sslist <- readRDS(sprintf("%s/%s/graph_param_SStab_full_%s_Family.rds",dir_02_06,fb,fb))
#each plant

samp_basin <- NULL
for(pl in 1:nrow(pltab)){#pl <- 4
  
  if(sum(pltab[pl,]==1)){
    plnam <- colnames(pltab)[which(pltab[pl,]==1)]
  }else{
    plnam <- colnames(sp_info)[ncol(sp_info)]
  }
  
  pocm <- ocmat[which(rownames(ocmat) %in% info[info$plant==plnam,"ID"]),]
  
  sa2p <- sa2params(sa,as.numeric(pltab[pl,]))
  
  hgestp <- sa2p[[4]]
  jestp <- sa2p[[2]]

  samp_ss <- foreach(i=1:nrow(pocm),
                     .packages = c("rELA","foreach"),.combine = "c")%do%{
                       Bi(pocm[i,],
                          hgestp,jestp)[[1]]
                     }
  
  samp_basin <- rbind(samp_basin,data.frame(plant=plnam,
                                               sample=rownames(pocm),
                                               basin=sslist[samp_ss,"rename_SS"]))
  
}
saveRDS(samp_basin,file=sprintf("%s/samples_basin_info_%s_Family.rds",save.dir,fb))

ramat2 <- ramat[samp_basin$sample,]

pcoa_bray <- cmdscale(vegdist(ramat2,method="bray"),k=2,eig=TRUE)

bp_b<- cbind(samp_basin,PC1=pcoa_bray$points[samp_basin$sample,1],
      PC2=pcoa_bray$points[samp_basin$sample,2])

g <- ggplot(bp_b,aes(x=PC1,y=PC2,fill=basin,starshape=sprintf("*%s*",plant)))+
  geom_star(size=3,alpha=0.7)+
  theme_bw()+
  labs(  x=paste0("PCoA 1 (",round(pcoa_bray$eig[1]/sum(pcoa_bray$eig)*100,2),"%)" ),
         y=paste0("PCoA 2 (",round(pcoa_bray$eig[2]/sum(pcoa_bray$eig)*100,2),"%)"),
         fill="Basins in Fungal landscapes",starshape="Host background")+
  theme(plot.title = element_text(hjust = 0.5,size=16),
        aspect.ratio=1,
        axis.title = element_text(size=14),
        legend.text = element_markdown(size=12),
        axis.text = element_text(size=12),
        legend.title = element_text(size=14))+
  scale_fill_manual(values=bcolor)+
  scale_starshape_manual(values=c(4,13,15,30,23,5))+
  guides(fill=guide_legend(override.aes=list(size=5,alpha=1,starshape=15)))

       
g       


pcoa_jac <- cmdscale(vegdist(ocmat,method="jaccard"),k=2,eig=TRUE)
bp_j<- cbind(samp_basin,PC1=pcoa_jac$points[samp_basin$sample,1],
      PC2=pcoa_jac$points[samp_basin$sample,2])
g2 <- ggplot(bp_j,aes(x=PC1,y=PC2,fill=basin,starshape=sprintf("*%s*",plant)))+
  geom_star(size=3,alpha=0.7,show.legend = F)+
  theme_bw()+
  labs(  x=paste0("PCoA 1 (",round(pcoa_jac$eig[1]/sum(pcoa_jac$eig)*100,2),"%)" ),
         y=paste0("PCoA 2 (",round(pcoa_jac$eig[2]/sum(pcoa_jac$eig)*100,2),"%)"),
         fill="Basins in Fungal landscapes",starshape="Host background")+
  theme(plot.title = element_text(hjust = 0.5,size=16),
        aspect.ratio=1,
        axis.title = element_text(size=14),
        legend.text = element_markdown(size=12),
        axis.text = element_text(size=12),
        legend.title = element_text(size=14))+
  scale_fill_manual(values=bcolor)+
  scale_starshape_manual(values=c(4,13,15,30,23,5))+
  guides(fill=guide_legend(override.aes=list(size=5,alpha=1,starshape=15)))

g2

g2+g
ggsave(sprintf("%s/PCoA_%s_Family_samples_basin_host.pdf",save.dir,fb),
       g2+g,
       width=13,height=7)

##################
fb <- "Prokaryota"

ocmat <- readRDS(sprintf("%s/%s/input_ocmat_%s.rds",dir_02_06,fb,fb))
ramat <- ra_l[[fb]]

info <- readRDS(sprintf("%s/comp_sample_info_plant2.rds",dir_02_01))
sp_info <- readRDS(sprintf("%s/ELA_input_plant.rds",dir_02_01))
sa <- readRDS(sprintf("%s/%s/runSA_ST_full_%s_Family.rds",dir_02_06,fb,fb))
pltab <- unique(sp_info[rownames(ocmat),-c(1,ncol(sp_info))])

sslist <- readRDS(sprintf("%s/%s/graph_param_SStab_full_%s_Family.rds",dir_02_06,fb,fb))
#each plant

samp_basin <- NULL
for(pl in 1:nrow(pltab)){#pl <- 4
  
  if(sum(pltab[pl,]==1)){
    plnam <- colnames(pltab)[which(pltab[pl,]==1)]
  }else{
    plnam <- colnames(sp_info)[ncol(sp_info)]
  }
  
  pocm <- ocmat[which(rownames(ocmat) %in% info[info$plant==plnam,"ID"]),]
  
  sa2p <- sa2params(sa,as.numeric(pltab[pl,]))
  
  hgestp <- sa2p[[4]]
  jestp <- sa2p[[2]]
  
  samp_ss <- foreach(i=1:nrow(pocm),
                     .packages = c("rELA","foreach"),.combine = "c")%do%{
                       Bi(pocm[i,],
                          hgestp,jestp)[[1]]
                     }
  
  samp_basin <- rbind(samp_basin,data.frame(plant=plnam,
                                            sample=rownames(pocm),
                                            basin=ifelse(samp_ss %in% sslist[,1],
                                                         sslist[match(samp_ss,sslist[,1]),"rename_SS"],
                                                         "Minor basins")
                                            ))
  
}
saveRDS(samp_basin,file=sprintf("%s/samples_basin_info_%s_Family.rds",save.dir,fb))

maj_b <- setdiff(unique(samp_basin$basin),"Minor basins")

samp_basin$basin <- factor(samp_basin$basin,
                             levels=c(maj_b[order(as.numeric(gsub("P_B","",maj_b)))],
                                      "Minor basins"))
ramat2 <- ramat[samp_basin$sample,]

pcoa_bray <- cmdscale(vegdist(ramat2,method="bray"),k=2,eig=TRUE)

bp_b<- cbind(samp_basin,PC1=pcoa_bray$points[samp_basin$sample,1],
             PC2=pcoa_bray$points[samp_basin$sample,2])

g <- ggplot(bp_b,aes(x=PC1,y=PC2,fill=basin,starshape=sprintf("*%s*",plant)))+
  geom_star(size=3,alpha=0.7)+
  theme_bw()+
  labs(  x=paste0("PCoA 1 (",round(pcoa_bray$eig[1]/sum(pcoa_bray$eig)*100,2),"%)" ),
         y=paste0("PCoA 2 (",round(pcoa_bray$eig[2]/sum(pcoa_bray$eig)*100,2),"%)"),
         fill="Basins in prokaryotic landscapes",starshape="Host background")+
  theme(plot.title = element_text(hjust = 0.5,size=16),
        aspect.ratio=1,
        axis.title = element_text(size=14),
        legend.text = element_markdown(size=12),
        axis.text = element_text(size=12),
        legend.title = element_text(size=14))+
  scale_fill_manual(values=bcolor)+
  scale_starshape_manual(values=c(4,13,15,30,23,5))+
  guides(fill=guide_legend(override.aes=list(size=5,alpha=1,starshape=15)))


g       


pcoa_jac <- cmdscale(vegdist(ocmat,method="jaccard"),k=2,eig=TRUE)
bp_j<- cbind(samp_basin,PC1=pcoa_jac$points[samp_basin$sample,1],
             PC2=pcoa_jac$points[samp_basin$sample,2])
g2 <- ggplot(bp_j,aes(x=PC1,y=PC2,fill=basin,starshape=sprintf("*%s*",plant)))+
  geom_star(size=3,alpha=0.7,show.legend = F)+
  theme_bw()+
  labs(  x=paste0("PCoA 1 (",round(pcoa_jac$eig[1]/sum(pcoa_jac$eig)*100,2),"%)" ),
         y=paste0("PCoA 2 (",round(pcoa_jac$eig[2]/sum(pcoa_jac$eig)*100,2),"%)"),
         fill="Basins in prokaryotic landscapes",starshape="Host background")+
  theme(plot.title = element_text(hjust = 0.5,size=16),
        aspect.ratio=1,
        axis.title = element_text(size=14),
        legend.text = element_markdown(size=12),
        axis.text = element_text(size=12),
        legend.title = element_text(size=14))+
  scale_fill_manual(values=bcolor)+
  scale_starshape_manual(values=c(4,13,15,30,23,5))+
  guides(fill=guide_legend(override.aes=list(size=5,alpha=1,starshape=15)))

g2

g2+g
ggsave(sprintf("%s/PCoA_%s_Family_samples_basin_host.pdf",save.dir,fb),
       g2+g,
       width=13,height=7)
