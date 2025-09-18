
##########################################################################
#function
set.seed(1234)

####################

#install.packages("mgcv")
library(ggplot2)
library(ggtext)
library(ggrepel)
library(parallel)
library('RColorBrewer')
#library('scatterpie')
library('patchwork')
library("rELA")
library('tidyverse')
library("doParallel")
library(foreach)
#library(renv)
library("Rcpp")
library("RcppArmadillo")
library(vegan)
library(ggalluvial)
library(dplyr)
library(ComplexHeatmap)

library(circlize)

#######

save.dir <- "Output/03_10_02_graphics_Fullstates_flow_Spl_250813"
dir.create(save.dir)



dir_03_10 <- "Output/03_10_graphics_states_flow_flow_Spl_250508"


df_f <- readRDS(sprintf("%s/gradland_pcoa_all_Fungi.rds",dir_03_10))
df_p <- readRDS(sprintf("%s/gradland_pcoa_all_Prokaryote.rds",dir_03_10))
colvec <- readRDS(sprintf("%s/states_colvector.rds",dir_03_10))

df <- unique(df_f[grep("F_B",df_f$state_id),c("PCo1","PCo2","ssid2")])

g_full <- ggplot(df,aes(x=PCo1,y=PCo2))+
  geom_point(aes(fill=ssid2),shape=22,size=2,show.legend = FALSE)+
  geom_text_repel(aes(label=ssid2),size=3)+
  labs(fill="Stable states",shape="",x="PCoA 1", y="PCoA 2")+
  theme_bw()+
  theme(legend.position = "bottom",
        aspect.ratio = 1,
        legend.direction = "vertical")+
  scale_fill_manual(values=colvec)

g_full

df2 <- unique(df_p[grep("P_B",df_p$state_id),c("PCo1","PCo2","ssid2")])

g_full2 <- ggplot(df2,aes(x=PCo1,y=PCo2))+
  geom_point(aes(fill=ssid2),shape=22,size=2,show.legend = FALSE)+
  geom_text_repel(aes(label=ssid2),size=3)+
  labs(fill="Stable states",shape="",x="PCoA 1", y="PCoA 2")+
  theme_bw()+
  theme(legend.position = "bottom",
        aspect.ratio = 1,
        legend.direction = "vertical")+
  scale_fill_manual(values=colvec)

g_full2

g_full+g_full2

ggsave(sprintf("%s/PCoA_fullSS.pdf",save.dir),
       h=4,w=8)
