











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
setwd("/Users/ngch/Desktop/Now_using/data/sugaFB_ELA/revise_260513")
save.dir <- "Analysis/02_06_01_04_summarize_PerMANOVA_basins_samples"
dir.create(save.dir)

dir_02_06_01_03 <- "Analysis/02_06_01_03_PerMANOVA_basins_samples" 

res <- list.files(dir_02_06_01_03,pattern="adonis_basin_each",full.names = TRUE) %>%
  lapply(function(x){x2 <- readRDS(x)
  data.frame(x2$Df[[1]],x2$R2[[1]],x2$`Pr(>F)`[[1]])
  }) %>%
  bind_rows()

res2 <- list.files(dir_02_06_01_03,pattern="adonis_oc_each",full.names = TRUE) %>%
  lapply(function(x){x2 <- readRDS(x)
  data.frame(x2$Df[[1]],x2$R2[[1]],x2$`Pr(>F)`[[1]])
  }) %>%
  bind_rows()


res3 <- list.files(dir_02_06_01_03,pattern="permdisp_basin_each",full.names = TRUE) %>%
  lapply(function(x){x2 <- anova(readRDS(x))
  data.frame(x2$Df[[1]],x2$`F value`[[1]],x2$`Pr(>F)`[[1]])
  }) %>%
  bind_rows()

res4 <- list.files(dir_02_06_01_03,pattern="permdisp_oc_each",full.names = TRUE) %>%
  lapply(function(x){x2 <- anova(readRDS(x))
  data.frame(x2$Df[[1]],x2$`F value`[[1]],x2$`Pr(>F)`[[1]])
  }) %>%
  bind_rows()

merge_res <- cbind(c("Fungi","Prokaryotes","Fungi","Prokaryotes"),
      rbind(cbind(res,res3[,2:3]),cbind(res2,res4[,2:3])))
colnames(merge_res) <- c("FB","Df","R2","p_PerMANOVA","F_value_PerMDISP","p_PerMDISP")
merge_res$p_BH_PerMANOVA <- p.adjust(merge_res[[4]],"BH")
merge_res$p_BH_PerMDISP <- p.adjust(merge_res[[6]],"BH")

saveRDS(merge_res,sprintf("%s/merge_res.rds",save.dir))
write.csv(merge_res[,c(1,2,3,5,7,8)],sprintf("%s/merge_res.csv",save.dir),row.names = FALSE)

