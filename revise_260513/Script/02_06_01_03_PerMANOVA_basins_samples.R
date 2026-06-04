











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
save.dir <- "Analysis/02_06_01_03_PerMANOVA_basins_samples"
dir.create(save.dir)

dir_02_01 <- "Basedata/02_01_ELA_prep_abundance_threshold" 
dir_02_06 <- "Basedata/02_06_ELA"
dir_02_06_01_02 <- "Basedata/02_06_01_02_PCoA_basins_samples"

ra_l <- readRDS(sprintf("%s/Df_RA_Family.rds",dir_02_01))
info <- readRDS(sprintf("%s/comp_sample_info_plant2.rds",dir_02_01))

nperm <- 10000
##################
fb <- "Fungi"

#####
basin_group <- cbind(basin=paste0("F_B",1:7),
                     group=c("EcM","EcM","EcM","AM","AM","End","AM"))
rownames(basin_group) <- basin_group[,"basin"]
#####
ocmat <- readRDS(sprintf("%s/%s/input_ocmat_%s.rds",dir_02_06,fb,fb))
ramat <- ra_l[[fb]]

samp_basin <- readRDS(sprintf("%s/samples_basin_info_%s_Family.rds",dir_02_06_01_02,fb))

sb <- cbind(samp_basin,
            site=info[samp_basin$sample,"site2"],
            group=basin_group[samp_basin$basin,"group"])
rownames(sb) <- sb$sample
ramat2 <- ramat[samp_basin$sample,]

mat <- ramat2

d_mat <- vegdist(mat,method="bray")
ra_each <- adonis2(d_mat~basin,data=sb[rownames(mat),],
        permutations=how(
          within = Within(type = "free"),
          plots = Plots(type = "none"),
          blocks = paste0(
            sb[rownames(mat), c("plant")],
            sb[rownames(mat), c("site2")]
          ),
          nperm = nperm, # In original analysis, we set nperm = 10000
          observed = TRUE
        ),
        parallel = 8)
saveRDS(ra_each,sprintf("%s/adonis_basin_each_%s.rds",save.dir,fb))

permdisp_basin <- betadisper(d_mat,sb[rownames(mat),"basin"])
saveRDS(permdisp_basin,sprintf("%s/permdisp_basin_each_%s.rds",save.dir,fb))

sink(sprintf("%s/PerMANOVA_basin_samples_%s.txt",save.dir,fb))
ra_each
anova(permdisp_basin)
sink()


####
mat <- ocmat
d_mat <- vegdist(mat,method="jaccard")

oc_each <- adonis2(d_mat~basin,data=sb[rownames(mat),],
        permutations=how(
          within = Within(type = "free"),
          plots = Plots(type = "none"),
          blocks = paste0(
            sb[rownames(mat), c("plant")],
            sb[rownames(mat), c("site2")]
          ),
          nperm = nperm, # In original analysis, we set nperm = 10000
          observed = TRUE
        ),
        parallel = 8)

saveRDS(oc_each,sprintf("%s/adonis_oc_each_%s.rds",save.dir,fb))

permdisp_oc <- betadisper(d_mat,sb[rownames(mat),"basin"])
saveRDS(permdisp_oc,sprintf("%s/permdisp_oc_each_%s.rds",save.dir,fb))

sink(sprintf("%s/PerMANOVA_basin_samples_ocmat_%s.txt",save.dir,fb))
oc_each
anova(permdisp_oc)
sink()

##################
fb <- "Prokaryota"

#####
#####
ocmat <- readRDS(sprintf("%s/%s/input_ocmat_%s.rds",dir_02_06,fb,fb))
ramat <- ra_l[[fb]]

samp_basin <- readRDS(sprintf("%s/samples_basin_info_%s_Family.rds",dir_02_06_01_02,fb))

sb <- cbind(samp_basin,
            site=info[samp_basin$sample,"site2"])
rownames(sb) <- sb$sample
ramat2 <- ramat[samp_basin$sample,]

mat <- ramat2

d_mat <- vegdist(mat,method="bray")
ra_each <- adonis2(d_mat~basin,data=sb[rownames(mat),],
                   permutations=how(
                     within = Within(type = "free"),
                     plots = Plots(type = "none"),
                     blocks = paste0(
                       sb[rownames(mat), c("plant")],
                       sb[rownames(mat), c("site2")]
                     ),
                     nperm = nperm, # In original analysis, we set nperm = 10000
                     observed = TRUE
                   ),
                   parallel = 8)
saveRDS(ra_each,sprintf("%s/adonis_basin_each_%s.rds",save.dir,fb))

permdisp_basin <- betadisper(d_mat,sb[rownames(mat),"basin"])
saveRDS(permdisp_basin,sprintf("%s/permdisp_basin_each_%s.rds",save.dir,fb))

sink(sprintf("%s/PerMANOVA_basin_samples_%s.txt",save.dir,fb))
ra_each
anova(permdisp_basin)
sink()


####
mat <- ocmat
d_mat <- vegdist(mat,method="jaccard")

oc_each <- adonis2(d_mat~basin,data=sb[rownames(mat),],
                   permutations=how(
                     within = Within(type = "free"),
                     plots = Plots(type = "none"),
                     blocks = paste0(
                       sb[rownames(mat), c("plant")],
                       sb[rownames(mat), c("site2")]
                     ),
                     nperm = nperm, # In original analysis, we set nperm = 10000
                     observed = TRUE
                   ),
                   parallel = 8)

saveRDS(oc_each,sprintf("%s/adonis_oc_each_%s.rds",save.dir,fb))

permdisp_oc <- betadisper(d_mat,sb[rownames(mat),"basin"])
saveRDS(permdisp_oc,sprintf("%s/permdisp_oc_each_%s.rds",save.dir,fb))

sink(sprintf("%s/PerMANOVA_basin_samples_ocmat_%s.txt",save.dir,fb))
oc_each
anova(permdisp_oc)
sink()





