set.seed(1234)
library(rELA)
library(vegan)

library(foreach)
library(doParallel)
library(tidyr)
library(compositions)
library(gtools)







source("packages/01_1_function.R")
#########################################################################
save.dir <- "02_02_correspond_Comm_taxaOcc_Fungi"
dir.create(save.dir)
#####################
#ELA_prep_dir <- "/Users/ngch/Desktop/Now_using/data/sugaFB_ELA/analysis_series_SpCom/02_ELA_v060_241108/02_01_ELA_prep_abundance_threshold"
##select target taxa
unident_th <- 0.5

tb_m <- readRDS(sprintf("%s/Df_RA_Family.rds",ELA_prep_dir))

ocmat <- readRDS(sprintf("%s/ocmat_remove_M2SD_Fungi_Family.rds",ELA_prep_dir))


info <- readRDS(sprintf("%s/comp_sample_info_plant2.rds",ELA_prep_dir))
sp_info <- readRDS(sprintf("%s/ELA_input_plant.rds",ELA_prep_dir))

#######
fb <- "Fungi"
relf2 <- tb_m$Fungi/rowSums(tb_m$Fungi)

relf3 <- relf2[which(relf2[,"Unidentified"] < unident_th & rownames(relf2) %in% rownames(ocmat)),]

ocmat2 <- ocmat[rownames(relf3),]

ocmat3 <- ocmat2[,which(colSums(ocmat2)>30 & colSums(ocmat2) < nrow(ocmat2)-30)]

CTRL.t <- how(within = Within(type = "free"), 
              plots = Plots(type = "none"),
              blocks = paste0(info[rownames(ocmat3),c("plant")],
                              info[rownames(ocmat3),c("site")]), 
              nperm = 10000,
              observed = TRUE)

  exp_df <- data.frame(info[rownames(ocmat3),c("plant","site")],
                  Occ=ocmat3[,tx])
  perm <- adonis2(data=exp_df,
          relf3~plant+site+Occ,
          by = "margin",
          permutations=CTRL.t,
          parallel = n.core)
  
  R2 <- data.frame(taxa=colnames(ocmat3)[tx],
                       R2=perm$R2[3],raw_p=perm$`Pr(>F)`[3])

saveRDS(R2,sprintf("%s/R2_%s_%s.rds",save.dir,colnames(ocmat3)[tx],fb))