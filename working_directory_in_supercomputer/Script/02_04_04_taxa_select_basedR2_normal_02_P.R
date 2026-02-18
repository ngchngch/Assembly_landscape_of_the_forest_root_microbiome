set.seed(1234)
library(rELA)
library(vegan)

library(foreach)
library(doParallel)
library(tidyr)
library(compositions)

source("packages/01_1_function.R")
#########################################################################
save.dir <- "02_04_04_taxa_select_basedR2_normal_02"
dir.create(save.dir)
#####################
#parameters
taxa <- "Family"





##select target taxa
tb_m <- readRDS(sprintf("%s/Df_RA_Family.rds",ELA_prep_dir))

info <- readRDS(sprintf("%s/comp_sample_info_plant2.rds",ELA_prep_dir))
sp_info <- readRDS(sprintf("%s/ELA_input_plant.rds",ELA_prep_dir))

#######
fb <- "Prokaryote"

res_f <- NULL
for(rel_th in seq(0,0.1,by=0.01)){
ocmat <- readRDS(sprintf("%s/ocmat_remove_M2SD_Prokaryota_Family_relth%s.rds",ELA_prep_dir,rel_th))

df2 <- readRDS(sprintf("%s/Permanova_sig_taxa_R2_relth%s.rds",dir_02_03_04,rel_th))
#######
relf2 <- tb_m$Prokaryota/rowSums(tb_m$Prokaryota)

df1 <- df2[df2$fb==fb,]
fdf <- df1[order(df1[,"R2"],decreasing=TRUE),]

if(nrow(fdf) >= nSp){
mn <- c()
for(pl in 1:length(unique(info$plant))){
  rel_df <- relf2[rownames(relf2) %in% rownames(info)[info$plant == unique(info$plant)[pl]],
                  ]
  ocmat2 <- ocmat[rownames(rel_df),fdf[1:nSp,"taxa"]]
  
  bc_dist <- as.matrix(vegdist(rel_df,method = "bray"))
  bc_dist[is.na(bc_dist)] <- 0
  diag(bc_dist) <- NA
  bc_dist[upper.tri(bc_dist)] <- NA
  jac_dist <- as.matrix(vegdist(ocmat2,method = "jaccard"))
  jac_dist[is.na(jac_dist)] <- 0
  diag(jac_dist) <- NA
  jac_dist[upper.tri(jac_dist)] <- NA
  
  col <- mantel(as.dist(bc_dist),as.dist(jac_dist),method="kendall",permutations=FALSE)
  mn[pl] <- mantel(as.dist(bc_dist),as.dist(jac_dist),method="kendall",permutations=FALSE)$statistic
  
}
  
res_f <- rbind(res_f,cbind(rel_th,nSp,mean(mn),matrix(mn,nrow=1)))

}
}
  
colnames(res_f) <- c("rel_th","nSp","mean_tau",unique(info$plant))

saveRDS(res_f,sprintf("%s/mantel_tau_%s_%s.rds",save.dir,fb,nSp))

