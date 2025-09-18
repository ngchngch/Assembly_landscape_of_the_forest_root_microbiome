set.seed(1234)
library(rELA)
library(vegan)

library(foreach)
library(doParallel)
library(tidyr)
library(compositions)

source("packages/01_1_function.R")
#########################################################################
save.dir <- "02_03_summarize_basedR2"
dir.create(save.dir)
#####################
#parameters
taxa <- "Family"






##select target taxa
tb_m <- readRDS(sprintf("%s/Df_RA_Family.rds",ELA_prep_dir))

ocmat <- readRDS(sprintf("%s/ocmat_remove_M2SD_Fungi_Family.rds",ELA_prep_dir))
info <- readRDS(sprintf("%s/comp_sample_info_plant2.rds",ELA_prep_dir))
sp_info <- readRDS(sprintf("%s/ELA_input_plant.rds",ELA_prep_dir))

##R2 of each taxa

lf <- list.files(dir_02_02_F,pattern = "R2_",
                 full.names = TRUE,recursive = TRUE)

lp <- list.files(dir_02_02_P,pattern = "R2_",
                 full.names = TRUE,recursive = TRUE)

df <- rbind(cbind(fb="Fungi",do.call(rbind,lapply(lf,readRDS))),
            cbind(fb="Prokaryote",do.call(rbind,lapply(lp,readRDS))))

df$p_BH <- p.adjust(as.numeric(df$raw_p),method = "BH")

df2 <- df[df$p_BH<0.05,]

saveRDS(df2,sprintf("%s/Permanova_sig_taxa_R2.rds",save.dir))
write.csv(df2,sprintf("%s/Permanova_sig_taxa_R2.csv",save.dir))
