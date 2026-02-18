set.seed(1234)
library(rELA)
library(vegan)

library(foreach)
library(doParallel)
library(tidyr)
library(compositions)

source("packages/01_1_function.R")
#########################################################################
save.dir <- "02_03_04_summarize_basedR2_normal_02"
dir.create(save.dir)
#####################
#parameters
taxa <- "Family"






##select target taxa
tb_m <- readRDS(sprintf("%s/Df_RA_Family.rds",ELA_prep_dir))


info <- readRDS(sprintf("%s/comp_sample_info_plant2.rds",ELA_prep_dir))
sp_info <- readRDS(sprintf("%s/ELA_input_plant.rds",ELA_prep_dir))

##R2 of each taxa

for(rel_th in seq(0,0.1,by=0.01)){
#ocmat <- readRDS(sprintf("%s/ocmat_remove_M2SD_Fungi_Family_relth%s.rds",ELA_prep_dir,rel_th))
lf <- list.files(sprintf("%s/rel_th%s",dir_02_02_04_F,rel_th),pattern = "R2_",
                 full.names = TRUE,recursive = TRUE)

lp <- list.files(sprintf("%s/rel_th%s",dir_02_02_04_P,rel_th),pattern = "R2_",
                 full.names = TRUE,recursive = TRUE)

df <- rbind(cbind(fb="Fungi",do.call(rbind,lapply(lf[grep("Fungi.rds",lf)],readRDS))),
            cbind(fb="Prokaryote",do.call(rbind,lapply(lp[grep("Prokaryota.rds",lp)],readRDS))))

df$p_BH <- p.adjust(as.numeric(df$raw_p),method = "BH")

df2 <- df#[df$p_BH<0.05,]

saveRDS(df2,sprintf("%s/Permanova_sig_taxa_R2_relth%s.rds",save.dir,rel_th))
write.csv(df2,sprintf("%s/Permanova_sig_taxa_R2_relth%s.csv",save.dir,rel_th))

}


