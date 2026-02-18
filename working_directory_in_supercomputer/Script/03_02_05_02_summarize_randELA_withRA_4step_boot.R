##########################################################################
#install.packages("/Users/ngch/Desktop/Now_using/data/sugaFB_ELA/packages/rELAv0.51_fujita240702.tar.gz")
#setwd("/Users/ngch/Desktop/Now_using/data/sugaFB_ELA")
set.seed(1234)
#########################################################################
save.dir <- "03_02_05_02_summarize_randELA_withRA_4step_boot"
dir.create(save.dir)






########################################################################
#read original functions
source("packages/01_1_function.R")

target_mat <- matrix(c("Cladophialophora","Populus","Fungi","Fungi",
                       "Oidiodendron","Acer","Fungi","Fungi",
                       "Candidatus Udaeobacter","Juglans","Bacteria","Bacteria",
                       "Meliniomyces","Pinus","Bacteria","Fungi"),ncol=4,byrow=TRUE)

rea_all <- NULL

for(i in 1:nrow(target_mat)){
for(boot_id in 1:100){
  or_ls <- list.files(paste(dir_03_01_04,target_mat[i,3],target_mat[i,1],sep="/"),
                      pattern = sprintf("ELA_SSprob_diff_%s_%s_boot%s.rds",target_mat[i,3],target_mat[i,1],boot_id),
                      recursive = TRUE,full.names = TRUE)
  ls<- list.files(paste(dir_03_02_05,target_mat[i,3],target_mat[i,1],sep="/"),
                  pattern = sprintf("ELA_SSprob_diff_%s_%s_boot%s_rand",target_mat[i,3],target_mat[i,1],boot_id),
                  recursive = TRUE,full.names = TRUE)
  
  or_res <- readRDS(or_ls)
  rand_res <- do.call(rbind,lapply(ls,readRDS))
  
  or_res$z_dland <- NA
  or_res[or_res$ra=="perc25",]$z_dland <- (or_res[or_res$ra=="perc25",]$d_land-mean(rand_res[or_res$ra=="perc25",]$d_land))/sd(rand_res[or_res$ra=="perc25",]$d_land)
  or_res[or_res$ra=="median",]$z_dland <- (or_res[or_res$ra=="median",]$d_land-mean(rand_res[or_res$ra=="median",]$d_land))/sd(rand_res[or_res$ra=="median",]$d_land)
  or_res[or_res$ra=="perc75",]$z_dland <- (or_res[or_res$ra=="perc75",]$d_land-mean(rand_res[or_res$ra=="perc75",]$d_land))/sd(rand_res[or_res$ra=="perc75",]$d_land)
  
  or_res$z_deven <- NA
  or_res[or_res$ra=="perc25",]$z_deven <- (or_res[or_res$ra=="perc25",]$d_even-mean(rand_res[or_res$ra=="perc25",]$d_even))/sd(rand_res[or_res$ra=="perc25",]$d_even)
  or_res[or_res$ra=="median",]$z_deven <- (or_res[or_res$ra=="median",]$d_even-mean(rand_res[or_res$ra=="median",]$d_even))/sd(rand_res[or_res$ra=="median",]$d_even)
  or_res[or_res$ra=="perc75",]$z_deven <- (or_res[or_res$ra=="perc75",]$d_even-mean(rand_res[or_res$ra=="perc75",]$d_even))/sd(rand_res[or_res$ra=="perc75",]$d_even)
  
  or_res$p_dland <- NA
  or_res[or_res$ra=="perc25",]$p_dland <- mean(rand_res[or_res$ra=="perc25",]$d_land >= or_res[or_res$ra=="perc25",]$d_land)
  or_res[or_res$ra=="median",]$p_dland <- mean(rand_res[or_res$ra=="median",]$d_land >= or_res[or_res$ra=="median",]$d_land)
  or_res[or_res$ra=="perc75",]$p_dland <- mean(rand_res[or_res$ra=="perc75",]$d_land >= or_res[or_res$ra=="perc75",]$d_land)
  
  or_res$p_deven <- NA
  or_res[or_res$ra=="perc25",]$p_deven <- ifelse(or_res[or_res$ra=="perc25",]$d_even>0,
                                                 mean(rand_res[or_res$ra=="perc25",]$d_even >= or_res[or_res$ra=="perc25",]$d_even),
                                                 mean(rand_res[or_res$ra=="perc25",]$d_even <= or_res[or_res$ra=="perc25",]$d_even))
  or_res[or_res$ra=="median",]$p_deven <- ifelse(or_res[or_res$ra=="median",]$d_even>0,
                                                 mean(rand_res[or_res$ra=="median",]$d_even >= or_res[or_res$ra=="median",]$d_even),
                                                 mean(rand_res[or_res$ra=="median",]$d_even <= or_res[or_res$ra=="median",]$d_even))
  or_res[or_res$ra=="perc75",]$p_deven <- ifelse(or_res[or_res$ra=="perc75",]$d_even>0,
                                                 mean(rand_res[or_res$ra=="perc75",]$d_even >= or_res[or_res$ra=="perc75",]$d_even),
                                                 mean(rand_res[or_res$ra=="perc75",]$d_even <= or_res[or_res$ra=="perc75",]$d_even))
  
  rea_all <- rbind(rea_all,cbind(taxa=target_mat[i,1],boot_id=boot_id,or_res))
}  
}

saveRDS(rea_all,file=paste0(save.dir,"/","Zconv_ELA_withRA_4step_boot_summary.rds"))

