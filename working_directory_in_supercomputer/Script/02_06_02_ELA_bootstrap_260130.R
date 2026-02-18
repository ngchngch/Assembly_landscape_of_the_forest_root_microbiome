











##########################################################################
set.seed(1234)


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
save.dir <- "02_06_02_ELA_bootstrap_260130"
dir.create(save.dir)


########################################################################
#read original functions
source("packages/01_1_function.R")
Rcpp::sourceCpp("packages/ELA_functions_v060.cpp")

########################################################################
tx_f <- readRDS("Base_data/Fungi/taxa_list_mod.rds")
tx_b <- as.data.frame(readRDS("Base_data/OTU_basedata_set/taxonomy_list_dada2.rds"))

taxa_col <- readRDS("color/cols_Family_02_Barplot_240720.rds")

####################################################
taxa <- "Family"

#ELA_prep_dir <- "analysis_series_SpCom/02_ELA_241216/02_01_ELA_prep_abundance_threshold"
#dir_02_05 <- "analysis_series_SpCom/02_ELA_241216/02_05_summary_taxa_select"

#read data
rel_df <- readRDS(sprintf("%s/Df_RA_%s.rds",ELA_prep_dir,taxa))
ocmat <- list(Fungi=readRDS(sprintf("%s/ELA_input_ocmat_Fungi.rds",dir_02_05)),
              Prokaryota=readRDS(sprintf("%s/ELA_input_ocmat_Prokaryote.rds",dir_02_05)))

info <- readRDS(sprintf("%s/comp_sample_info_plant2.rds",ELA_prep_dir))
sp_info <- readRDS(sprintf("%s/ELA_input_plant.rds",ELA_prep_dir))


###################################################
#ELA parameters
qth <- 10^-5 # do not change!!
SS.itr <- 150000
fb <-"Fungi"
dir <- sprintf("%s/%s",save.dir,fb)
dir.create(dir)
#
set.seed(boot_id)

true_sa <- readRDS(sprintf("%s/%s/runSA_ST_full_%s_Family.rds",dir_02_06,fb,fb))
true_basin <- readRDS(sprintf("%s/%s/detected_SS_%s.rds",dir_02_06,fb,fb))

  boot_order <- sample(1:nrow(ocmat[[fb]]),replace = TRUE)
  ocmatf <- ocmat[[fb]][boot_order,]
  
  dim(ocmatf)
  
  
  sa <- runSA(ocmat=as.matrix(ocmatf),
              enmat = sp_info[boot_order,-c(1,ncol(sp_info))],
              qth=qth, rep=1280, threads=n.core)
  
  pltab <- unique(sp_info[rownames(ocmatf),-c(1,ncol(sp_info))])
  
  
  plnam <- c()
  ene_df <- NULL
  bcomp_freq <- NULL

   for(pl in 1:nrow(pltab)){#pl <- 4
  
    if(sum(pltab[pl,]==1)){
    plnam[pl] <- colnames(pltab)[which(pltab[pl,]==1)]
  }else{
    plnam[pl] <- colnames(sp_info)[ncol(sp_info)]
  }
    
      
    ela_true <- readRDS(sprintf("%s/%s/ELA_full_%s_%s.rds",dir_02_06,fb,fb,plnam[pl]))
   full_ss <-  ela_true[[1]]
   fulltp <- unlist(ela_true[[3]]);fulltp <- fulltp[fulltp!="Inf"]
   
  pocm <- ocmatf[which(rownames(ocmatf) %in% info[info$plant==plnam[pl],"ID"]),]
  
  sa2p <- sa2params(sa,as.numeric(pltab[pl,]))
  
  hgestp <- sa2p[[4]]
  jestp <- sa2p[[2]]
  
  ss_energy <- c()
  for(ssid in 1:length(full_ss)){
    ss_energy[ssid] <- cEnergy(id2bin(full_ss[ssid],ncol(pocm)),hgestp,jestp)
  }
  
  ene_df <- rbind(ene_df,data.frame(plant=plnam[pl],state="Bottoms",
                                    state_id=full_ss,energy=ss_energy))
  
  if(length(fulltp)>0){
    tp_energy <- c()
    for(tpid in 1:length(fulltp)){
      tp_energy[tpid] <- cEnergy(id2bin(fulltp[tpid],ncol(pocm)),
                                 hgestp,jestp)
    }
    ene_df <- rbind(ene_df,data.frame(plant=plnam[pl],state="Boundaries",
                                      state_id=fulltp,energy=tp_energy))
    
  }

###################
sa2p_t <- sa2params(true_sa,as.numeric(pltab[pl,]))
  
  hgestp_t <- sa2p_t[[4]]
  jestp_t <- sa2p_t[[2]]

      state <- foreach(i=1:SS.itr,.combine="rbind")%do%{
        st <- runif(length(hgestp), 0, 2) |> as.integer()
      }
      rownames(state) <- sprintf("Start_%05d",1:SS.itr)

    tss <- SSestimate_given(hgestp_t, jestp_t, state)
    bss <- SSestimate_given(hgestp, jestp, state)

    mst <- tss[,-ncol(tss)]
    rownames(mst) <- rownames(state)
                     
    id <- apply(mst, 1, paste, collapse='')
    uid <- unique(id)

for(i_id in 1:length(uid)){
bcomp_freq <- rbind(bcomp_freq,c(plant=plnam[pl],ssid=bin2id(unique(mst)[i_id,]),comp_seq=uid[i_id],
                                 colMeans(bss[id==uid[i_id],-ncol(bss)])))
}
    

}

  saveRDS(ene_df,sprintf("%s/energies_%s_boot%s.rds",dir,fb,boot_id))
 saveRDS(bcomp_freq,sprintf("%s/basincomposition_freq_%s_boot%s.rds",dir,fb,boot_id))

#######################################################################
fb <-"Prokaryota"
  dir <- sprintf("%s/%s",save.dir,fb)
  dir.create(dir)
  #
  set.seed(boot_id)
  
  true_sa <- readRDS(sprintf("%s/%s/runSA_ST_full_%s_Family.rds",dir_02_06,fb,fb))
true_basin <- readRDS(sprintf("%s/%s/detected_SS_%s.rds",dir_02_06,fb,fb))

  boot_order <- sample(1:nrow(ocmat[[fb]]),replace = TRUE)
  ocmatf <- ocmat[[fb]][boot_order,]
  
  dim(ocmatf)
  
  
  sa <- runSA(ocmat=as.matrix(ocmatf),
              enmat = sp_info[boot_order,-c(1,ncol(sp_info))],
              qth=qth, rep=1280, threads=n.core)
  
  pltab <- unique(sp_info[rownames(ocmatf),-c(1,ncol(sp_info))])
  
  
  plnam <- c()
  ene_df <- NULL
  bcomp_freq <- NULL

   for(pl in 1:nrow(pltab)){#pl <- 4
  
    if(sum(pltab[pl,]==1)){
    plnam[pl] <- colnames(pltab)[which(pltab[pl,]==1)]
  }else{
    plnam[pl] <- colnames(sp_info)[ncol(sp_info)]
  }
    
      
    ela_true <- readRDS(sprintf("%s/%s/ELA_full_%s_%s.rds",dir_02_06,fb,fb,plnam[pl]))
   full_ss <-  ela_true[[1]]
   fulltp <- unlist(ela_true[[3]]);fulltp <- fulltp[fulltp!="Inf"]
   
  pocm <- ocmatf[which(rownames(ocmatf) %in% info[info$plant==plnam[pl],"ID"]),]
  
  sa2p <- sa2params(sa,as.numeric(pltab[pl,]))
  
  hgestp <- sa2p[[4]]
  jestp <- sa2p[[2]]
  
  ss_energy <- c()
  for(ssid in 1:length(full_ss)){
    ss_energy[ssid] <- cEnergy(id2bin(full_ss[ssid],ncol(pocm)),hgestp,jestp)
  }
  
  ene_df <- rbind(ene_df,data.frame(plant=plnam[pl],state="Bottoms",
                                    state_id=full_ss,energy=ss_energy))
  
  if(length(fulltp)>0){
    tp_energy <- c()
    for(tpid in 1:length(fulltp)){
      tp_energy[tpid] <- cEnergy(id2bin(fulltp[tpid],ncol(pocm)),
                                 hgestp,jestp)
    }
    ene_df <- rbind(ene_df,data.frame(plant=plnam[pl],state="Boundaries",
                                      state_id=fulltp,energy=tp_energy))
    
  }

###################
sa2p_t <- sa2params(true_sa,as.numeric(pltab[pl,]))
  
  hgestp_t <- sa2p_t[[4]]
  jestp_t <- sa2p_t[[2]]

      state <- foreach(i=1:SS.itr,.combine="rbind")%do%{
        st <- runif(length(hgestp), 0, 2) |> as.integer()
      }
      rownames(state) <- sprintf("Start_%05d",1:SS.itr)

    tss <- SSestimate_given(hgestp_t, jestp_t, state)
    bss <- SSestimate_given(hgestp, jestp, state)

    mst <- tss[,-ncol(tss)]
    rownames(mst) <- rownames(state)
                     
    id <- apply(mst, 1, paste, collapse='')
    uid <- unique(id)

for(i_id in 1:length(uid)){
bcomp_freq <- rbind(bcomp_freq,c(plant=plnam[pl],ssid=bin2id(unique(mst)[i_id,]),comp_seq=uid[i_id],
                                 colMeans(bss[id==uid[i_id],-ncol(bss)])))
}
    

}

  saveRDS(ene_df,sprintf("%s/energies_%s_boot%s.rds",dir,fb,boot_id))
 saveRDS(bcomp_freq,sprintf("%s/basincomposition_freq_%s_boot%s.rds",dir,fb,boot_id))
