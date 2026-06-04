











##########################################################################
#install.packages("/Users/ngch/Desktop/Now_using/data/sugaFB_ELA/packages/rELAv0.51_fujita240702.tar.gz")
#setwd("/Users/ngch/Desktop/Now_using/data/sugaFB_ELA")



library(ggplot2)
library(ggstar)
library(parallel)
library(foreach)
library(vegan)
library("Rcpp")
library("RcppArmadillo")
library("doParallel")
library('tidyverse')
library('gtools')
library('igraph')
library('RColorBrewer')
library("stringdist")
#library('scatterpie')
library("rELA")
library("ggstar")



#ELA_prep_dir='02_01_ELA_prep_abundance_threshold'
#n.core=8
#dir_02_03='02_03_summarize_occ_taxa_th'





#########################################################################
save.dir <- "05_03_randELA_withRA_4step_Nsp80_nonoise_260527"
dir.create(save.dir)

########################################################################
#read original functions
source("packages/01_1_function.R")

#from fujita
Rcpp::sourceCpp("packages/ELA_functions_v060.cpp")

##############################################

SSchange <- function(state,sa,env_grad,env_cat=NA,steps=4,eq_steps=TRUE,
                     start=NA,
                     RA_label=NA,
                     range=NA,SS.itr=20000,threads=1,reporting=TRUE){
  #env_grad <- enmatf[,1];env_cat <- sp_info[1,-c(1,ncol(sp_info))]#Betula
  #range <- enmatf[,1]
  cluster = makeCluster(threads)
  registerDoParallel(cluster)
  on.exit(stopCluster(cluster))
  ##===============================================##
  ##GradELA
  s2 <- proc.time()[3]
  if(reporting){cat('Start SSestimation with gradient factor\n')}
  
  if(is.na(RA_label[1])){
    RA_label <- seq(1,step,1)
  }
  
  if(eq_steps){
    
    if(is.na(range[1])){mi <- min(env_grad)}else{mi <- range[1]}
    if(is.na(range[1])){ma <- max(env_grad)}else{ma <- range[2]}
    
    dm <- (ma - mi)/(steps - 1)
    de <- seq(mi, ma, dm)
    
  }else{
    if(length(range)<steps){steps <- length(range)}
    ur <- ifelse(range>start,range,NA)
    RA_label <- RA_label[!is.na(ur)]
    
    ur2 <- ur[!is.na(ur)]
    if(length(ur2)>0){
      de <- c(start,ur2[order(ur2)])  
      RA_label <- RA_label[order(ur2)]
    }else{
      de <- NA
    }
  }
  
  if(!is.na(de[1])){
    ssprop <- foreach(i = 1:length(de),.packages = c("rELA","tidyr","doParallel","vegan"),
                      .combine="c")%do%{
                        #i <- 2
                        if(!is.na(env_cat[1])){
                          ee <- c(de[i],as.numeric(env_cat))  
                        }else{
                          ee <- de[i]
                        }
                        
                        param <- sa2params(sa,ee)
                        
                        hge <- param[[4]]
                        je <- param[[2]]
                        
                        ss <- SSestimate_given(hge, je, state)
                        
                        mst <- ss[,-ncol(ss)]
                        rownames(mst) <- rownames(state)
                        
                        id <- apply(mst, 1, paste, collapse='')
                        ssid <- apply(unique(mst),1,bin2id)
                        
                        umst <- data.frame(env=ee[1],id2=id[names(ssid)],as.numeric(table(id)),unique(mst))
                        dimnames(umst) <- list(ssid,c("env","id2","freq",names(hge)))
                        
                        names(ssid) <- apply(unique(mst), 1, paste, collapse='')
                        prop <- as.numeric(table(id)[names(ssid)])/SS.itr
                        return(list(list(df=data.frame(env_grad=de[i],ssid=ssid,
                                                       h1=vegan::diversity(prop),
                                                       prop=prop),
                                         simulation=id,
                                         ss_structure=umst)))
                      }
    
    sstr1 <- foreach(l=1:length(ssprop),.combine = "rbind")%do%{
      return(ssprop[[l]]$ss_structure)
    }
    sstr2 <- sstr1
    
    sstr <- unique(sstr1[,-c(1,3)])
    rownames(sstr) <- sstr[,"id2"]
    sdf <- foreach(l=1:length(ssprop),.combine = "rbind")%do%{
      return(ssprop[[l]]$df)
    }
    ssim <- foreach(l=1:length(ssprop),.combine = "cbind")%do%{
      return(ssprop[[l]]$simulation)
    }
    
    minss <- unique(ssim[,1])
    d_land <- c()
    for(i in 2:ncol(ssim)){#i <- 2
      ssp <- NULL
      for(j in 1:length(minss)){
        ssp_tb <- table(ssim[which(ssim[,1] %in% minss[j]),i])
        ssp <-  rbind(ssp,data.frame(ss1=minss[j],ss2=names(ssp_tb),
                                     count=as.numeric(ssp_tb)))
      }
      ssp_d <-foreach(k = 1:nrow(ssp),.packages = "vegan",.combine = "c")%dopar%{#k <- 1
        res <- ssp[k,3]*vegdist(rbind(sstr[ssp[k,1],-1],
                                      sstr[ssp[k,2],-1]),
                                method="jaccard")[1]/(SS.itr)
        if(is.na(res)){res <- 0}
        return(res)
      }
      d_land[i-1] <- sum(ssp_d)
    }
    
    # foreach(l=1:length(ssprop),.combine = "rbind")%dopar%{
    #   return(ssprop[[l]]$df)
    # }
    
    res <- spread(sdf,key=ssid,value=prop)
    
    res[is.na(res)] <- 0
    
    #delta evenness
    d_even <- c(res$h1[2:length(res$h1)]-res$h1[1])
    
    if(reporting){cat(sprintf("Elapsed time %.2f sec\n", proc.time()[3] - s2))}
    
    return(list(skip=FALSE,
                SStable=cbind(rownames(sstr2),sstr2),
                df=res,
                result=data.frame(ra=RA_label,
                                  d_land=d_land,
                                  d_even=d_even)))
    
  }else{
    return(list(skip=TRUE,
                SStable=NA,
                df=NA,
                result=data.frame(ra=NA,
                                  d_land=NA,
                                  d_even=NA)))
  }
  
}


########################################################################
qt_seq <- c(0.5)

sc <- readRDS(sprintf("%s/top3_seed_across.rds",dir_05_01_00))

for(r in 1:15){
seed <- sc[r,2]
conn <- sc[r,1]

set.seed(cur_seed)

cpath0 <- list.files(sprintf("%s/seed%s",dir_05_01,seed),pattern = "community_matrix",full.names=TRUE)
cpath <- cpath0[grep(paste0("conn",conn),cpath0)]
cnam <- sapply(strsplit(cpath,"_"),function(x){gsub(".rds","",x[length(x)])})
dir.create(sprintf("%s/seed%s",save.dir,seed))


for(sp in 1:80){
for(conn in 1:length(cpath)){
  comm_list <- readRDS(cpath[conn])
  cnam1 <- cnam[conn]
  dir.create(sprintf("%s/seed%s/%s",save.dir,seed,cnam1))
  
  bin_mat <- comm_list$bin_mat;dimnames(bin_mat) <- list(paste0("S",1:nrow(bin_mat)),paste0("Sp",1:ncol(bin_mat)))
  ab_mat <- comm_list$ab;dimnames(ab_mat) <- list(paste0("S",1:nrow(ab_mat)),paste0("Sp",1:ncol(ab_mat)))
  #Rcpp::sourceCpp("for_Github/packages/ELA_functions_v060.cpp")
  #sp <- 1
  
  ###################################################
  #ELA parameters
  qth <- 10^-5 # do not change!!
  SS.itr <- 20000
  
  part <- sprintf("%s/seed%s/%s/%s",save.dir,seed,cnam1,paste0("ID",sp))
  dir.create(part)
  
  cat(sprintf("%s\n",sp))
  
  #list[ocmatf, abmatf, enmatf, samplelabelf, specieslabelf, factorlabelf] <-
  r_seq <- sample(1:nrow(ab_mat))
  enmat_ns <- ab_mat[r_seq,sp]
  enmatf <- matrix(scale(ab_mat[r_seq,sp]),ncol=1)
  rownames(enmatf) <- rownames(ab_mat)
  #saveRDS(enmatf,sprintf("%s/check_enmatf.rds",save.dir))
  colnames(enmatf) <- "Abundance"
  
  bin_sp <- bin_mat[r_seq,sp]
  bin0 <- as.matrix(bin_mat[,-sp])
  bin <- bin0[,colMeans(bin0)>0.1&colMeans(bin0)<0.9]
  
  if(sum(bin_sp==1)>0&sum(bin_sp==0)>0){
sa <- runSA(ocmat= bin,enmat = enmatf,
              qth=qth, rep=16, threads=n.core)
  saveRDS(sa, file=sprintf("%s/runSA_wo%s_%s.rds",part,sp,cnam1))
  
  #sa <- readRDS("/Users/ngch/Desktop/Now_using/data/sugaFB_ELA/analysis_series/03_03_same_samp_ELA_withRA_Genus_240703/Fungi/Acidibacter/runSA_Fungi_woAcidibacter.rds")
  #make input start sets
  hg <- sa2params(sa)[[4]]
  state <- foreach(i=1:SS.itr,.combine="rbind")%do%{
    st <- runif(length(hg), 0, 2) |> as.integer()
  }
  rownames(state) <- sprintf("Start_%05d",1:SS.itr)
  #saveRDS(state, file=sprintf("%s/start_%s.rds",part,cnam1))
  
  ############
  

  ra_perc <- quantile(enmatf[bin_sp==1,],qt_seq)
  ran_perc <- quantile(enmat_ns[bin_sp==1],qt_seq)
  
    sprop <- SSchange(state=state,
                      sa=sa,
                      steps=length(qt_seq),
                      RA_label=paste0("perc",qt_seq*100),
                      reporting = FALSE,
                      start=min(enmatf),
                      range=ra_perc,
                      eq_steps = FALSE,
                      SS.itr=SS.itr,threads=n.core)
    
    
    if(!sprop$skip){
      sstabl <- sprop$SStable
      md_sprop <- cbind(id=sp,ab=ran_perc,
                        scale_ab=ra_perc,
                        scale_ab0=c(min(enmatf),rep(NA,length(qt_seq)-1)),
                        sd_ab=c(sd(enmat_ns),
                                rep(NA,length(qt_seq)-1)),
                        mean_ab=c(mean(enmat_ns),
                                  rep(NA,length(qt_seq)-1)),
                        sprop[["result"]])
      
    }else{
      sstabl <- sprop$SStable
      md_sprop <- data.frame(id=sp,
                             ab=rep(NA,length(qt_seq)),
                             scale_ab=rep(NA,length(qt_seq)),
                             scale_ab0=rep(NA,length(qt_seq)),
                             sd_ab=rep(NA,length(qt_seq)),
                             mean_ab=rep(NA,length(qt_seq)),
                             d_land=rep(NA,length(qt_seq)),d_even=rep(NA,length(qt_seq)))
      
    }
  
  }else{
md_sprop <- data.frame(id=sp,
                             ab=rep(NA,length(qt_seq)),
                             scale_ab=rep(NA,length(qt_seq)),
                             scale_ab0=rep(NA,length(qt_seq)),
                             sd_ab=rep(NA,length(qt_seq)),
                             mean_ab=rep(NA,length(qt_seq)),
                             d_land=rep(NA,length(qt_seq)),d_even=rep(NA,length(qt_seq)))
      
  }
    
  
  
  mdp <- md_sprop
  saveRDS(mdp,sprintf("%s/ELA_SSprob_diff_ID%s_%s_rand%s.rds",part,sp,cnam1,cur_seed))
   #write.table(mdp,sprintf("%s/check_SSprob_diff_ID%s.txt",part,nam))
  
  cat("|\n")  
}}

}
