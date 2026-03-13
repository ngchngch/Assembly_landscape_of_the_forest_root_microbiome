











##########################################################################
#install.packages("/Users/ngch/Desktop/Now_using/data/sugaFB_ELA/packages/rELAv0.51_fujita240702.tar.gz")
#setwd("/Users/ngch/Desktop/Now_using/data/sugaFB_ELA")
set.seed(1234)


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
library("ggforce")



#ELA_prep_dir='02_01_ELA_prep_abundance_threshold'
#n.core=8
#dir_02_03='02_03_summarize_occ_taxa_th'





#########################################################################
save.dir <- "03_02_05_randELA_withRA_4step_boot"
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
  #cluster = makeCluster(threads)
  #registerDoParallel(cluster)
  #on.exit(stopCluster(cluster))
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
      ssp_d <-foreach(k = 1:nrow(ssp),.packages = "vegan",.combine = "c")%do%{#k <- 1
        res <- ssp[k,3]*vegdist(rbind(sstr[ssp[k,1],-1],
                                      sstr[ssp[k,2],-1]),
                                method="jaccard")[1]/(SS.itr)
        if(is.na(res)){res <- 0}
        return(res)
      }
      d_land[i-1] <- sum(ssp_d)
    }
    
    # foreach(l=1:length(ssprop),.combine = "rbind")%do%{
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


Bin_2sd <- function(df){
  binmat <- NULL
  for(i in 1:ncol(df)){#i <- 1
    x <- df[,i]
    x2 <- log(x[x>0])
    lth <- exp(mean(x2)-2*sd(x2))
    binmat <- cbind(binmat,matrix(ifelse(x<lth,0,1),ncol=1))
  }
  return(binmat)
}


blockSample <- function(mat,el,name){
  block_rand <- c()
  names <- c()
  randnames <- c()
  for(i in 1:length(unique(el))){#i <- 1
    names <- c(names,name[which(el %in% unique(el)[i])])
    randnames <- c(randnames,sample(name[which(el %in% unique(el)[i])]))
  }
  if(is.vector(mat)){
    block_rand <- mat[names,]
    names(block_rand) <- randnames
    block_rand <- block_rand[names(mat),]
  }else{
    block_rand <- mat[names,]
    rownames(block_rand) <- randnames
    block_rand <- block_rand[rownames(mat),]
  }
  
  return(list(matrix=block_rand,rownames=randnames))
}

########################################################################
tx_f <- readRDS("Base_data/Fungi/taxa_list_mod.rds")
tx_b <- as.data.frame(readRDS("Base_data/OTU_basedata_set/taxonomy_list_dada2.rds"))

####################################################
taxa <- "Family"
#read data
tb_gns <- readRDS(sprintf("%s/NoCLR_ExpVar_clrRA_target_taxa.rds",ELA_prep_dir))
ocmat <- list(Fungi=readRDS(sprintf("%s/ELA_input_ocmat_Fungi.rds",dir_02_05)),
              Prokaryota=readRDS(sprintf("%s/ELA_input_ocmat_Prokaryote.rds",dir_02_05)))

info <- readRDS(sprintf("%s/comp_sample_info_plant2.rds",ELA_prep_dir))
sp_info <- readRDS(sprintf("%s/ELA_input_plant.rds",ELA_prep_dir))

###################################################
target_mat <- matrix(c("Cladophialophora","Populus","Fungi","Fungi",
                      "Oidiodendron","Acer","Fungi","Fungi",
                      "Candidatus Udaeobacter","Juglans","Bacteria","Bacteria",
                      "Meliniomyces","Pinus","Bacteria","Fungi"),ncol=4,byrow=TRUE)


inputdir <- "threadlipper_analysis_series/240430_caverage_rarefaction"

ls<- c(list.files(paste(inputdir,target_mat[1,4],sep="/"), pattern = sprintf("matrix_list_%s",target_mat[1,1]),
                recursive = TRUE,full.names = TRUE),
                list.files(paste(inputdir,target_mat[2,4],sep="/"), pattern = sprintf("matrix_list_%s",target_mat[2,1]),
                recursive = TRUE,full.names = TRUE),
                list.files(paste(inputdir,target_mat[3,4],sep="/"), pattern = sprintf("matrix_list_%s",target_mat[3,1]),
                recursive = TRUE,full.names = TRUE),
                list.files(paste(inputdir,target_mat[4,4],sep="/"), pattern = sprintf("matrix_list_%s",target_mat[4,1]),
                recursive = TRUE,full.names = TRUE))
names(ls) <- sapply(strsplit(ls,split="/"),function(x)x[4])
cat(ls)
#select ELA target OTU
##OTU in more.than 100 samples (no relative abundance threshold) & targeted in full ELA

#setdiff(tgOTU,txs)
#list <- ls[txs]
#ELA parameters
tb_g <- readRDS(sprintf("%s/ExpVar_clrRA_target_taxa.rds",ELA_prep_dir))
tgOTU <- setdiff(c(colnames(tb_g$Fungi),colnames(tb_g$Prokaryota)),"Unidentified")

list <- ls#[tgOTU]

###################################################
set.seed(seed)

plmat <- unique(sp_info[,-c(1,ncol(sp_info))])
pnam <-c()
      ##for function test only check Betula Pinus Acer
      for(pl in 1:nrow(unique(sp_info[,-1]))){#pl <-6
        cat("=")
        plmat0 <- plmat[pl,]
        
        if(sum(plmat0[1,]==1)==1){
          pl_nam <- colnames(plmat0)[which(plmat0[1,]==1)]
        }else{
          if(all(plmat0[1,]==0)){
            pl_nam <-colnames(sp_info)[ncol(sp_info)]
          }
        }
        
        pnam[pl] <- pl_nam
        
      }
rownames(plmat) <- pnam

for(boot_id in 1:100){
for(sp in 1:length(list)){
#ELA parameters
qth <- 10^-5 # do not change!!
SS.itr <- 20000
pl_nam <- target_mat[sp,2]
#i = 1 finished
#for(i in 1:length(list)){#i <-1
      s2 <- proc.time()[3]
      
      nam <- names(list)[sp]
      pt <- list[sp]
      
      if(str_detect(pattern = "Fungi",string = pt)){
        tag_fb <- "Fungi"
        abmat <- tb_g[[tag_fb]]
        ramat <- tb_gns[[tag_fb]]}
      if(str_detect(pattern = "Bacteria",string = pt)){
        tag_fb <- "Prokaryota"
        abmat <- tb_g[[tag_fb]]
        ramat <- tb_gns[[tag_fb]]
      }
     
      ## -- ELA
      fb <- target_mat[sp,3]
      if(fb=="Fungi"){
       #ocmatf0 <- f2[,which(colSums(f2)>30)]
      }else{
        #ocmatf0 <- b2[,which(colSums(b2)>30)]
        fb=="Prokaryota"
      }
      
      dir.create(sprintf("%s/%s",save.dir,fb))
      part <- sprintf("%s/%s/%s",save.dir,fb,nam)
      dir.create(part)
      
      cat(sprintf("%s->%s\n",nam,fb))
      
      input <- readRDS(sprintf("%s/%s/%s/inputdata_%s_wo%s_boot%s.rds",dir_03_01_04,fb,nam,fb,nam,boot_id))
      ocmatf <- input[[1]]
      enmatf0 <- input[[2]]
      which_pres0 <- input[[3]]
      ramat0 <- input[[4]]
      info <- input[[5]]
      base_set <- cbind(enmatf0,ramat0,which_pres0)
     rownames(base_set) <-  make.unique(rownames(base_set))
      cat(rownames(base_set)[1:10])
      
      rab <- blockSample(base_set,
                     info[,"plant"],
                     rownames(base_set))$matrix


      enmatf <- cbind(RA=rab[,1],enmatf0[,-1])
      #rownames(enmatf) <- rownames(base_set)
      ramat <- rab[,1+ncol(enmatf0)]#;names(ramat)<- rownames(base_set)
      which_pres <- rab[,1+ncol(enmatf0)+1:ncol(which_pres0)]
      #rownames(which_pres) <- rownames(base_set)
        sa <- runSA(ocmat=as.matrix(ocmatf),enmat = enmatf,
                    qth=qth, rep=16, threads=n.core)
        #saveRDS(sa, file=sprintf("%s/runSA_%s_wo%s.rds",part,fb,nam))
     
      #sa <- readRDS("/Users/ngch/Desktop/Now_using/data/sugaFB_ELA/analysis_series/03_03_same_samp_ELA_withRA_Genus_240703/Fungi/Acidibacter/runSA_Fungi_woAcidibacter.rds")
      #make input start sets
      hg <- sa2params(sa)[[4]]
      state <- foreach(i=1:SS.itr,.combine="rbind")%do%{
        st <- runif(length(hg), 0, 2) |> as.integer()
      }
      rownames(state) <- sprintf("Start_%05d",1:SS.itr)
      #saveRDS(state, file=sprintf("%s/start_%s.rds",part,fb))
      #######################
      
      cat("processing(")
      cat(nrow(unique(sp_info[,-1])))
      cat(") |")
      
        psamp <- make.unique(info[info$plant == pl_nam,"ID"])
        ra <- enmatf[which(rownames(enmatf) %in% psamp),"RA"]
        ra_noclr <- ramat[which(rownames(enmatf) %in% psamp)]
        wpres <- which_pres[rownames(enmatf)[which(rownames(enmatf) %in% psamp)],nam]
        cat(ra[1:10])

        
        if(sum(ra_noclr>0) > 4){
          ra_perc <- quantile(ra[wpres==1],c(0.25,0.5,0.75))
          ran_perc <- quantile(ra_noclr[wpres==1],c(0.25,0.5,0.75))
          
          sprop <- SSchange(state=state,
                                  sa=sa,
                                  steps=3,
                                  RA_label=c("perc25","median","perc75"),
                                  env_cat=plmat[pl_nam,],reporting = FALSE,
                                  start=mean(ra[wpres==0]),
                                  range=c(ra_perc[1],ra_perc[2],ra_perc[3]),
                                  eq_steps = FALSE,
                                  SS.itr=SS.itr,threads=n.core)
          
          if(!sprop$skip){
            sstabl <- sprop$SStable
            md_sprop <- cbind(plant=pl_nam,ab=c(ran_perc[1],
                                                               ran_perc[2],
                                                               ran_perc[3]),
                                             std_ab_clr=c(ra_perc[1],
                                                      ra_perc[2],
                                                      ra_perc[3]),
                                             std_ab_clr0=c(mean(ra[ra_noclr==0]),NA,NA),
                                             sd_ab_clr=c(sd(abmat[rownames(ocmatf),nam]),
                                                        NA,NA),
                                             mean_ab_clr=c(mean(abmat[rownames(ocmatf),nam]),
                                                           NA,NA),
                                             sprop[["result"]])
            
          }else{
            sstabl <-NA
            md_sprop <- data.frame(plant=pl_nam,ab=c(NA,NA,NA),ra=c("perc25","median","perc75"),
                                                  std_ab_clr=c(NA,NA,NA),
                                                  std_ab_clr0=NA,
                                                  sd_ab_clr=c(NA,NA,NA),
                                                  mean_ab_clr=c(NA,NA,NA),
                                                  d_land=c(NA,NA,NA),d_even=c(NA,NA,NA))
            
          }
          
        }else{
          sprop <- NA
          sstabl <-NA
          md_sprop <- data.frame(plant=pl_nam,ab=c(NA,NA,NA),ra=c("perc25","median","perc75"),
                                                std_ab_clr=c(NA,NA,NA),
                                                std_ab_clr0=NA,
                                                sd_ab_clr=c(NA,NA,NA),
                                                mean_ab_clr=c(NA,NA,NA),
                                                d_land=c(NA,NA,NA),d_even=c(NA,NA,NA))
        }
     mdp <- cbind(Taxa=nam,md_sprop)
      saveRDS(mdp,sprintf("%s/ELA_SSprob_diff_%s_%s_boot%s_rand%s.rds",part,fb,nam,boot_id,seed))
      #write.table(mdp,sprintf("%s/check_SSprob_diff_%s_%s.txt",part,fb,nam))
      
      cat("|\n")
    
    cat(sprintf("Elapsed time %.2f sec\n", proc.time()[3] - s2))
  #}
  
}


}
