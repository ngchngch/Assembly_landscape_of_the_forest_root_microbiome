











##########################################################################
#function
set.seed(1234)


SSflow <-  function(state,sa,env_grad,env_cat=NA,steps=31,eq_steps=TRUE,
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
    RA_label <- seq(1,steps-1,1)
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
                                         energy=data.frame(env=de[i],
                                                           ssid=apply(unique(mst),1,bin2id),
                                                           energy=apply(unique(mst),1,cEnergy,hge,je)),
                                         simulation=id,
                                         ss_structure=umst)))
                      }
    
    ene <- foreach(l=1:length(ssprop),.combine = "rbind")%do%{
      return(ssprop[[l]]$energy)
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
    
  ssp <- NULL
  for(i in 1:(ncol(ssim)-1)){#i <- 2
    
    minss <- unique(ssim[,i])
    
    for(j in 1:length(minss)){
    ssp_tb <- table(ssim[which(ssim[,i] %in% minss[j]),i+1])
      ssp <-  rbind(ssp,data.frame(step=paste0(i,"->",i+1),ss1=minss[j],ss2=names(ssp_tb),
                                   count=as.numeric(ssp_tb)))
    }}
  
  
  
  return(list(sankey=ssim,alluvial=ssp,
              gradELA=ene,
              df=res,
              result=data.frame(ra=RA_label,
                                d_land=d_land,
                                d_even=d_even)))
  }else{
    return(list(sankey=NA,alluvial=NA,
                gradELA=NA,df=NA,
                result=NA))
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

####################
#setwd("/Users/ngch/Desktop/Now_using/data/sugaFB_ELA")
#install.packages("mgcv")
library(ggplot2)
library(ggtext)
library(ggnetwork)
library(parallel)
library('RColorBrewer')
#library('scatterpie')
library('patchwork')
library("rELA")
library('tidyverse')
library("doParallel")
library(foreach)
#library(renv)
library("Rcpp")
library("RcppArmadillo")
library(igraph)
  library(ggalluvial)
  library(dplyr)
  library(ggsankey)
  
  
#######





#ir_03_07 <- "03_07_graphics_Zconv_landchanges_biplot"



#dir_03_01 <- "03_01_ELA_withRA"









#########################################################################
save.dir <- "03_05_states_flow_diagram_250307"
dir.create(save.dir)

########################################################################
#read original functions
source("packages/01_1_function.R")

#from fujita
Rcpp::sourceCpp("packages/ELA_functions_v060.cpp")


####################################################
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
inputdir <- "From_localcomputer/01_caverage_rarefaction"

ls<- list.files(inputdir, pattern = "matrix_list",
                recursive = TRUE,full.names = TRUE)
names(ls) <- sapply(strsplit(ls,split="/"),function(x)x[4])

#select ELA target OTU
##OTU in more.than 100 samples (no relative abundance threshold) & targeted in full ELA

sa_list <- list.files(dir_03_01,pattern = "runSA_",full.names = TRUE,recursive = TRUE)
#setdiff(tgOTU,txs)
#list <- ls[txs]
#ELA parameters
tb_g <- readRDS(sprintf("%s/ExpVar_clrRA_target_taxa.rds",ELA_prep_dir))
tgOTU <- setdiff(c(colnames(tb_g$Fungi),colnames(tb_g$Prokaryota)),"Unidentified")

list <- ls[tgOTU]
###################################################
#ELA parameters
qth <- 10^-5 # do not change!!
SS.itr <- 20000

#i = 1 finished
#for(i in 1:length(list)){#i <-1
s2 <- proc.time()[3]

nam <- names(list)[spi]
pt <- list[spi]
df <- readRDS(pt)

f <- Taxa.mat(df$Fungi,tx_f,taxa);f <- f[which(rownames(f) %in% rownames(ocmat$Fungi)),
                                         which(colnames(f) %in% colnames(ocmat$Fungi))]
f1 <- f/rowSums(f)
f2 <- Bin_2sd(f1);dimnames(f2) <- dimnames(f1)

b <- Taxa.mat(df$Bacteria,tx_b,taxa);b <- b[which(rownames(b) %in% rownames(ocmat$Prokaryota)),
                                            which(colnames(b) %in% colnames(ocmat$Prokaryota))]
b1 <- b/rowSums(b)

b2 <- Bin_2sd(b1);dimnames(b2) <- dimnames(b1)


if(str_detect(pattern = "Fungi",string = pt)){
  tag_fb <- "Fungi"
  abmat <- tb_g[[tag_fb]]
  ramat <- tb_gns[[tag_fb]]}
if(str_detect(pattern = "Bacteria",string = pt)){
  tag_fb <- "Prokaryota"
  abmat <- tb_g[[tag_fb]]
  ramat <- tb_gns[[tag_fb]]
}
which_pres <- Bin_2sd(ramat[,colSums(ramat>0)>30]);dimnames(which_pres) <- dimnames(ramat[,colSums(ramat>0)>30])

## -- ELA
fb <- "Fungi"
ocmatf <- f2

dir.create(sprintf("%s/%s",save.dir,fb))
part <- sprintf("%s/%s/%s",save.dir,fb,nam)

dir.create(part)
dir.create(sprintf("%s/Sankey",part))
dir.create(sprintf("%s/Flow",part))


cat(sprintf("%s->%s\n",nam,fb))

#list[ocmatf, abmatf, enmatf, samplelabelf, specieslabelf, factorlabelf] <-
enmatf <- cbind(RA=scale(abmat[rownames(ocmatf),nam]),
                sp_info[rownames(ocmatf),-c(1,ncol(sp_info))])


# sa <- runSA(ocmat=as.matrix(ocmatf),enmat = enmatf,
#             qth=qth, rep=256, threads=n.core)
# saveRDS(sa, file=sprintf("%s/runSA_%s_wo%s.rds",part,fb,nam))
# 
sa <- readRDS(sa_list[intersect(grep(nam,sa_list,fixed = TRUE),grep(fb,sa_list))])
#make input start sets
hg <- sa2params(sa)[[4]]
state <- foreach(i=1:SS.itr,.combine="rbind")%do%{
  st <- runif(length(hg), 0, 2) |> as.integer()
}
rownames(state) <- sprintf("Start_%05d",1:SS.itr)
saveRDS(state, file=sprintf("%s/start_%s.rds",part,fb))
#######################

cat("processing(")
cat(nrow(unique(sp_info[,-1])))
cat(") |")

pnam <-c()
mxmin_ab <- list(NULL)
sflow <- list(NULL)
ssims <- list(NULL)
flow_data <- list(NULL)
gradland <- NULL
gela <- NULL

##for function test only check Betula Pinus Acer
for(pl in 1:nrow(unique(sp_info[,-1]))){#pl <-3
  cat("=")
  plmat <- unique(sp_info[,-c(1,ncol(sp_info))])[pl,]
  
  if(sum(plmat[1,]==1)==1){
    pl_nam <- colnames(plmat)[which(plmat[1,]==1)]
  }else{
    if(all(plmat[1,]==0)){
      pl_nam <-colnames(sp_info)[ncol(sp_info)]
    }
  }
  
  pnam[pl] <- pl_nam
  
  psamp <- info[info$plant == pl_nam,"ID"]
  ra <- enmatf[which(rownames(enmatf) %in% psamp),"RA"]
  ra_noclr <- ramat[rownames(enmatf)[which(rownames(enmatf) %in% psamp)],
                    nam]
  
  wpres <- which_pres[rownames(enmatf)[which(rownames(enmatf) %in% psamp)],
                      nam]
  
  
  if(sum(ra_noclr>0) > 4){
    ra_perc <- quantile(ra[wpres==1],c(0.5))
    
    sflow[[pl]] <- SSflow(state=state,
                          sa=sa,
                          steps=32,
                          RA_label=NA,
                          env_cat=plmat,reporting = FALSE,
                          range=c(mean(ra[wpres==0]),max(ra)),
                          eq_steps = TRUE,
                          SS.itr=SS.itr,threads=n.core)
    
    mxmin_ab[[pl]] <- (c(mean(ra[wpres==0]),max(ra))*sd(abmat[rownames(ocmatf),nam]))+mean(abmat[rownames(ocmatf),nam])
    
    
    dat <- sflow[[pl]]$alluvial
    
    clrab <- (sflow[[pl]]$df[,1]*sd(abmat[rownames(ocmatf),nam]))+mean(abmat[rownames(ocmatf),nam])
    
    gradland <- rbind(gradland,
                      cbind(plant=pl_nam,ab=clrab,
                            rbind(c(NA,NA,NA),
                                  sflow[[pl]][["result"]])))
    
    sp <- sflow[[pl]]$gradELA
    sp$env <- (sp$env*sd(abmat[rownames(ocmatf),nam]))+mean(abmat[rownames(ocmatf),nam])
    gela <- rbind(gela,
                  cbind(plant=pl_nam,sp))
    
    
    step <- 32
    # re-format into long data format
    dat2 = NULL
    for (s in 1:(step-1)) {
      tem <- dat %>%
        filter(step == paste0(s,"->", s+1)) %>%
        {data.frame(
          alluvium = paste(s, c(seq(nrow(.)), 
                                seq(nrow(.)))), # same alluvium for each pair
          time = c(rep(s+ifelse(s>1,0.0001,0), nrow(.)), 
                   rep(s+1, nrow(.))), # +0.0001 to avoid overlapping with previous section
          count = c(.$count, .$count),
          group = c(.$ss1, .$ss2)
        )
        }
      dat2 <- rbind(dat2, tem)
    }
    
    ga <- dat2 %>%
      filter(count>0) %>%
      ggplot(aes(x = time, y = count, stratum = group, alluvium = alluvium)) +
      geom_flow(aes(fill = group)) +
      geom_stratum(data = dat2[dat2$time%%1==0,], aes(fill = group), alpha = 0.3) +
      #geom_text(data = dat2[dat2$time%%1==0,], stat = "stratum", aes(label = after_stat(stratum))) +
      theme_bw()+
      theme(legend.position = "none")
    
    ggsave(sprintf("%s/Flow/Flow_diagram1_%s_%s_%s.png",part,fb,nam,pl_nam),ga)
    
    flow_data[[pl]] <- dat2
    ############################
    
    ssim <- sflow[[pl]]$sankey
    long_ssim <- NULL
    for(i in 1:(ncol(ssim)-1)){
      long_ssim <- rbind(long_ssim,data.frame(now_x=i,next_x=i+1,
                                              now_node=ssim[,i],
                                              next_node=ssim[,i+1]))
    }
    
    gs <- ggplot(as.data.frame(long_ssim), aes(x = now_x, 
                                               next_x = next_x, 
                                               node = now_node, 
                                               next_node = next_node,
                                               fill = factor(now_node))) +
      geom_sankey(flow.alpha = 0.75, node.color = 0.1) +
      scale_fill_viridis_d( alpha = 0.95) +
      theme_sankey(base_size = 16)+
      theme(legend.position = "none")
    
    ggsave(sprintf("%s/Sankey/Sankey_diagram_%s_%s_%s.png",part,fb,nam,pl_nam),gs)
    ssims[[pl]] <- long_ssim
  }else{
    sflow[[pl]] <- NA
    ssims[[pl]] <- NA
    mxmin_ab[[pl]] <-NA
    flow_data[[pl]] <- NA
  }
}

names(mxmin_ab) <- pnam
names(sflow) <- pnam
names(ssims) <- pnam
names(flow_data) <- pnam

saveRDS(mxmin_ab,sprintf("%s/MaxMin_ab_%s_%s.rds",part,fb,nam))
saveRDS(sflow,sprintf("%s/SSFlow_%s_%s.rds",part,fb,nam))
saveRDS(ssims,sprintf("%s/Sankey_data_%s_%s.rds",part,fb,nam))
saveRDS(flow_data,sprintf("%s/Flow_data_%s_%s.rds",part,fb,nam))

saveRDS(gela,sprintf("%s/GradELA_data_%s_%s.rds",part,fb,nam))

saveRDS(gradland,sprintf("%s/landscape_change_data_%s_%s.rds",part,fb,nam))


g <- ggplot(gela,aes(x=env,y=energy))+
  geom_line(aes(color=ssid),size=1)+
  geom_point(aes(color=ssid),size=2,shape=1)+
  theme_bw()+
  theme(aspect.ratio = 1)+
  facet_wrap(~plant,scales="free",ncol=2)

ggsave(sprintf("%s/gradELA_%s_%s.png",part,fb,nam),plot=g,w=7,h=9,dpi=300)

g1 <- ggplot(gradland,aes(x=ab,y=d_land))+
  geom_line(aes(color=plant),size=1)+
  geom_point(aes(color=plant),size=2,shape=1)+
  theme_bw()+
  theme(aspect.ratio = 1)
g2 <- ggplot(gradland,aes(x=ab,y=d_even))+
  geom_line(aes(color=plant),size=1)+
  geom_point(aes(color=plant),size=2,shape=1)+
  theme_bw()+
  theme(aspect.ratio = 1)

ggsave(plot=g1/g2 + plot_layout(guides = "collect") ,
       sprintf("%s/landscape_change_%s_%s.png",part,fb,nam),w=6,h=12,dpi=300)


cat("|\n")

## -- ELA
fb <- "Prokaryota"
ocmatf <- b2

dir.create(sprintf("%s/%s",save.dir,fb))
part <- sprintf("%s/%s/%s",save.dir,fb,nam)
dir.create(part)
dir.create(sprintf("%s/Sankey",part))
dir.create(sprintf("%s/Flow",part))

cat(sprintf("%s->%s\n",nam,fb))

#list[ocmatf, abmatf, enmatf, samplelabelf, specieslabelf, factorlabelf] <-
enmatf <- cbind(RA=scale(abmat[rownames(ocmatf),nam]),
                sp_info[rownames(ocmatf),-c(1,ncol(sp_info))])


# sa <- runSA(ocmat=as.matrix(ocmatf),enmat = enmatf,
#             qth=qth, rep=256, threads=n.core)
# saveRDS(sa, file=sprintf("%s/runSA_%s_wo%s.rds",part,fb,nam))
# 
sa <- readRDS(sa_list[intersect(grep(nam,sa_list,fixed = TRUE),grep(fb,sa_list))])
#make input start sets
hg <- sa2params(sa)[[4]]
state <- foreach(i=1:SS.itr,.combine="rbind")%do%{
  st <- runif(length(hg), 0, 2) |> as.integer()
}
rownames(state) <- sprintf("Start_%05d",1:SS.itr)
saveRDS(state, file=sprintf("%s/start_%s.rds",part,fb))
#######################

cat("processing(")
cat(nrow(unique(sp_info[,-1])))
cat(") |")

pnam <-c()
mxmin_ab <- list(NULL)
sflow <- list(NULL)
ssims <- list(NULL)
flow_data <- list(NULL)
gradland <- NULL
gela <- NULL

##for function test only check Betula Pinus Acer
for(pl in 1:nrow(unique(sp_info[,-1]))){#pl <-3
  cat("=")
  plmat <- unique(sp_info[,-c(1,ncol(sp_info))])[pl,]
  
  if(sum(plmat[1,]==1)==1){
    pl_nam <- colnames(plmat)[which(plmat[1,]==1)]
  }else{
    if(all(plmat[1,]==0)){
      pl_nam <-colnames(sp_info)[ncol(sp_info)]
    }
  }
  
  pnam[pl] <- pl_nam
  
  psamp <- info[info$plant == pl_nam,"ID"]
  ra <- enmatf[which(rownames(enmatf) %in% psamp),"RA"]
  ra_noclr <- ramat[rownames(enmatf)[which(rownames(enmatf) %in% psamp)],
                    nam]
  wpres <- which_pres[rownames(enmatf)[which(rownames(enmatf) %in% psamp)],
                      nam]
  
  if(sum(ra_noclr>0) > 4){
    ra_perc <- quantile(ra[wpres==1],c(0.5))
    
    sflow[[pl]] <- SSflow(state=state,
                          sa=sa,
                          steps=32,
                          RA_label=NA,
                          env_cat=plmat,reporting = FALSE,
                          range=c(mean(ra[wpres==0]),max(ra)),
                          eq_steps = TRUE,
                          SS.itr=SS.itr,threads=n.core)
    
    mxmin_ab[[pl]] <- (c(mean(ra[wpres==0]),max(ra))*sd(abmat[rownames(ocmatf),nam]))+mean(abmat[rownames(ocmatf),nam])
    
    
    dat <- sflow[[pl]]$alluvial
    
    clrab <- (sflow[[pl]]$df[,1]*sd(abmat[rownames(ocmatf),nam]))+mean(abmat[rownames(ocmatf),nam])
    
    gradland <- rbind(gradland,
                      cbind(plant=pl_nam,ab=clrab,
                            rbind(c(NA,NA,NA),
                                  sflow[[pl]][["result"]])))
    
    sp <- sflow[[pl]]$gradELA
    sp$env <- (sp$env*sd(abmat[rownames(ocmatf),nam]))+mean(abmat[rownames(ocmatf),nam])
    gela <- rbind(gela,
                  cbind(plant=pl_nam,sp))
    
    
    step <- 32
    # re-format into long data format
    dat2 = NULL
    for (s in 1:(step-1)) {
      tem <- dat %>%
        filter(step == paste0(s,"->", s+1)) %>%
        {data.frame(
          alluvium = paste(s, c(seq(nrow(.)), 
                                seq(nrow(.)))), # same alluvium for each pair
          time = c(rep(s+ifelse(s>1,0.0001,0), nrow(.)), 
                   rep(s+1, nrow(.))), # +0.0001 to avoid overlapping with previous section
          count = c(.$count, .$count),
          group = c(.$ss1, .$ss2)
        )
        }
      dat2 <- rbind(dat2, tem)
    }
    
    ga <- dat2 %>%
      filter(count>0) %>%
      ggplot(aes(x = time, y = count, stratum = group, alluvium = alluvium)) +
      geom_flow(aes(fill = group)) +
      geom_stratum(data = dat2[dat2$time%%1==0,], aes(fill = group), alpha = 0.3) +
      #geom_text(data = dat2[dat2$time%%1==0,], stat = "stratum", aes(label = after_stat(stratum))) +
      theme_bw()+
      theme(legend.position = "none")
    
    ggsave(sprintf("%s/Flow/Flow_diagram1_%s_%s_%s.png",part,fb,nam,pl_nam),ga)
    
    flow_data[[pl]] <- dat2
    ############################
    
    ssim <- sflow[[pl]]$sankey
    long_ssim <- NULL
    for(i in 1:(ncol(ssim)-1)){
      long_ssim <- rbind(long_ssim,data.frame(now_x=i,next_x=i+1,
                                              now_node=ssim[,i],
                                              next_node=ssim[,i+1]))
    }
    
    gs <- ggplot(as.data.frame(long_ssim), aes(x = now_x, 
                                               next_x = next_x, 
                                               node = now_node, 
                                               next_node = next_node,
                                               fill = factor(now_node))) +
      geom_sankey(flow.alpha = 0.75, node.color = 0.1) +
      scale_fill_viridis_d( alpha = 0.95) +
      theme_sankey(base_size = 16)+
      theme(legend.position = "none")
    
    ggsave(sprintf("%s/Sankey/Sankey_diagram_%s_%s_%s.png",part,fb,nam,pl_nam),gs)
    ssims[[pl]] <- long_ssim
  }else{
    sflow[[pl]] <- NA
    ssims[[pl]] <- NA
    mxmin_ab[[pl]] <-NA
    flow_data[[pl]] <- NA
  }
}

names(mxmin_ab) <- pnam
names(sflow) <- pnam
names(ssims) <- pnam
names(flow_data) <- pnam

saveRDS(mxmin_ab,sprintf("%s/MaxMin_ab_%s_%s.rds",part,fb,nam))
saveRDS(sflow,sprintf("%s/SSFlow_%s_%s.rds",part,fb,nam))
saveRDS(ssims,sprintf("%s/Sankey_data_%s_%s.rds",part,fb,nam))
saveRDS(flow_data,sprintf("%s/Flow_data_%s_%s.rds",part,fb,nam))

saveRDS(gela,sprintf("%s/GradELA_data_%s_%s.rds",part,fb,nam))

saveRDS(gradland,sprintf("%s/landscape_change_data_%s_%s.rds",part,fb,nam))

g <- ggplot(gela,aes(x=env,y=energy))+
  geom_line(aes(color=ssid),size=1)+
  geom_point(aes(color=ssid),size=2,shape=1)+
  theme_bw()+
  theme(aspect.ratio = 1)+
  facet_wrap(~plant,scales="free",ncol=2)

ggsave(sprintf("%s/gradELA_%s_%s.png",part,fb,nam),plot=g,w=7,h=9,dpi=300)

g1 <- ggplot(gradland,aes(x=ab,y=d_land))+
  geom_line(aes(color=plant),size=1)+
  geom_point(aes(color=plant),size=2,shape=1)+
  theme_bw()+
  theme(aspect.ratio = 1)
g2 <- ggplot(gradland,aes(x=ab,y=d_even))+
  geom_line(aes(color=plant),size=1)+
  geom_point(aes(color=plant),size=2,shape=1)+
  theme_bw()+
  theme(aspect.ratio = 1)

ggsave(plot=g1/g2 + plot_layout(guides = "collect") ,
       sprintf("%s/landscape_change_%s_%s.png",part,fb,nam),w=6,h=12,dpi=300)


cat("|\n")
