











##########################################################################
set.seed(1234)

library(ggplot2)
library(ggstar)
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
library("ggtext")


library("ggforce")







#########################################################################
save.dir <- "02_06_ELA"
dir.create(save.dir)


########################################################################
#read original functions
source("packages/01_1_function.R")


showDG_mod <- function(ela, oc, label="",SS_colmat,
                       na.color="black",minor.color="gray50",fontsize=5,
                       annot_adj=c(0.75, 2.00)){
  if(length(ela[[1]])>1){
    s <- ncol(oc)
    grobj <- GraphObj(ela)
    
    DG_mod(grobj, s,ss_list=ela[[1]], DG_sample=label,
           annot_adj=c(annot_adj[1], annot_adj[2]),SS_colmat=SS_colmat,
           na.color=na.color,minor.color=minor.color,fontsize=fontsize)
  }else{
    return(cat("only one stable state found\n"))
  }}

DG_mod <- function(
    grobj, s,ss_list, DG_sample, annot_adj=c(0.75, 2.00),
    SS_colmat,na.color="black",minor.color="gray90",fontsize=5){
  # grobj <- GraphObj(ela);s <- ncol(ocmat[[fb]]);SS_colmat=gpr;label=sprintf("%s %s",fb,pl[i]);na.color="gray90"
  range_ <- c(min(grobj$energy), max(grobj$energy))
  grobj_pre <- grobj
  grobj_pre$len_cee <- lapply(
    grobj_pre$cee, function(x){length(x)}) %>% as.integer
  #d_list <- unlist(map(1: s, function(x){paste0('d', x)}))
  grobj_ <- cbind(
    grobj_pre,
    do.call(
      cbind,
      grobj_pre$point %>% map(function(x){id2bin(x, s)})) %>% t #%>%
    #'colnames<-'(d_list)
  )
  jun_pre <- grobj_[grobj_$len_cee == 1,]
  jun <- jun_pre[order(jun_pre$nodes2xposi),]
  jen_pre <- grobj_[grobj_$len_cee > 1,]
  jen <- jen_pre[order(jen_pre$nodes2xposi),]
  
  # grobj to plot
  grobj_to_plot <- data.frame()
  if(nrow(grobj_) != 1){
    for(i in 1: (nrow(grobj_))){
      aa <- grobj_[i,]
      aa$point_str <- SS_colmat[aa$point,"rename_SS"]
      bb_pre <- grobj_[grobj_$cee %>% map(
        function(x){all(unlist(aa$cee) %in% unlist(x))}) %>% unlist(),]
      bb <- bb_pre[bb_pre$ccc_str != aa$ccc_str,][1,]
      bb$point_str <- rep("", length(bb$point))
      # bb$point_str <- paste0('C', bb$point)
      between_aa_bb <- aa
      between_aa_bb$energy <- bb$energy
      between_aa_bb$point_str <- ''
      line_break <- between_aa_bb
      line_break$nodes2xposi <- NA
      line_break$energy <- NA
      grobj_to_plot_ <- rbind(aa, between_aa_bb, bb, line_break)
      grobj_to_plot_$id <- i
      grobj_to_plot <- grobj_to_plot %>% rbind(grobj_to_plot_)
    }
  }
  
  grobj_to_plot$point_str[grepl("c(",grobj_to_plot$ccc_str,fixed = TRUE)] <- ""
  grobj_to_plot$point_str[is.na(grobj_to_plot$point_str)] <- "otherstate"
  
  pick <- which(grobj_to_plot$point_str != "" )
  # pick <- intersect(which(grobj_to_plot$point_str != "" ),
  #                    grep("c(",grobj_to_plot$ccc_str,invert=TRUE,fixed=TRUE))
  # 
  #col[is.na(col)] <- na.color
  
  minor <- intersect(which(grobj_to_plot$point_str == "otherstate" ),pick)
  grobj_to_plot$shape <-  ""
  grobj_to_plot$shape[pick] <- "SS"
  grobj_to_plot$shape[minor] <- "MinorSS"
  
  grobj_to_plot$size <-  ""
  grobj_to_plot$size[pick] <- "SS"
  grobj_to_plot$size[minor] <- "MinorSS"
  
  grobj_to_plot$size <- factor(grobj_to_plot$size, levels=c("SS", "MinorSS", ""))
  grobj_to_plot$shape <- factor(grobj_to_plot$shape, levels=c("SS", "MinorSS", ""))
  
  grobj_to_plot$point_str[minor] <- " "
  grobj_to_plot$point_str[grobj_to_plot$point_str=="otherstate"] <- ""
  grobj_to_plot$point_str2 <- factor(grobj_to_plot$point_str,
                                     levels=c(setdiff(unique(grobj_to_plot[,"point_str"]),
                                                      c(" ", "")),
                                              " ", ""))
  
  col <- c(SS_colmat[match(setdiff(unique(grobj_to_plot[,"point_str"]),
                                   c(" ", "")),SS_colmat$rename_SS),"color"],minor.color,na.color)
  # scatter with annot and line
  g <- ggplot(
    grobj_to_plot,
    aes(x=nodes2xposi, y=energy))
  g <- g + geom_path() +
    geom_point(aes(color=point_str2,shape=shape,size=size),stroke=2,show.legend = F) +
    xlab("") +
    ylab("Energy")+
    #labs(title=DG_sample)+
    geom_text(hjust=annot_adj[1], vjust=annot_adj[2],size=fontsize,
              aes(fontface=2, label=point_str)) +
    theme_bw() +
    theme(axis.title = element_text(size=18),
          axis.text.y = element_text(size=12),
          aspect.ratio = 1,
          panel.background = element_rect(fill = "white", color = "black", size = 2))+
    scale_color_manual(values=col[which(levels(grobj_to_plot$point_str2) %in% unique(grobj_to_plot$point_str2))])+
    scale_shape_manual(values=c(16,16, 18)[which(levels(grobj_to_plot$shape) %in% unique(grobj_to_plot$shape))]) +
    scale_size_manual(values=c(5,3, 1)[which(levels(grobj_to_plot$size) %in% unique(grobj_to_plot$size))]) +
    coord_cartesian(ylim=c(min(grobj_to_plot$energy,na.rm = TRUE)-0.05*(max(grobj_to_plot$energy,na.rm = TRUE)-min(grobj_to_plot$energy,na.rm = TRUE)),NA))
  
  # pie
  #g <- g + geom_scatterpie(
  #  aes(x=nodes2xposi, y=energy),
  #  data=grobj_, cols=d_list)
  
  # plot
  g <- g + ggtitle(DG_sample)
  plot(g)
  return(NULL)}

########################################################################
tx_f <- readRDS("Base_data/Fungi/taxa_list_mod.rds")
tx_b <- as.data.frame(readRDS("Base_data/OTU_basedata_set/taxonomy_list_dada2.rds"))

taxa_col <- readRDS("color/cols_Family_02_Barplot_240720.rds")

####################################################
taxa <- "Family"



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
maj_th <- 0.01


  ocmatf <- ocmat[[fb]]
  
  dim(ocmatf)
  saveRDS(ocmatf,sprintf("%s/input_ocmat_%s.rds",dir,fb))
  
  
  sa <- runSA(ocmat=as.matrix(ocmatf),
              enmat = sp_info[rownames(ocmatf),-c(1,ncol(sp_info))],
              qth=qth, rep=1280, threads=n.core)
  
  saveRDS(sa, file=sprintf("%s/runSA_ST_full_%s_%s.rds",dir,fb,taxa))
  
  pltab <- unique(sp_info[rownames(ocmatf),-c(1,ncol(sp_info))])
  
  enbr <- NULL
  SS_samp <- NULL
  detected_ss <- NULL
  detected_ss_all <- NULL
  ela <- list(NULL)
  
  for(pl in 1:nrow(pltab)){#pl <- 4
  
    if(sum(pltab[pl,]==1)){
    plnam <- colnames(pltab)[which(pltab[pl,]==1)]
  }else{
    plnam <- colnames(sp_info)[ncol(sp_info)]
  }
    
  pocm <- ocmatf[which(rownames(ocmatf) %in% info[info$plant==plnam,"ID"]),]
  
  sa2p <- sa2params(sa,as.numeric(pltab[pl,]))
  
  hgestp <- sa2p[[4]]
  jestp <- sa2p[[2]]
  
  png(sprintf("%s/h_histgram_%s_%s.png",dir,fb,plnam),width=10,height=10,units="in",res=300)
  hist(sa2p[[1]], breaks=100, col="grey", main="intercept", xlab="h", ylab="Frequency")
  dev.off()
  
  png(sprintf("%s/J_histgram_%s_%s.png",dir,fb,plnam),width=10,height=10,units="in",res=300)
  hist(jestp, breaks=100, col="grey", main="Associations", xlab="J", ylab="Frequency")
  dev.off()
  
  ela[[pl]] <- ELA(sa, env=as.numeric(pltab[pl,]),
             SS.itr=SS.itr, FindingTip.itr=10000, # <- the number of steps for finding stable states and tipping points (basically no need to change)
             threads=n.core, reporting=TRUE)
  
  saveRDS(ela[[pl]], file=sprintf("%s/ELA_full_%s_%s.rds",dir,fb,plnam))
  
  ss <- ela[[pl]][[1]]
  eb <- c()
  
  for(i in 1:length(ss)){
    e.stable.np <- ela[[pl]][[2]][i]
    stablestates.np <- ela[[pl]][[1]]
    tippingen.np <- ela[[pl]][[4]]
    e.tipping.np <- min(unlist(c(as.data.frame(tippingen.np)[stablestates.np==ss[i],],
                                 as.data.frame(tippingen.np)[,stablestates.np==ss[i]])))
    eb[i] <- e.tipping.np - e.stable.np
  };rm(i)
  
  
  if(length(eb)>1){
    
    pdf(sprintf("%s/Bdepth_histgram_%s_%s.pdf",dir,fb,plnam),width=100,height=100)
    hist(eb/max(eb), breaks=100, col="grey", xlab="Depth of basins",
         ylab="Frequency")
    dev.off()
    
  }
  
  enbr <- rbind(enbr,data.frame(plant=plnam,SS=ss,EB=eb,ncol=ncol(ocmatf)))
  
  
  cluster = makeCluster(n.core)
  registerDoParallel(cluster)
  
  samp_ss <- foreach(i=1:nrow(pocm),
                     .packages = c("rELA","foreach"),.combine = "c")%dopar%{
                       Bi(pocm[i,],
                          hgestp,jestp)[[1]]
                     }
  
  stopCluster(cluster)
  
  detected_ss_all <- rbind(detected_ss_all,
                       cbind(plant=plnam,
                             SS=ss,
                             ncol=ncol(pocm)))
  SS_samp <- rbind(SS_samp,
                   cbind(ID=rownames(pocm),
                                 plant=plnam,SS=samp_ss,ncol=ncol(pocm)))
  
  
  
  detected_ss <- rbind(detected_ss,
                           cbind(plant=plnam,
                                 SS=ss[ss %in% names(table(samp_ss))[table(samp_ss)/length(samp_ss)>maj_th]],
                                 ncol=ncol(pocm)))
  
  print(sum(ss %in% names(table(samp_ss))[table(samp_ss)/length(samp_ss)>maj_th]))
}

  saveRDS(detected_ss,sprintf("%s/detected_SS_%s.rds",dir,fb))
  # gstbe <- gStability(sa, ocmatf, enmat=sp_info[rownames(ocmatf),-c(1,ncol(sp_info))],
  #                     th=0., threads=n.core)
  # 
  # saveRDS(gstbe,sprintf("%s/%s/gStability_%s.rds",save.dir,fb,plnam))
  # 
  maj_s <- unique(detected_ss[,"SS"])
  
  #SS taxa binary heatmap
  sscf <- t(sapply(maj_s,id2bin,ncol(ocmatf)))
  colnames(sscf) <- colnames(ocmatf)
  rownames(sscf) <- maj_s
  
  sscf2 <- unique(sscf[,colSums(sscf)>0])
  
  dss <- vegdist(sscf2, method='jaccard')
  hcd <- hclust(dss)
  
  dss_y <- vegdist(t(sscf2), method='jaccard')
  hcd_y <- hclust(dss_y)
  
  sscf3 <- sscf2[hcd$order,]
  rownames(sscf3) <- sprintf("F_B%s",1:nrow(sscf2))
  
  saveRDS(sscf3,sprintf("%s/%s/SScomposition_%s.rds",save.dir,fb,fb))
  
  SStab <- cbind(rownames(sscf2)[hcd$order],
                 rename_SS=sprintf("F_B%s",1:nrow(sscf2)),
                 color=color[-4][1:length(rownames(sscf2)[hcd$order])])
  
  rownames(SStab) <- SStab[,1]
  
  saveRDS(SStab, file=sprintf("%s/graph_param_SStab_full_%s_%s.rds",dir,fb,taxa))
  
  
  gsscf <- gather(cbind(rename_ss=sprintf("F_B%s",1:nrow(sscf2)),
                        ssid=rownames(sscf2)[hcd$order],
                        as.data.frame(sscf2[hcd$order,])),
                  key=Taxa,value=pres,-c(1,2))
  
  
  gsscf$ss <- factor(gsscf$rename_ss,levels=sprintf("F_B%s",1:nrow(sscf2)))
  gsscf$Taxa2 <- factor(gsscf$Taxa,levels=colnames(sscf2)[hcd_y$order])
  
  #all basin
  gotu <- ggplot(gsscf, aes(x=ss,y=Taxa2,fill=as.factor(pres)))+
    geom_tile(show.legend = F)+
    labs(fill='Presence',y=sprintf("Fungi (%s)",taxa),x='Basins of attraction')+
    theme_bw()+
    theme(aspect.ratio=40/length(unique(gsscf$ss)),
          axis.title = element_text(size=15),
          legend.text = element_markdown(size=12),
          axis.text.x = element_markdown(size=12,angle=45,hjust=1,vjust=1),
          axis.text.y = element_markdown(size=12))+
    scale_fill_manual(values=c("white","darkorange"))
  #gotu
  
  
  ggsave(sprintf("%s/SSOTUcomp_heatmap_full_All_%s_%s.pdf",dir,fb,taxa),
         plot=gotu,w=13,h=9)
  
  #each plant
  for(pl in 1:nrow(pltab)){#pl <- 4
    
    if(sum(pltab[pl,]==1)){
      plnam <- colnames(pltab)[which(pltab[pl,]==1)]
    }else{
      plnam <- colnames(sp_info)[ncol(sp_info)]
    }
    
    ss_pl <- unique(detected_ss[detected_ss[,"plant"] == plnam,"SS"])
  gotu <- ggplot(gsscf[gsscf$ssid %in% ss_pl,], aes(x=ss,
                            y=Taxa2,
                            fill=as.factor(pres)))+
    geom_tile(show.legend = F)+
    labs(fill='Presence',y=sprintf("Fungi (%s)",taxa),x='Basin of attraction')+
    theme_bw()+
    theme(aspect.ratio=6/sum(detected_ss[,"plant"] == plnam),
      axis.title = element_text(size=15),
          legend.text = element_markdown(size=12),
          axis.text.x = element_markdown(size=12,angle=45,hjust=1,vjust=1),
          axis.text.y = element_markdown(size=12))+
    scale_fill_manual(values=c("white","darkorange"))
  #gotu
  
  
  ggsave(sprintf("%s/SSOTUcomp_heatmap_full_%s_%s_%s.pdf",dir,plnam,fb,taxa),
         plot=gotu,w=6,h=4)
  saveRDS(gsscf[gsscf$ssid %in% ss_pl,],
          sprintf("%s/SSOTUcomp_heatmap_full_%s_%s_%s.rds",dir,plnam,fb,taxa))
  
  pdf(sprintf("%s/DG_%s_%s.pdf",dir,plnam,fb),
      width=6,height=4)
  showDG_mod(ela[[pl]],ocmatf,SS_colmat= as.data.frame(matrix(SStab[ss_pl,],ncol=3,
                                                              dimnames = list(ss_pl,c("SS","rename_SS","color")))),
             label=sprintf("%s %s",fb,plnam),
             na.color="black",
             minor.color = "gray50",
             annot_adj=c(0.5, 2.00))
  dev.off()
  
  }
  ####

#######################################################################
fb <-"Prokaryota"
dir <- sprintf("%s/%s",save.dir,fb)
dir.create(dir)


ocmatf <- ocmat[[fb]]

dim(ocmatf)
saveRDS(ocmatf,sprintf("%s/input_ocmat_%s.rds",dir,fb))


sa <- runSA(ocmat=as.matrix(ocmatf),
            enmat = sp_info[rownames(ocmatf),-c(1,ncol(sp_info))],
            qth=qth, rep=1280, threads=n.core)

saveRDS(sa, file=sprintf("%s/runSA_ST_full_%s_%s.rds",dir,fb,taxa))

pltab <- unique(sp_info[rownames(ocmatf),-c(1,ncol(sp_info))])

ela <- list(NULL)
enbr <- NULL
SS_samp <- NULL
detected_ss <- NULL
detected_ss_all <- NULL

for(pl in 1:nrow(pltab)){#pl <- 4
  
  if(sum(pltab[pl,]==1)){
    plnam <- colnames(pltab)[which(pltab[pl,]==1)]
  }else{
    plnam <- colnames(sp_info)[ncol(sp_info)]
  }
  
  pocm <- ocmatf[which(rownames(ocmatf) %in% info[info$plant==plnam,"ID"]),]
  
  sa2p <- sa2params(sa,as.numeric(pltab[pl,]))
  
  hgestp <- sa2p[[4]]
  jestp <- sa2p[[2]]
  
  
  png(sprintf("%s/h_histgram_%s_%s.png",dir,fb,plnam),width=10,height=10,units="in",res=300)
  hist(sa2p[[1]], breaks=100, col="grey", main="intercept", xlab="h", ylab="Frequency")
  dev.off()
  
  png(sprintf("%s/J_histgram_%s_%s.png",dir,fb,plnam),width=10,height=10,units="in",res=300)
  hist(jestp, breaks=100, col="grey", main="Associations", xlab="J", ylab="Frequency")
  dev.off()
  
  ela[[pl]] <- ELA(sa, env=as.numeric(pltab[pl,]),
             SS.itr=SS.itr, FindingTip.itr=10000, # <- the number of steps for finding stable states and tipping points (basically no need to change)
             threads=n.core, reporting=TRUE)
  
  saveRDS(ela[[pl]], file=sprintf("%s/ELA_full_%s_%s.rds",dir,fb,plnam))
  
  ss <- ela[[pl]][[1]]
  eb <- c()
  
  for(i in 1:length(ss)){
    e.stable.np <- ela[[pl]][[2]][i]
    stablestates.np <- ela[[pl]][[1]]
    tippingen.np <- ela[[pl]][[4]]
    e.tipping.np <- min(unlist(c(as.data.frame(tippingen.np)[stablestates.np==ss[i],],
                                 as.data.frame(tippingen.np)[,stablestates.np==ss[i]])))
    eb[i] <- e.tipping.np - e.stable.np
  };rm(i)
  
  
  if(length(eb)>1){
    
    pdf(sprintf("%s/Bdepth_histgram_%s_%s.pdf",dir,fb,plnam),width=100,height=100)
    hist(eb/max(eb), breaks=100, col="grey", xlab="Depth of basins",
         ylab="Frequency")
    dev.off()
    
  }
  
  enbr <- rbind(enbr,data.frame(plant=plnam,SS=ss,EB=eb,ncol=ncol(ocmatf)))
  
  
  cluster = makeCluster(n.core)
  registerDoParallel(cluster)
  
  samp_ss <- foreach(i=1:nrow(pocm),
                     .packages = c("rELA","foreach"),.combine = "c")%dopar%{
                       Bi(pocm[i,],
                          hgestp,jestp)[[1]]
                     }
  
  stopCluster(cluster)
  
  
  detected_ss_all <- rbind(detected_ss_all,
                           cbind(plant=plnam,
                                 SS=ss,
                                 ncol=ncol(pocm)))
  SS_samp <- rbind(SS_samp,
                   cbind(ID=rownames(pocm),
                         plant=plnam,SS=samp_ss,ncol=ncol(pocm)))
  
  
  
  detected_ss <- rbind(detected_ss,
                       cbind(plant=plnam,
                             SS=ss[ss %in% names(table(samp_ss))[table(samp_ss)/length(samp_ss)>maj_th]],
                             ncol=ncol(pocm)))
  
  print(sum(ss %in% names(table(samp_ss))[table(samp_ss)/length(samp_ss)>maj_th]))
}

saveRDS(detected_ss,sprintf("%s/detected_SS_%s.rds",dir,fb))
# gstbe <- gStability(sa, ocmatf, enmat=sp_info[rownames(ocmatf),-c(1,ncol(sp_info))],
#                     th=0., threads=n.core)
# 
# saveRDS(gstbe,sprintf("%s/%s/gStability_%s.rds",save.dir,fb,plnam))
# 
maj_s <- unique(detected_ss[,"SS"])

#SS taxa binary heatmap
sscf <- t(sapply(maj_s,id2bin,ncol(ocmatf)))
colnames(sscf) <- colnames(ocmatf)
rownames(sscf) <- maj_s

sscf2 <- unique(sscf[,colSums(sscf)>0])

dss <- vegdist(sscf2, method='jaccard')
hcd <- hclust(dss)

dss_y <- vegdist(t(sscf2), method='jaccard')
hcd_y <- hclust(dss_y)

sscf3 <- sscf2[hcd$order,]
rownames(sscf3) <- sprintf("P_B%s",1:nrow(sscf2))

saveRDS(sscf3,sprintf("%s/%s/SScomposition_%s.rds",save.dir,fb,fb))


SStab <- cbind(rownames(sscf2)[hcd$order],
               rename_SS=sprintf("P_B%s",1:nrow(sscf2)),
               color=color[-4][1:length(rownames(sscf2)[hcd$order])])

rownames(SStab) <- SStab[,1]

saveRDS(SStab, file=sprintf("%s/graph_param_SStab_full_%s_%s.rds",dir,fb,taxa))


gsscf <- gather(cbind(rename_ss=sprintf("P_B%s",1:nrow(sscf2)),
                      ssid=rownames(sscf2)[hcd$order],
                      as.data.frame(sscf2[hcd$order,])),
                key=Taxa,value=pres,-c(1,2))


gsscf$ss <- factor(gsscf$rename_ss,levels=sprintf("P_B%s",1:nrow(sscf2)))
gsscf$Taxa2 <- factor(gsscf$Taxa,levels=colnames(sscf2)[hcd_y$order])

#all basin
gotu <- ggplot(gsscf, aes(x=ss,y=Taxa2,fill=as.factor(pres)))+
  geom_tile(show.legend = F)+
  labs(fill='Presence',y=sprintf("Prokaryotes (%s)",taxa),x='Basin of attraction')+
  theme_bw()+
  theme(aspect.ratio=40/length(unique(gsscf$ss)),
        axis.title = element_text(size=15),
        legend.text = element_markdown(size=12),
        axis.text.x = element_markdown(size=12,angle=45,hjust=1,vjust=1),
        axis.text.y = element_markdown(size=12))+
  scale_fill_manual(values=c("white","darkgreen"))
#gotu


ggsave(sprintf("%s/SSOTUcomp_heatmap_full_All_%s_%s.pdf",dir,fb,taxa),
       plot=gotu,w=13,h=9)

##each plant
for(pl in 1:nrow(pltab)){#pl <- 4
  
  if(sum(pltab[pl,]==1)){
    plnam <- colnames(pltab)[which(pltab[pl,]==1)]
  }else{
    plnam <- colnames(sp_info)[ncol(sp_info)]
  }
  
  ss_pl <- unique(detected_ss[detected_ss[,"plant"] == plnam,"SS"])
  gotu <- ggplot(gsscf[gsscf$ssid %in% ss_pl,], aes(x=ss,
                                                    y=Taxa2,
                                                    fill=as.factor(pres)))+
    geom_tile(show.legend = F)+
    labs(fill='Presence',y=sprintf("Prokaryotes (%s)",taxa),x='Basins of attraction')+
    theme_bw()+
    theme(aspect.ratio=40/sum(detected_ss[,"plant"] == plnam),
          axis.title = element_text(size=15),
          legend.text = element_markdown(size=12),
          axis.text.x = element_markdown(size=12,angle=45,hjust=1,vjust=1),
          axis.text.y = element_markdown(size=12))+
    scale_fill_manual(values=c("white","darkgreen"))
  #gotu
  
  
  ggsave(sprintf("%s/SSOTUcomp_heatmap_full_%s_%s_%s.pdf",dir,plnam,fb,taxa),
         plot=gotu,w=12,h=8)
  saveRDS(gsscf[gsscf$ssid %in% ss_pl,],
          sprintf("%s/SSOTUcomp_heatmap_full_%s_%s_%s.rds",dir,plnam,fb,taxa))
  
  
  pdf(sprintf("%s/DG_%s_%s.pdf",dir,plnam,fb),
      width=6,height=4)
  showDG_mod(ela[[pl]],ocmatf,SS_colmat= as.data.frame(matrix(SStab[ss_pl,],ncol=3,
                                                              dimnames = list(ss_pl,c("SS","rename_SS","color")))),
             label=sprintf("%s %s",fb,plnam),
             na.color="black",
             minor.color = "gray50",
             annot_adj=c(0.5, 2.00))
  dev.off()
}
####
