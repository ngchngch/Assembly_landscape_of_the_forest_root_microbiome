
##########################################################################
set.seed(1234)

library(ggplot2)
library(ggstar)
library(parallel)
library(purrr)
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


########################################################################
#read original functions
#source("packages/01_1_function.R")


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
    SS_colmat,na.color="black",minor.color="gray90",fontsize=8){
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
    xlab("Basins & Tipping points") +
    ylab("Energy")+
    geom_text_repel(hjust=annot_adj[1], vjust=annot_adj[2],size=fontsize,
              aes(fontface=2, label=point_str)) +
    theme_bw() +
    theme(axis.title = element_text(size=18),
          axis.text.y = element_text(size=13),
          axis.text.x = element_blank(),
          plot.title = element_markdown(size=20),
          aspect.ratio = 1,
          panel.background = element_rect(fill = "white", color = "black", size = 2))+
    scale_color_manual(values=col[which(levels(grobj_to_plot$point_str2) %in% unique(grobj_to_plot$point_str2))])+
    scale_shape_manual(values=c(15,15, 16)[which(levels(grobj_to_plot$shape) %in% unique(grobj_to_plot$shape))]) +
    scale_size_manual(values=c(5,3, 1)[which(levels(grobj_to_plot$size) %in% unique(grobj_to_plot$size))]) +
    coord_cartesian(ylim=c(min(grobj_to_plot$energy,na.rm = TRUE)-0.08*(max(grobj_to_plot$energy,na.rm = TRUE)-min(grobj_to_plot$energy,na.rm = TRUE)),NA),
                    xlim=c(min(grobj_to_plot$nodes2xposi,na.rm = TRUE)-0.06*(max(grobj_to_plot$nodes2xposi,na.rm = TRUE)-min(grobj_to_plot$nodes2xposi,na.rm = TRUE)),
                           max(grobj_to_plot$nodes2xposi,na.rm = TRUE)+0.06*(max(grobj_to_plot$nodes2xposi,na.rm = TRUE)-min(grobj_to_plot$nodes2xposi,na.rm = TRUE))))
  
  # pie
  #g <- g + geom_scatterpie(
  #  aes(x=nodes2xposi, y=energy),
  #  data=grobj_, cols=d_list)
  
  # plot
  g <- g + ggtitle(DG_sample)
  plot(g)
  return(NULL)}


####
save.dir <- "Output/03_12_DG_fullELA_recolor_250501"
dir.create(save.dir)

########################################################################
dir_03_10 <- "Output/03_10_graphics_states_flow_flow_spl_250501"
ELA_prep_dir <- "Output_supercomputer/02_01_ELA_prep_abundance_threshold"
dir_02_05 <- "Output_supercomputer/02_05_summary_taxa_select"
dir_02_06 <- "Output_supercomputer/02_06_ELA"

###
ocmat <- list(Fungi=readRDS(sprintf("%s/ELA_input_ocmat_Fungi.rds",dir_02_05)),
              Prokaryota=readRDS(sprintf("%s/ELA_input_ocmat_Prokaryote.rds",dir_02_05)))

sp_info <- readRDS(sprintf("%s/ELA_input_plant.rds",ELA_prep_dir))

###
sscolor <- readRDS(sprintf("%s/states_colvector.rds",dir_03_10))

ela_list <- list.files(pattern="^ELA_full",
                       dir_02_06,full.names = TRUE,recursive = TRUE)

sstab_list <- list.files(pattern="^graph_param_SStab",
                        dir_02_06,full.names = TRUE,recursive = TRUE)

dSS_list <- list.files(pattern="^detected_SS",
                        dir_02_06,full.names = TRUE,recursive = TRUE)
fb <-"Fungi"
dir <- sprintf("%s/%s",save.dir,fb)

dir.create(dir)

ocmatf <- ocmat[[fb]]
pltab <- unique(sp_info[rownames(ocmatf),-c(1,ncol(sp_info))])

for(pl in 1:nrow(pltab)){#pl <- 4
  
  if(sum(pltab[pl,]==1)){
    plnam <- colnames(pltab)[which(pltab[pl,]==1)]
  }else{
    plnam <- colnames(sp_info)[ncol(sp_info)]
  }
  
  ela <- readRDS(ela_list[grep(paste(fb,plnam,sep="_"),ela_list)])

  SStab <- readRDS(sstab_list[grep(fb,sstab_list)])
  
  detected_ss <- readRDS(dSS_list[grep(fb,dSS_list)])
  
  ss_pl <- unique(detected_ss[detected_ss[,"plant"] == plnam,"SS"])
  
  pdf(sprintf("%s/DG_%s_%s.pdf",dir,plnam,fb),
      width=6,height=4)
  
showDG_mod(ela,ocmatf,
           SS_colmat= as.data.frame(matrix(cbind(SStab[ss_pl,c(1,2)],
                                                 sscolor[SStab[ss_pl,2]]),ncol=3,
                                           dimnames = list(ss_pl,c("SS","rename_SS","color")))),
           label=sprintf("%s in *%s*",ifelse(fb=="Fungi","Fungal community","Prokaryotic community")
                                        ,plnam),
           na.color="black",
           minor.color = "gray50",
           annot_adj=c(0.5, 2.00))
dev.off()
}

fb <-"Prokaryota"
dir <- sprintf("%s/%s",save.dir,fb)
dir.create(dir)

ocmatf <- ocmat[[fb]]
pltab <- unique(sp_info[rownames(ocmatf),-c(1,ncol(sp_info))])

for(pl in 1:nrow(pltab)){#pl <- 4
  
  if(sum(pltab[pl,]==1)){
    plnam <- colnames(pltab)[which(pltab[pl,]==1)]
  }else{
    plnam <- colnames(sp_info)[ncol(sp_info)]
  }
  
  ela <- readRDS(ela_list[grep(paste(fb,plnam,sep="_"),ela_list)])
  
  SStab <- readRDS(sstab_list[grep(fb,sstab_list)])
  
  detected_ss <- readRDS(dSS_list[grep(fb,dSS_list)])
  
  ss_pl <- unique(detected_ss[detected_ss[,"plant"] == plnam,"SS"])
  
  pdf(sprintf("%s/DG_%s_%s.pdf",dir,plnam,fb),
      width=6,height=4)
  
  showDG_mod(ela,ocmatf,
             SS_colmat= as.data.frame(matrix(cbind(SStab[ss_pl,c(1,2)],
                                                   sscolor[SStab[ss_pl,2]]),ncol=3,
                                             dimnames = list(ss_pl,c("SS","rename_SS","color")))),
             label=sprintf("%s in *%s*",ifelse(fb=="Fungi","Fungal community","Prokaryotic community")
                           ,plnam),
             na.color="black",
             minor.color = "gray50",
             annot_adj=c(0.5, 2.00))
  dev.off()
}


