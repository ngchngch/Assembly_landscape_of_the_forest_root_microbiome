set.seed(1234)

# install.packages("devtools")
# library(devtools)
# devtools::install_github("zdk123/SpiecEasi",force=TRUE)

####################

#install.packages("mgcv")
library(ggplot2)
library(ggtext)
library(parallel)
library('RColorBrewer')
#library('scatterpie')
library(bipartite)
library("rELA")
library('tidyverse')
library("doParallel")
library(foreach)
#library(renv)
library("Rcpp")
library("RcppArmadillo")
library(igraph)

library(ggnetwork)

#######

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

#######








ELA_prep_dir <- "Output_supercomputer/02_01_ELA_prep_abundance_threshold"

dir_02_05 <- "Output_supercomputer/02_05_summary_taxa_select"

#########################################################################
save.dir <- "Output/02_08_hostpreference_Family"
dir.create(save.dir)

########################################################################
#read original functions
source("packages/01_1_function.R")

########################################################################
tx_f <- readRDS("Base_data/Fungi/taxa_list_mod.rds")
tx_p <- as.data.frame(readRDS("Base_data/OTU_basedata_set/taxonomy_list_dada2.rds"))

info <- readRDS(sprintf("%s/comp_sample_info_plant2.rds",ELA_prep_dir))


tb_gns <- readRDS(sprintf("%s/NoCLR_ExpVar_clrRA_target_taxa.rds",ELA_prep_dir))
tb_g <- readRDS(sprintf("%s/ExpVar_clrRA_target_taxa.rds",ELA_prep_dir))


oc_f2 <- readRDS(sprintf("%s/ELA_input_ocmat_Fungi.rds",dir_02_05))

oc_p2 <- readRDS(sprintf("%s/ELA_input_ocmat_Prokaryote.rds",dir_02_05))
    

####################################################
n.core <- 8

rand <- 10000

##Fungi
#original

or_f <- Taxa.mat(t(oc_f2),info,"plant")
or_f2 <- or_f

#dprime
or_dpr_p <- dfun(t(or_f2))$dprime
or_dpr_f <- dfun(or_f2)$dprime

#randamization

cluster = makeCluster(n.core)
registerDoParallel(cluster)

ror_f <- foreach(r = 1:rand,.packages = "bipartite")%dopar%{
  set.seed(r)
  roc_f <- blockSample(oc_f2,info[rownames(oc_f2),"site"],rownames(oc_f2))$matrix
  mat <-  Taxa.mat(t(roc_f),info,"plant")[rownames(or_f2),]
  dpr_f <- dfun(mat)$dprime
  dpr_p <- dfun(t(mat))$dprime
  return(list(mat=mat,dpr_f=dpr_f,dpr_p=dpr_p))
}

rand_dpr_f <- foreach(r = 1:rand,.combine = "rbind")%dopar%{
  return(ror_f[[r]][["dpr_f"]])
}

z_dpr_f <- (or_dpr_f-apply(rand_dpr_f,2,mean))/apply(rand_dpr_f,2,sd)

p_dpr_f <- p.adjust(foreach(r = 1:length(or_dpr_f),.combine = "c")%do%{
  p <- sum(rand_dpr_f[,r]>=or_dpr_f[r])/rand
  ifelse(p!=0,p,0.0000999999)},method = "BH")
  

rand_dpr_p <- foreach(r = 1:rand,.combine = "rbind")%dopar%{
  return(ror_f[[r]][["dpr_p"]])
}

z_dpr_p <- (or_dpr_p-apply(rand_dpr_p,2,mean))/apply(rand_dpr_p,2,sd)

p_dpr_p <- p.adjust(foreach(r = 1:length(or_dpr_p),.combine = "c")%do%{
  p <- sum(rand_dpr_p[,r]>=or_dpr_p[r])/rand
  ifelse(p!=0,p,0.0000999999)},method = "BH")

z_mat <- matrix(NA,nrow = nrow(or_f2),ncol = ncol(or_f2),dimnames = dimnames(or_f2))
p_mat <- matrix(NA,nrow = nrow(or_f2),ncol = ncol(or_f2),dimnames = dimnames(or_f2))

  for(row in 1:nrow(or_f2)){
    for(col in 1:ncol(or_f2)){
      cat(row,col,"\n")
      #row <- 1;col <- 1
      rand_val <- foreach(r = 1:rand,.combine = "c")%dopar%{
      ror_f[[r]][["mat"]][row,col]
      }
      z_mat[row,col] <- (or_f2[row,col]-mean(rand_val))/sd(rand_val)
      
      pm <- ifelse(is.na(z_mat[row,col]),0,z_mat[row,col])
      if(pm>0){
        p_mat[row,col] <- sum(or_f2[row,col]<=rand_val)/rand
      }else{
        p_mat[row,col] <- sum(or_f2[row,col]>=rand_val)/rand
      }
      
  }
}

stopCluster(cluster)

saveRDS(z_mat,sprintf("%s/z_mat_Fungi.rds",save.dir))
saveRDS(p_mat,sprintf("%s/p_mat_Fungi.rds",save.dir))

saveRDS(list(dprime_plant=z_dpr_p,
              dprime_fungus=z_dpr_f,
              p_dprime_plant=p_dpr_p,
              p_dprime_fungus=p_dpr_f),
         sprintf("%s/dprime_Fungi.rds",save.dir))

p_mat_BH <- matrix(p.adjust(as.vector(p_mat),method = "BH"),nrow = nrow(p_mat),ncol = ncol(p_mat),
                   dimnames = dimnames(p_mat))
col_dpl <- colorRamp2(c(min(z_dpr_p,na.rm=T),
                        max(z_dpr_p,na.rm=T)),
                      c("white", "chartreuse"))

ha1 <- HeatmapAnnotation(
  show_annotation_name = F, 
  p_value = anno_text(sapply(p_dpr_p,function(x)sig(x,star=TRUE,bin=FALSE)),which = "column",
                      just="top",location=0.5,
                      gp=gpar(fontsize=17)),
  dprime_host=z_dpr_p,
  col=list(dprime_host=col_dpl),
  gp = gpar(col = "black"))


col_dpf <- colorRamp2(c(min(z_dpr_f,na.rm=T),
                        max(z_dpr_f,na.rm=T)),
                      c("white", "deeppink2"))

rha <- rowAnnotation(p_value = anno_text(sapply(p_dpr_f,function(x)sig(x,star=TRUE,bin=FALSE)),
                                         which = "row",just="top",location=0.5,
                                         gp=gpar(fontsize=17)),
                     dprime_symbiont=z_dpr_f,
                     show_annotation_name = F,
                     col=list(dprime_symbiont=col_dpf),
                     
                     gp = gpar(col = "black"))


col1 = colorRamp2(c(min(z_mat,na.rm=T), 0, max(z_mat,na.rm=T)), c("blue", "white", "red"))

pdf(sprintf("%s/Fungi_2DP.pdf",save.dir),width=9,height=11)

Heatmap(z_mat, 
        col=col1,show_row_names = TRUE,
        show_column_names = TRUE,
        show_column_dend = FALSE, show_row_dend = FALSE,
        row_names_gp = gpar(fontsize = 14),
        column_names_gp = gpar(fontsize = 14,fontface = "italic"),
        #top_annotation = tha,
        bottom_annotation = ha1,right_annotation = rha,
        #row_labels = gt_render(d),
        cluster_rows = TRUE, cluster_columns = TRUE,
        clustering_method_columns = "ward.D2",
        clustering_method_rows = "ward.D2",
        na_col="gray90",
        rect_gp = gpar(col = "black"),
        heatmap_legend_param = list(title="Two dimensional preference (2DP)"),
        cell_fun = function(j, i, x, y, w, h, fill) {
          if(p_mat_BH[i, j] < 0.0005) {
            grid.text("***", x, y,rot=0,vjust = 0.8,
                      gp=gpar(col="white", fontsize=18))
          } else {if(p_mat_BH[i, j] < 0.005) {
            grid.text("**", x, y,rot=0,vjust = 0.8,
                      gp=gpar(col="white", fontsize=18))
          }else{if(p_mat_BH[i, j] < 0.025) {
            grid.text("*", x, y,rot=0,vjust = 0.8,
                      gp=gpar(col="white", fontsize=18))
          }
          }
          }
        }
)

dev.off()

##Prokaryote
#original

or_f <- Taxa.mat(t(oc_p2),info,"plant")
or_f2 <- or_f

#dprime
or_dpr_p <- dfun(t(or_f2))$dprime
or_dpr_f <- dfun(or_f2)$dprime

#randamization

cluster = makeCluster(n.core)
registerDoParallel(cluster)

ror_f <- foreach(r = 1:rand,.packages = "bipartite")%dopar%{
  set.seed(r)
  roc_f <- blockSample(oc_p2,info[rownames(oc_p2),"site"],rownames(oc_p2))$matrix
  mat <-  Taxa.mat(t(roc_f),info,"plant")[rownames(or_f2),]
  dpr_f <- dfun(mat)$dprime
  dpr_p <- dfun(t(mat))$dprime
  return(list(mat=mat,dpr_f=dpr_f,dpr_p=dpr_p))
}

rand_dpr_f <- foreach(r = 1:rand,.combine = "rbind")%dopar%{
  return(ror_f[[r]][["dpr_f"]])
}

z_dpr_f <- (or_dpr_f-apply(rand_dpr_f,2,mean))/apply(rand_dpr_f,2,sd)

p_dpr_f <- p.adjust(foreach(r = 1:length(or_dpr_f),.combine = "c")%do%{
  p <- sum(rand_dpr_f[,r]>=or_dpr_f[r])/rand
  ifelse(p!=0,p,0.0000999999)},method = "BH")


rand_dpr_p <- foreach(r = 1:rand,.combine = "rbind")%dopar%{
  return(ror_f[[r]][["dpr_p"]])
}

z_dpr_p <- (or_dpr_p-apply(rand_dpr_p,2,mean))/apply(rand_dpr_p,2,sd)

p_dpr_p <- p.adjust(foreach(r = 1:length(or_dpr_p),.combine = "c")%do%{
  p <- sum(rand_dpr_p[,r]>=or_dpr_p[r])/rand
  ifelse(p!=0,p,0.0000999999)},method = "BH")

z_mat <- matrix(NA,nrow = nrow(or_f2),ncol = ncol(or_f2),dimnames = dimnames(or_f2))
p_mat <- matrix(NA,nrow = nrow(or_f2),ncol = ncol(or_f2),dimnames = dimnames(or_f2))

for(row in 1:nrow(or_f2)){
  for(col in 1:ncol(or_f2)){
    cat(row,col,"\n")
    #row <- 1;col <- 1
    rand_val <- foreach(r = 1:rand,.combine = "c")%dopar%{
      ror_f[[r]][["mat"]][row,col]
    }
    z_mat[row,col] <- (or_f2[row,col]-mean(rand_val))/sd(rand_val)
    
    pm <- ifelse(is.na(z_mat[row,col]),0,z_mat[row,col])
    if(pm>0){
      p_mat[row,col] <- sum(or_f2[row,col]<=rand_val)/rand
    }else{
      p_mat[row,col] <- sum(or_f2[row,col]>=rand_val)/rand
    }
    
  }
}

stopCluster(cluster)

saveRDS(z_mat,sprintf("%s/z_mat_Prokaryote.rds",save.dir))
saveRDS(p_mat,sprintf("%s/p_mat_Prokaryote.rds",save.dir))
write.csv(z_mat,sprintf("%s/z_mat_Prokaryote.csv",save.dir))
write.csv(p_mat,sprintf("%s/p_mat_Prokaryote.csv",save.dir))

p_mat_BH <- matrix(p.adjust(as.vector(p_mat),method = "BH"),nrow = nrow(p_mat),ncol = ncol(p_mat),
                   dimnames = dimnames(p_mat))
col_dpl <- colorRamp2(c(min(z_dpr_p,na.rm=T),
                        max(z_dpr_p,na.rm=T)),
                      c("white", "chartreuse"))

ha1 <- HeatmapAnnotation(
  show_annotation_name = F, 
  p_value = anno_text(sapply(p_dpr_p,function(x)sig(x,star=TRUE,bin=FALSE)),which = "column",
                      just="top",location=0.5,
                      gp=gpar(fontsize=17)),
  dprime_host=z_dpr_p,
  col=list(dprime_host=col_dpl),
  gp = gpar(col = "black"))


col_dpf <- colorRamp2(c(min(z_dpr_f,na.rm=T),
                        max(z_dpr_f,na.rm=T)),
                      c("white", "deeppink2"))

rha <- rowAnnotation(p_value = anno_text(sapply(p_dpr_f,function(x)sig(x,star=TRUE,bin=FALSE)),
                                         which = "row",just="top",location=0.5,
                                         gp=gpar(fontsize=17)),
                     dprime_symbiont=z_dpr_f,
                     show_annotation_name = F,
                     col=list(dprime_symbiont=col_dpf),
                     
                     gp = gpar(col = "black"))


col1 = colorRamp2(c(min(z_mat,na.rm=T), 0, max(z_mat,na.rm=T)), c("blue", "white", "red"))

pdf(sprintf("%s/Prokaryote_2DP.pdf",save.dir),width=9,height=11)

Heatmap(z_mat, 
        col=col1,show_row_names = TRUE,
        show_column_names = TRUE,
        show_column_dend = FALSE, show_row_dend = FALSE,
        row_names_gp = gpar(fontsize = 14),
        column_names_gp = gpar(fontsize = 14,fontface = "italic"),
        #top_annotation = tha,
        bottom_annotation = ha1,right_annotation = rha,
        #row_labels = gt_render(d),
        cluster_rows = TRUE, cluster_columns = TRUE,
        clustering_method_columns = "ward.D2",
        clustering_method_rows = "ward.D2",
        na_col="gray90",
        rect_gp = gpar(col = "black"),
        heatmap_legend_param = list(title="Two dimensional preference (2DP)"),
        cell_fun = function(j, i, x, y, w, h, fill) {
          if(p_mat_BH[i, j] < 0.0005) {
            grid.text("***", x, y,rot=0,vjust = 0.8,
                      gp=gpar(col="white", fontsize=18))
          } else {if(p_mat_BH[i, j] < 0.005) {
            grid.text("**", x, y,rot=0,vjust = 0.8,
                      gp=gpar(col="white", fontsize=18))
          }else{if(p_mat_BH[i, j] < 0.025) {
            grid.text("*", x, y,rot=0,vjust = 0.8,
                      gp=gpar(col="white", fontsize=18))
          }
          }
          }
        }
)

dev.off()

