
##########################################################################
#function

#setwd("/Users/ngch/Desktop/Now_using/data/sugaFB_ELA/revise_251212/output_spacom/result_for_revise/analysis_in_local_computer")
set.seed(1234)

####################

#install.packages("mgcv")
library(ggplot2)
library(ggtext)
library(ggrepel)
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
library(vegan)
library(ggalluvial)
library(dplyr)
library(ComplexHeatmap)

library(circlize)

#######

save.dir <- "Output/02_06_07_02_graphics_ELA_genuslevel"
dir.create(save.dir)

#########
dir_02_06_07 <- "Output_supercomputer/02_06_07_ELA_genuslevel"

#########
tx <- list(Fungi=readRDS("Base_data/Fungi/taxa_list_mod.rds"),
           Prokaryota=readRDS("Base_data/Bacteria/240404_Megre/merged/taxonomy_list_dada2.rds"))

phy <- c("Phylum","Class","Order","Family","Genus")
#alpha_ss_col <- 0.8
ss_initial <- "FG"
fb <- "Fungi"

sa <- readRDS(sprintf("%s/%s/runSA_ST_full_Fungi_Genus.rds",dir_02_06_07,fb))
df_f <- readRDS(sprintf("%s/%s/detected_SS_Fungi.rds",dir_02_06_07,fb))
tx2 <- tx[[fb]]

com_tx <- unique(tx2[tx2[,"Genus"] %in% rownames(sa[[1]]),phy])
rownames(com_tx) <- com_tx[,"Genus"]
write.csv(com_tx[order(com_tx[,"Phylum"],com_tx[,"Family"]),],sprintf("%s/sscomp_taxalist_%s.csv",save.dir,fb))

ssnams <- unique(df_f[,2])


mss <- T(apply(unique(df_f[,-1]),1,function(x){
  id2bin(x[1],as.numeric(x[2]))
}))

dimnames(mss) <- list(ssnams,rownames(sa[[1]]))

ssc_f <- mss[,colSums(mss)>0]

pss_f <- as.data.frame(df_f[,1:2])
wid_f <- spread(cbind(unique(pss_f[,c("SS","plant")]),val=1),plant,val)
wid_f[is.na(wid_f)] <- 0
rownames(wid_f) <- wid_f[,1]

jac_dist <- function(x) {
  vegdist(x, method = "jaccard")
}

saveRDS(list(host_plant=wid_f,genus_comp=ssc_f),
        sprintf("%s/genuscomp_for_heatmap_%s.rds",save.dir,fb))
plant_order <- c("Pinus","Larix","Betula","Populus","Acer","Juglans")
col_fun = colorRamp2(c(0, 1), c("white", "white"))

plss <- Heatmap(wid_f[,plant_order], row_names_side = "left", column_names_side = "bottom", 
                col = col_fun,
                cell_fun = function(j, i,x,y, width, height, fill) {
                  if(wid_f[rownames(ssc_f),plant_order][i,j] == 1){
                    grid.circle(x = x, y = y, r = min(unit.c(width, height))/2.2,
                                gp = gpar(fill = "darkred"))
                  }
                },
                show_heatmap_legend = FALSE,
                cluster_columns = FALSE,
                show_row_dend = FALSE,
                show_column_dend = FALSE,
                rect_gp = gpar(col = "black"), 
                #row_names_gp = gpar(cex=0.9, fontface = "bold"), 
                column_names_gp = gpar(cex=0.9, fontface = "bold.italic"), 
                row_dend_width = unit(2, "cm"),
                column_dend_height = unit(2, "cm"), 
                column_title = "Host plant",
                column_names_rot = 45, #row_title = "Basins of attraction", 
                #clustering_distance_rows = jac_dist,
                clustering_distance_columns = jac_dist, 
                #clustering_method_rows = "ward.D2",
                clustering_method_columns = "ward.D2")





ht_ord <- draw(Heatmap(ssc_f,clustering_distance_rows = jac_dist,
                       clustering_distance_columns = jac_dist, 
                       clustering_method_rows = "ward.D2",
                       clustering_method_columns = "ward.D2"))

rord <- row_order(ht_ord)
cord <- column_order(ht_ord)

##set new ids
new_ssid <- sprintf("%s_B%s",ss_initial,1:length(ssnams))

###########manual###################
sscolor <- c(rep("darkorange",length(new_ssid)))

####################################

names(sscolor) <- new_ssid

names(new_ssid) <- ssnams[rord]

rownames(ssc_f) <- new_ssid[rownames(ssc_f)]
####
ssc_f2 <- ssc_f[rord,cord]
cust_col <- foreach(i = 1:nrow(ssc_f2), .combine = rbind) %do% {
  ifelse(ssc_f2[i, ]==1,sscolor[rownames(ssc_f2)[i]],"white")
} 
rownames(cust_col) <- rownames(ssc_f2)


sscomp <- Heatmap(ssc_f, row_names_side = "left", column_names_side = "bottom", 
                  col = col_fun,
                  show_column_dend = FALSE,
                  
                  show_heatmap_legend = FALSE,
                  row_dend_side = "left", 
                  rect_gp = gpar(fill = as.vector(cust_col),
                                 col = "black"), 
                  row_names_gp = gpar(cex=0.9, fontface = "bold"), 
                  column_names_gp = gpar(cex=0.9, fontface = "bold.italic"), 
                  row_dend_width = unit(2, "cm"),
                  column_dend_height = unit(2, "cm"), 
                  column_title = sprintf("%s (Genus)",ifelse(fb=="Fungi",fb,"Prokaryote")),
                  column_names_rot = 45, row_title = "Basins bottoms", 
                  clustering_distance_rows = jac_dist,
                  clustering_distance_columns = jac_dist, 
                  clustering_method_rows = "ward.D2",
                  clustering_method_columns = "ward.D2")

mh <- draw(sscomp+plss)

#check row order
all(rord==row_order(mh))

pdf(sprintf("%s/SSheatmap_%s.pdf",save.dir,fb),width = 8,h=4.5)
mh

dev.off()
###################
ss_initial <- "PG"
fb <- "Prokaryota"

sa <- readRDS(sprintf("%s/%s/runSA_ST_full_%s_Genus.rds",dir_02_06_07,fb,fb))
df_f <- readRDS(sprintf("%s/%s/detected_SS_%s.rds",dir_02_06_07,fb,fb))
tx2 <- tx[[fb]]

com_tx <- unique(tx2[tx2[,"Genus"] %in% rownames(sa[[1]]),phy])
rownames(com_tx) <- com_tx[,"Genus"]
write.csv(com_tx[order(com_tx[,"Phylum"],com_tx[,"Family"]),],sprintf("%s/sscomp_taxalist_%s.csv",save.dir,fb))

ssnams <- unique(df_f[,2])

mss <- T(apply(unique(df_f[,-1]),1,function(x){
  id2bin(x[1],as.numeric(x[2]))
}))

dimnames(mss) <- list(ssnams,rownames(sa[[1]]))

ssc_f <- mss[,colSums(mss)>0]

pss_f <- as.data.frame(df_f[,1:2])
wid_f <- spread(cbind(unique(pss_f[,c("SS","plant")]),val=1),plant,val)
wid_f[is.na(wid_f)] <- 0
rownames(wid_f) <- wid_f[,1]

jac_dist <- function(x) {
  vegdist(x, method = "jaccard")
}

saveRDS(list(host_plant=wid_f,genus_comp=ssc_f),
        sprintf("%s/genuscomp_for_heatmap_%s.rds",save.dir,fb))

plant_order <- c("Pinus","Larix","Betula","Populus","Acer","Juglans")
col_fun = colorRamp2(c(0, 1), c("white", "white"))

plss <- Heatmap(wid_f[,plant_order], row_names_side = "left", column_names_side = "bottom", 
                col = col_fun,
                cell_fun = function(j, i,x,y, width, height, fill) {
                  if(wid_f[rownames(ssc_f),plant_order][i,j] == 1){
                    grid.circle(x = x, y = y, r = min(unit.c(width, height))/2.1,
                                gp = gpar(fill = "darkred"))
                  }
                },
                show_heatmap_legend = FALSE,
                cluster_columns = FALSE,
                show_row_dend = FALSE,
                show_column_dend = FALSE,
                rect_gp = gpar(col = "black"), 
                #row_names_gp = gpar(cex=0.9, fontface = "bold"), 
                column_names_gp = gpar(cex=0.9, fontface = "bold.italic"), 
                row_dend_width = unit(2, "cm"),
                column_dend_height = unit(2, "cm"), 
                column_title = "Host plant",
                column_names_rot = 45, #row_title = "Basins of attraction", 
                #clustering_distance_rows = jac_dist,
                clustering_distance_columns = jac_dist, 
                #clustering_method_rows = "ward.D2",
                clustering_method_columns = "ward.D2")





ht_ord <- draw(Heatmap(ssc_f,clustering_distance_rows = jac_dist,
                       clustering_distance_columns = jac_dist, 
                       clustering_method_rows = "ward.D2",
                       clustering_method_columns = "ward.D2"))

rord <- row_order(ht_ord)
cord <- column_order(ht_ord)

##set new ids
new_ssid <- sprintf("%s_B%s",ss_initial,1:length(ssnams))

###########manual###################
sscolor <- c(rep("darkgreen",length(new_ssid)))

####################################

names(sscolor) <- new_ssid

names(new_ssid) <- ssnams[rord]

rownames(ssc_f) <- new_ssid[rownames(ssc_f)]
####
ssc_f2 <- ssc_f[rord,cord]
cust_col <- foreach(i = 1:nrow(ssc_f2), .combine = rbind) %do% {
  ifelse(ssc_f2[i, ]==1,sscolor[rownames(ssc_f2)[i]],"white")
} 
rownames(cust_col) <- rownames(ssc_f2)


sscomp <- Heatmap(ssc_f, row_names_side = "left", column_names_side = "bottom", 
                  col = col_fun,
                  show_column_dend = FALSE,
                  
                  show_heatmap_legend = FALSE,
                  row_dend_side = "left", 
                  rect_gp = gpar(fill = as.vector(cust_col),
                                 col = "black"), 
                  row_names_gp = gpar(cex=0.9, fontface = "bold"), 
                  column_names_gp = gpar(cex=0.9, fontface = "bold.italic"), 
                  row_dend_width = unit(2, "cm"),
                  column_dend_height = unit(2, "cm"), 
                  column_title = sprintf("%s (Genus)",ifelse(fb=="Fungi",fb,"Prokaryote")),
                  column_names_rot = 45, row_title = "Basins bottoms", 
                  clustering_distance_rows = jac_dist,
                  clustering_distance_columns = jac_dist, 
                  clustering_method_rows = "ward.D2",
                  clustering_method_columns = "ward.D2")

mh <- draw(sscomp+plss)

#check row order
all(rord==row_order(mh))

pdf(sprintf("%s/SSheatmap_%s.pdf",save.dir,fb),width = 8,h=4.5)
mh

dev.off()

