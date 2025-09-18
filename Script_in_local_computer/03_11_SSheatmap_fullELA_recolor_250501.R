
library(ComplexHeatmap)
library(cluster)

dir_03_10_02 <- "Output/03_10_graphics_states_flow_flow_spl_250508"

####
save.dir <- "Output/03_11_SSheatmap_fullELA_recolor_250501"
dir.create(save.dir)

###
sscolor <- readRDS(sprintf("%s/states_colvector.rds",dir_03_10_02))


ssc_f <- readRDS(sprintf("%s/Fungi/SScomposition_Fungi.rds",dir_02_06))


llf <- list.files(sprintf("%s/Fungi",dir_02_06),
                  pattern = "SSOTUcomp_heatmap_full",
           full.names = TRUE,recursive = TRUE)

pss_f <- do.call(rbind,lapply(llf[grep("rds",llf)],function(x){
  #x <- llf[1]
  nam <- strsplit(x,"_")[[1]][12]
  return(cbind(plant=nam,readRDS(x)))
  
}))

wid_f <- spread(cbind(unique(pss_f[,c("ss","plant")]),val=1),plant,val)
wid_f[is.na(wid_f)] <- 0
rownames(wid_f) <- wid_f[,1]

jac_dist <- function(x) {
  vegdist(x, method = "jaccard")
}

library(circlize)

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
        column_names_gp = gpar(cex=0.9, fontface = "bold"), 
        row_dend_width = unit(2, "cm"),
        column_dend_height = unit(2, "cm"), 
        column_title = "Fungi ( Family )",
        column_names_rot = 45, row_title = "Stable states", 
        clustering_distance_rows = jac_dist,
        clustering_distance_columns = jac_dist, 
        clustering_method_rows = "ward.D2",
        clustering_method_columns = "ward.D2")

mh <- draw(sscomp+plss)

#check row order
all(rord==row_order(mh))

pdf(sprintf("%s/SSheatmap_Fungi.pdf",save.dir),width = 8,h=4.5)
mh

dev.off()

###
ssc_f <- readRDS(sprintf("%s/Prokaryota/SScomposition_Prokaryota.rds",dir_02_06))


llf <- list.files(sprintf("%s/Prokaryota",dir_02_06),
                  pattern = "SSOTUcomp_heatmap_full",
                  full.names = TRUE,recursive = TRUE)

pss_f <- do.call(rbind,lapply(llf[grep("rds",llf)],function(x){
  #x <- llf[1]
  nam <- strsplit(x,"_")[[1]][12]
  return(cbind(plant=nam,readRDS(x)))
  
}))

wid_f <- spread(cbind(unique(pss_f[,c("ss","plant")]),val=1),plant,val)
wid_f[is.na(wid_f)] <- 0
rownames(wid_f) <- wid_f[,1]

jac_dist <- function(x) {
  vegdist(x, method = "jaccard")
}


ht_ord <- draw(Heatmap(ssc_f,clustering_distance_rows = jac_dist,
                       clustering_distance_columns = jac_dist, 
                       clustering_method_rows = "ward.D2",
                       clustering_method_columns = "ward.D2"))

rord <- row_order(ht_ord)
cord <- column_order(ht_ord)

ssc_f2 <- ssc_f[rord,cord]

saveRDS(colnames(ssc_f2),sprintf("%s/Heatmap_taxa_order_prokaryote.rds",save.dir))

cust_col <- foreach(i = 1:nrow(ssc_f2), .combine = rbind) %do% {
  ifelse(ssc_f2[i, ]==1,sscolor[rownames(ssc_f2)[i]],"white")
} 
rownames(cust_col) <- rownames(ssc_f2)


col_fun = colorRamp2(c(0, 1), c("white", "white"))

plss <- Heatmap(wid_f[rownames(ssc_f),plant_order], row_names_side = "left", column_names_side = "bottom", 
                col = col_fun,
                rect_gp = gpar(col="black"),
                cell_fun = function(j, i,x,y, width, height, fill) {
                  if(wid_f[rownames(ssc_f),plant_order][i,j] == 1){
                    grid.circle(x = x, y = y, r = 1.5*min(unit.c(width, height)),
                                gp = gpar(fill = "darkred"))
                  }
                },
                show_heatmap_legend = FALSE,
                cluster_columns = FALSE,
                cluster_rows = FALSE,
                show_row_dend = FALSE,
                show_column_dend = FALSE,
                #row_names_gp = gpar(cex=0.9, fontface = "bold"), 
                column_names_gp = gpar(cex=0.9, fontface = "italic"), 
                row_dend_width = unit(2, "cm"),
                column_dend_height = unit(2, "cm"), 
                column_title = "Host plant",
                column_names_rot = 45, #row_title = "Basins of attraction", 
                )







sscomp <- Heatmap(ssc_f, row_names_side = "left", column_names_side = "bottom", 
                  col = col_fun,
                  show_column_dend = FALSE,
                  
                  show_heatmap_legend = FALSE,
                  row_dend_side = "left", 
                  rect_gp = gpar(fill = as.vector(cust_col),
                                 col = "black"), 
                  row_names_gp = gpar(cex=0.9), 
                  column_names_gp = gpar(cex=0.9), 
                  row_dend_width = unit(2, "cm"),
                  column_dend_height = unit(2, "cm"), 
                  column_title = "Prokaryotes ( Family )",
                  column_names_rot = 45, row_title = "Stable states", 
                  clustering_distance_rows = jac_dist,
                  clustering_distance_columns = jac_dist, 
                  clustering_method_rows = "ward.D2",
                  clustering_method_columns = "ward.D2")

mh <- draw(sscomp+plss)

#check row order
all(rord==row_order(mh))

pdf(sprintf("%s/SSheatmap_Prok.pdf",save.dir),width = 12,h=7.3)
mh

dev.off()

