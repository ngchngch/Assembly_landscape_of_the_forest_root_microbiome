library(ComplexHeatmap)
library(cluster)

dir_02_06 <- "Output_supercomputer/02_06_ELA"

####
save.dir <- "Output/02_06_02_Order_SSheatmap_250918"
dir.create(save.dir)

##
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


ht_ord <- draw(Heatmap(ssc_f,clustering_distance_rows = jac_dist,
                       clustering_distance_columns = jac_dist, 
                       clustering_method_rows = "ward.D2",
                       clustering_method_columns = "ward.D2"))

rord <- row_order(ht_ord)
cord <- column_order(ht_ord)

ssc_f2 <- ssc_f[rord,cord]

saveRDS(colnames(ssc_f2),sprintf("%s/Heatmap_taxa_order_Fungi.rds",save.dir))

####
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