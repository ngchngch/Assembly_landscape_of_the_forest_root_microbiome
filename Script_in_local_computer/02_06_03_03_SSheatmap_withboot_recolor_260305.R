set.seed(1234)
library(ComplexHeatmap)
library(cluster)
library(circlize)
library(colorspace)


dir_03_10_02 <- "Output/03_10_graphics_states_flow_flow_Spl_250508"
dir_02_06 <- "Output_supercomputer/02_06_ELA"
dir_02_06_03_02 <- "Output_supercomputer/02_06_03_02_summarize_bootELA_bcom_DG"
####
save.dir <- "Output/02_06_03_03_SSheatmap_withboot_recolor_260305"
dir.create(save.dir)


##
#color scale

# 共通の break
brk <- seq(0, 1, 0.25)

# 明度・彩度を完全に共通化
L_vals <- c(100, 85, 70, 55, 40)   # 明度（高→低）
C_vals <- c(0,   30, 45, 55, 60)   # 彩度（白→鮮やか）

# 緑（Hue ≈ 130）
col_green <- hcl(
  h = 130,
  c = C_vals,
  l = L_vals
)

# オレンジ（Hue ≈ 40）
col_orange <- hcl(
  h = 40,
  c = C_vals*1.5,
  l = L_vals
)

col_fun_green <- colorRamp2(brk, col_green)
col_fun_orange <- colorRamp2(brk, col_orange)

###
fb <- "Fungi"

sscolor <- readRDS(sprintf("%s/states_colvector.rds",dir_03_10_02))

jac_dist <- function(x) {
  vegdist(x, method = "jaccard")
}

ssc_f <- readRDS(sprintf("%s/%s/SScomposition_%s.rds",dir_02_06,fb,fb))
boot_ss_freq <- readRDS(sprintf("%s/basincomposition_freq_%s_summary.rds",dir_02_06_03_02,fb))

ht_ord <- draw(Heatmap(ssc_f,clustering_distance_rows = jac_dist,
                       clustering_distance_columns = jac_dist, 
                       clustering_method_rows = "ward.D2",
                       clustering_method_columns = "ward.D2"))

cord <- colnames(ssc_f)[column_order(ht_ord)]
rord <- rownames(ssc_f)[row_order(ht_ord)]

ls_hmap <- list(NULL)
plnam <- c()
#####
pl <- 1
  
  if(all(colnames(boot_ss_freq) %in% cord)){
    ssc_f2 <- boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],cord]  
  }else{
    ssc_f2 <- cbind(boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],cord],
                    boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],
                                 setdiff(which(!colnames(boot_ss_freq) %in% cord),c(1:3))])
  }
  
  rownames(ssc_f2) <- boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],"new_ssid"]
  
  saveRDS(ssc_f2[na.omit(match(rord,rownames(ssc_f2))),],
          sprintf("%s/SSOTUcomp_heatmap_full_%s_%s.rds",save.dir,unique(boot_ss_freq$plant)[pl],fb)
  )
  
  #make full ss comp
  fullss <- t(sapply(boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],"ssid"],id2bin,34))
  dimnames(fullss) <- list(boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],"new_ssid"],
                           colnames(boot_ss_freq)[-c(1:3)])
  
  if(ncol(fullss) == length(cord)){
    binary_mat <-   fullss[na.omit(match(rord,rownames(fullss))),cord]
  }else{
    binary_mat <-    cbind(matrix(fullss[na.omit(match(rord,rownames(fullss))),cord],
                                  nrow=length(na.omit(match(rord,rownames(fullss))))),
                           matrix( fullss[na.omit(match(rord,rownames(fullss))),which(!colnames(fullss) %in% cord)],
                                   nrow=length(na.omit(match(rord,rownames(fullss)))))
    )
  }
  
  
  plnam <- c(plnam,unique(boot_ss_freq$plant)[pl])
  ls_hmap[[pl]] <- list(bootsterap=ssc_f2[na.omit(match(rord,rownames(ssc_f2))),],
                        full_membership=binary_mat)
  
  sscomp <- Heatmap(
    ssc_f2[na.omit(match(rord,rownames(ssc_f2))),],
    row_names_side = "left",
    column_names_side = "bottom",
    col = col_fun_orange,
    show_column_dend = FALSE,
    show_row_dend = FALSE,
    show_heatmap_legend = TRUE,
    
    rect_gp = gpar(col = "black"),
    row_names_gp = gpar(cex = 0.9, fontface = "bold"),
    column_names_gp = gpar(cex = 0.9, fontface = "bold"),
    row_dend_width = unit(2, "cm"),
    column_dend_height = unit(2, "cm"),
    name = "Bootstrap\nsupport",
    column_title = sprintf("%s (Family) in %s", fb, unique(boot_ss_freq$plant)[pl]),
    column_names_rot = 45,
    row_title = "Basin bottoms",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    
    cell_fun = function(j, i, x, y, width, height, fill) {
      if(binary_mat[i, j] == 1){
        grid.text(
          "*",
          x = x,
          y = y,
          gp = gpar(fontsize = 14, col = "black", fontface = "bold")
        )
      }
    }
  )
  pdf(sprintf("%s/SSOTUcomp_heatmap_full_%s_%s.pdf",save.dir,unique(boot_ss_freq$plant)[pl],fb),
      width = 11, height = 3+nrow(ssc_f2)/4)
  sscomp
  dev.off()
#########
  pl <- 2
  
  if(all(colnames(boot_ss_freq) %in% cord)){
    ssc_f2 <- boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],cord]  
  }else{
    ssc_f2 <- cbind(boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],cord],
                    boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],
                                 setdiff(which(!colnames(boot_ss_freq) %in% cord),c(1:3))])
  }
  
  rownames(ssc_f2) <- boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],"new_ssid"]
  
  #make full ss comp
  fullss <- t(sapply(boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],"ssid"],id2bin,34))
  dimnames(fullss) <- list(boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],"new_ssid"],
                           colnames(boot_ss_freq)[-c(1:3)])
  
  if(ncol(fullss) == length(cord)){
    binary_mat <-   fullss[na.omit(match(rord,rownames(fullss))),cord]
  }else{
    binary_mat <-    cbind(matrix(fullss[na.omit(match(rord,rownames(fullss))),cord],
                                  nrow=length(na.omit(match(rord,rownames(fullss))))),
                           matrix( fullss[na.omit(match(rord,rownames(fullss))),which(!colnames(fullss) %in% cord)],
                                   nrow=length(na.omit(match(rord,rownames(fullss)))))
    )
  }
  
  
  plnam <- c(plnam,unique(boot_ss_freq$plant)[pl])
  ls_hmap[[pl]] <- list(bootsterap=ssc_f2[na.omit(match(rord,rownames(ssc_f2))),],
                        full_membership=binary_mat)
  
  sscomp <- Heatmap(
    ssc_f2[na.omit(match(rord,rownames(ssc_f2))),],
    row_names_side = "left",
    column_names_side = "bottom",
    col = col_fun_orange,
    show_column_dend = FALSE,
    show_row_dend = FALSE,
    show_heatmap_legend = TRUE,
    
    rect_gp = gpar(col = "black"),
    row_names_gp = gpar(cex = 0.9, fontface = "bold"),
    column_names_gp = gpar(cex = 0.9, fontface = "bold"),
    row_dend_width = unit(2, "cm"),
    column_dend_height = unit(2, "cm"),
    name = "Bootstrap\nsupport",
    column_title = sprintf("%s (Family) in %s", fb, unique(boot_ss_freq$plant)[pl]),
    column_names_rot = 45,
    row_title = "Basin bottoms",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    
    cell_fun = function(j, i, x, y, width, height, fill) {
      if(binary_mat[i, j] == 1){
        grid.text(
          "*",
          x = x,
          y = y,
          gp = gpar(fontsize = 14, col = "black", fontface = "bold")
        )
      }
    }
  )
  pdf(sprintf("%s/SSOTUcomp_heatmap_full_%s_%s.pdf",save.dir,unique(boot_ss_freq$plant)[pl],fb),
      width = 11, height = 3+nrow(ssc_f2)/4)
  sscomp
  dev.off()
  #####
  pl <- 3
  
  if(all(colnames(boot_ss_freq) %in% cord)){
    ssc_f2 <- boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],cord]  
  }else{
    ssc_f2 <- cbind(boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],cord],
                    boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],
                                 setdiff(which(!colnames(boot_ss_freq) %in% cord),c(1:3))])
  }
  
  rownames(ssc_f2) <- boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],"new_ssid"]
  
  
  #make full ss comp
  fullss <- t(sapply(boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],"ssid"],id2bin,34))
  dimnames(fullss) <- list(boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],"new_ssid"],
                           colnames(boot_ss_freq)[-c(1:3)])
  
  if(ncol(fullss) == length(cord)){
    binary_mat <-   fullss[na.omit(match(rord,rownames(fullss))),cord]
  }else{
    binary_mat <-    cbind(matrix(fullss[na.omit(match(rord,rownames(fullss))),cord],
                                  nrow=length(na.omit(match(rord,rownames(fullss))))),
                           matrix( fullss[na.omit(match(rord,rownames(fullss))),which(!colnames(fullss) %in% cord)],
                                   nrow=length(na.omit(match(rord,rownames(fullss)))))
    )}
  
  
  plnam <- c(plnam,unique(boot_ss_freq$plant)[pl])
  ls_hmap[[pl]] <- list(bootsterap=ssc_f2[na.omit(match(rord,rownames(ssc_f2))),],
                        full_membership=binary_mat)
  
  sscomp <- Heatmap(
    ssc_f2[na.omit(match(rord,rownames(ssc_f2))),],
    row_names_side = "left",
    column_names_side = "bottom",
    col = col_fun_orange,
    show_column_dend = FALSE,
    show_row_dend = FALSE,
    show_heatmap_legend = TRUE,
    
    rect_gp = gpar(col = "black"),
    row_names_gp = gpar(cex = 0.9, fontface = "bold"),
    column_names_gp = gpar(cex = 0.9, fontface = "bold"),
    row_dend_width = unit(2, "cm"),
    column_dend_height = unit(2, "cm"),
    name = "Bootstrap\nsupport",
    column_title = sprintf("%s (Family) in %s", fb, unique(boot_ss_freq$plant)[pl]),
    column_names_rot = 45,
    row_title = "Basin bottoms",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    
    cell_fun = function(j, i, x, y, width, height, fill) {
      if(binary_mat[i, j] == 1){
        grid.text(
          "*",
          x = x,
          y = y,
          gp = gpar(fontsize = 14, col = "black", fontface = "bold")
        )
      }
    }
  )
  pdf(sprintf("%s/SSOTUcomp_heatmap_full_%s_%s.pdf",save.dir,unique(boot_ss_freq$plant)[pl],fb),
      width = 11, height = 3+nrow(ssc_f2)/4)
  sscomp
  dev.off()
  #####
  pl <- 4
  
  if(all(colnames(boot_ss_freq) %in% cord)){
    ssc_f2 <- boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],cord]  
  }else{
    ssc_f2 <- cbind(boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],cord],
                    boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],
                                 setdiff(which(!colnames(boot_ss_freq) %in% cord),c(1:3))])
  }
  
  rownames(ssc_f2) <- boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],"new_ssid"]
  
  #make full ss comp
  fullss <- t(sapply(boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],"ssid"],id2bin,34))
  dimnames(fullss) <- list(boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],"new_ssid"],
                           colnames(boot_ss_freq)[-c(1:3)])
  
  if(ncol(fullss) == length(cord)){
    binary_mat <-   fullss[na.omit(match(rord,rownames(fullss))),cord]
  }else{
    binary_mat <-   cbind(matrix(fullss[na.omit(match(rord,rownames(fullss))),cord],
                                 nrow=length(na.omit(match(rord,rownames(fullss))))),
                          matrix( fullss[na.omit(match(rord,rownames(fullss))),which(!colnames(fullss) %in% cord)],
                                  nrow=length(na.omit(match(rord,rownames(fullss)))))
    )
  }
  
  
  plnam <- c(plnam,unique(boot_ss_freq$plant)[pl])
  ls_hmap[[pl]] <- list(bootsterap=ssc_f2[na.omit(match(rord,rownames(ssc_f2))),],
                        full_membership=binary_mat)
  
  sscomp <- Heatmap(
    ssc_f2[na.omit(match(rord,rownames(ssc_f2))),],
    row_names_side = "left",
    column_names_side = "bottom",
    col = col_fun_orange,
    show_column_dend = FALSE,
    show_row_dend = FALSE,
    show_heatmap_legend = TRUE,
    
    rect_gp = gpar(col = "black"),
    row_names_gp = gpar(cex = 0.9, fontface = "bold"),
    column_names_gp = gpar(cex = 0.9, fontface = "bold"),
    row_dend_width = unit(2, "cm"),
    column_dend_height = unit(2, "cm"),
    name = "Bootstrap\nsupport",
    column_title = sprintf("%s (Family) in %s", fb, unique(boot_ss_freq$plant)[pl]),
    column_names_rot = 45,
    row_title = "Basin bottoms",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    
    cell_fun = function(j, i, x, y, width, height, fill) {
      if(binary_mat[i, j] == 1){
        grid.text(
          "*",
          x = x,
          y = y,
          gp = gpar(fontsize = 14, col = "black", fontface = "bold")
        )
      }
    }
  )
  pdf(sprintf("%s/SSOTUcomp_heatmap_full_%s_%s.pdf",save.dir,unique(boot_ss_freq$plant)[pl],fb),
      width = 11, height = 3+nrow(ssc_f2)/4)
  sscomp
  dev.off()
  #####
  pl <- 5
  
  if(all(colnames(boot_ss_freq) %in% cord)){
    ssc_f2 <- boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],cord]  
  }else{
    ssc_f2 <- cbind(boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],cord],
                    boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],
                                 setdiff(which(!colnames(boot_ss_freq) %in% cord),c(1:3))])
  }
  
  rownames(ssc_f2) <- boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],"new_ssid"]
  
  #make full ss comp
  fullss <- t(sapply(boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],"ssid"],id2bin,34))
  dimnames(fullss) <- list(boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],"new_ssid"],
                           colnames(boot_ss_freq)[-c(1:3)])
  
  if(ncol(fullss) == length(cord)){
    binary_mat <-   fullss[na.omit(match(rord,rownames(fullss))),cord]
  }else{
    binary_mat <-   cbind(matrix(fullss[na.omit(match(rord,rownames(fullss))),cord],
                                 nrow=length(na.omit(match(rord,rownames(fullss))))),
                          matrix( fullss[na.omit(match(rord,rownames(fullss))),which(!colnames(fullss) %in% cord)],
                                 nrow=length(na.omit(match(rord,rownames(fullss)))))
                         )
  }
  
  
  plnam <- c(plnam,unique(boot_ss_freq$plant)[pl])
  ls_hmap[[pl]] <- list(bootsterap=ssc_f2[na.omit(match(rord,rownames(ssc_f2))),],
                        full_membership=binary_mat)
  
  sscomp <- Heatmap(
    ssc_f2[na.omit(match(rord,rownames(ssc_f2))),],
    row_names_side = "left",
    column_names_side = "bottom",
    col = col_fun_orange,
    show_column_dend = FALSE,
    show_row_dend = FALSE,
    show_heatmap_legend = TRUE,
    
    rect_gp = gpar(col = "black"),
    row_names_gp = gpar(cex = 0.9, fontface = "bold"),
    column_names_gp = gpar(cex = 0.9, fontface = "bold"),
    row_dend_width = unit(2, "cm"),
    column_dend_height = unit(2, "cm"),
    name = "Bootstrap\nsupport",
    column_title = sprintf("%s (Family) in %s", fb, unique(boot_ss_freq$plant)[pl]),
    column_names_rot = 45,
    row_title = "Basin bottoms",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    
    cell_fun = function(j, i, x, y, width, height, fill) {
      if(binary_mat[i, j] == 1){
        grid.text(
          "*",
          x = x,
          y = y,
          gp = gpar(fontsize = 14, col = "black", fontface = "bold")
        )
      }
    }
  )
  pdf(sprintf("%s/SSOTUcomp_heatmap_full_%s_%s.pdf",save.dir,unique(boot_ss_freq$plant)[pl],fb),
      width = 11, height = 3+nrow(ssc_f2)/4)
  sscomp
  dev.off()
  #####
  pl <- 6
  
  if(all(colnames(boot_ss_freq) %in% cord)){
    ssc_f2 <- boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],cord]  
  }else{
    ssc_f2 <- cbind(boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],cord],
                    boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],
                                 setdiff(which(!colnames(boot_ss_freq) %in% cord),c(1:3))])
  }
  
  rownames(ssc_f2) <- boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],"new_ssid"]
  
  #make full ss comp
  fullss <- t(sapply(boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],"ssid"],id2bin,34))
  dimnames(fullss) <- list(boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],"new_ssid"],
                           colnames(boot_ss_freq)[-c(1:3)])
  
  if(ncol(fullss) == length(cord)){
    binary_mat <-   fullss[na.omit(match(rord,rownames(fullss))),cord]
  }else{
    binary_mat <-   cbind(matrix(fullss[na.omit(match(rord,rownames(fullss))),cord],
                                 nrow=length(na.omit(match(rord,rownames(fullss))))),
                          matrix( fullss[na.omit(match(rord,rownames(fullss))),which(!colnames(fullss) %in% cord)],
                                  nrow=length(na.omit(match(rord,rownames(fullss)))))
    )
  }
  
  
  plnam <- c(plnam,unique(boot_ss_freq$plant)[pl])
  ls_hmap[[pl]] <- list(bootsterap=ssc_f2[na.omit(match(rord,rownames(ssc_f2))),],
                                   full_membership=binary_mat)
  
  sscomp <- Heatmap(
    ssc_f2[na.omit(match(rord,rownames(ssc_f2))),],
    row_names_side = "left",
    column_names_side = "bottom",
    col = col_fun_orange,
    show_column_dend = FALSE,
    show_row_dend = FALSE,
    show_heatmap_legend = TRUE,
    
    rect_gp = gpar(col = "black"),
    row_names_gp = gpar(cex = 0.9, fontface = "bold"),
    column_names_gp = gpar(cex = 0.9, fontface = "bold"),
    row_dend_width = unit(2, "cm"),
    column_dend_height = unit(2, "cm"),
    name = "Bootstrap\nsupport",
    column_title = sprintf("%s (Family) in %s", fb, unique(boot_ss_freq$plant)[pl]),
    column_names_rot = 45,
    row_title = "Basin bottoms",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    
    cell_fun = function(j, i, x, y, width, height, fill) {
      if(binary_mat[i, j] == 1){
        grid.text(
          "*",
          x = x,
          y = y,
          gp = gpar(fontsize = 14, col = "black", fontface = "bold")
        )
      }
    }
  )
  pdf(sprintf("%s/SSOTUcomp_heatmap_full_%s_%s.pdf",save.dir,unique(boot_ss_freq$plant)[pl],fb),
      width = 11, height = 3+nrow(ssc_f2)/4)
  sscomp
  dev.off()
##
names(ls_hmap) <- plnam
saveRDS(ls_hmap, sprintf("%s/SSOTUcomp_heatmap_full_list.rds",save.dir))
  
###
fb <- "Prokaryota"

sscolor <- readRDS(sprintf("%s/states_colvector.rds",dir_03_10_02))

jac_dist <- function(x) {
  vegdist(x, method = "jaccard")
}

ssc_f <- readRDS(sprintf("%s/%s/SScomposition_%s.rds",dir_02_06,fb,fb))
boot_ss_freq <- readRDS(sprintf("%s/basincomposition_freq_%s_summary.rds",dir_02_06_03_02,fb))

ht_ord <- draw(Heatmap(ssc_f,clustering_distance_rows = jac_dist,
                       clustering_distance_columns = jac_dist, 
                       clustering_method_rows = "ward.D2",
                       clustering_method_columns = "ward.D2"))



sscolor <- readRDS(sprintf("%s/states_colvector.rds",dir_03_10_02))

jac_dist <- function(x) {
  vegdist(x, method = "jaccard")
}

ssc_f <- readRDS(sprintf("%s/%s/SScomposition_%s.rds",dir_02_06,fb,fb))
boot_ss_freq <- readRDS(sprintf("%s/basincomposition_freq_%s_summary.rds",dir_02_06_03_02,fb))

ht_ord <- draw(Heatmap(ssc_f,clustering_distance_rows = jac_dist,
                       clustering_distance_columns = jac_dist, 
                       clustering_method_rows = "ward.D2",
                       clustering_method_columns = "ward.D2"))

cord <- colnames(ssc_f)[column_order(ht_ord)]
rord <- rownames(ssc_f)[row_order(ht_ord)]

ls_hmap <- list(NULL)
plnam <- c()

#####
pl <- 1

if(all(colnames(boot_ss_freq) %in% cord)){
  ssc_f2 <- boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],cord]  
}else{
  ssc_f2 <- cbind(boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],cord],
                  boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],
                               setdiff(which(!colnames(boot_ss_freq) %in% cord),c(1:3))])
}

rownames(ssc_f2) <- boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],"new_ssid"]
#make full ss comp
fullss <- t(sapply(boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],"ssid"],id2bin,45))
dimnames(fullss) <- list(boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],"new_ssid"],
                         colnames(boot_ss_freq)[-c(1:3)])

if(ncol(fullss) == length(cord)){
  binary_mat <-   fullss[na.omit(match(rord,rownames(fullss))),cord]
}else{
  binary_mat <-   cbind(matrix(fullss[na.omit(match(rord,rownames(fullss))),cord],
                               nrow=length(na.omit(match(rord,rownames(fullss))))),
                        matrix( fullss[na.omit(match(rord,rownames(fullss))),which(!colnames(fullss) %in% cord)],
                                nrow=length(na.omit(match(rord,rownames(fullss)))))
  )
}


plnam <- c(plnam,unique(boot_ss_freq$plant)[pl])
ls_hmap[[pl]] <- list(bootsterap=ssc_f2[na.omit(match(rord,rownames(ssc_f2))),],
                      full_membership=binary_mat)

sscomp <- Heatmap(
  ssc_f2[na.omit(match(rord,rownames(ssc_f2))),],
  row_names_side = "left",
  column_names_side = "bottom",
  col = col_fun_green,
  show_column_dend = FALSE,
  show_row_dend = FALSE,
  show_heatmap_legend = TRUE,
  
  rect_gp = gpar(col = "black"),
  row_names_gp = gpar(cex = 0.9, fontface = "bold"),
  column_names_gp = gpar(cex = 0.9, fontface = "bold"),
  row_dend_width = unit(2, "cm"),
  column_dend_height = unit(2, "cm"),
  name = "Bootstrap\nsupport",
  column_title = sprintf("%s (Family) in %s", fb, unique(boot_ss_freq$plant)[pl]),
  column_names_rot = 45,
  row_title = "Basin bottoms",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    if(binary_mat[i, j] == 1){
      grid.text(
        "*",
        x = x,
        y = y,
        gp = gpar(fontsize = 14, col = "black", fontface = "bold")
      )
    }
  }
)
pdf(sprintf("%s/SSOTUcomp_heatmap_full_%s_%s.pdf",save.dir,unique(boot_ss_freq$plant)[pl],fb),
    width = 11, height = 3+nrow(ssc_f2)/4)
sscomp
dev.off()
#########
pl <- 2

if(all(colnames(boot_ss_freq) %in% cord)){
  ssc_f2 <- boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],cord]  
}else{
  ssc_f2 <- cbind(boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],cord],
                  boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],
                               setdiff(which(!colnames(boot_ss_freq) %in% cord),c(1:3))])
}

rownames(ssc_f2) <- boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],"new_ssid"]
#make full ss comp
fullss <- t(sapply(boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],"ssid"],id2bin,45))
dimnames(fullss) <- list(boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],"new_ssid"],
                         colnames(boot_ss_freq)[-c(1:3)])

if(ncol(fullss) == length(cord)){
  binary_mat <-   fullss[na.omit(match(rord,rownames(fullss))),cord]
}else{
  binary_mat <-   cbind(matrix(fullss[na.omit(match(rord,rownames(fullss))),cord],
                               nrow=length(na.omit(match(rord,rownames(fullss))))),
                        matrix( fullss[na.omit(match(rord,rownames(fullss))),which(!colnames(fullss) %in% cord)],
                                nrow=length(na.omit(match(rord,rownames(fullss)))))
  )
}


plnam <- c(plnam,unique(boot_ss_freq$plant)[pl])
ls_hmap[[pl]] <- list(bootsterap=ssc_f2[na.omit(match(rord,rownames(ssc_f2))),],
                      full_membership=binary_mat)

sscomp <- Heatmap(
  ssc_f2[na.omit(match(rord,rownames(ssc_f2))),],
  row_names_side = "left",
  column_names_side = "bottom",
  col = col_fun_green,
  show_column_dend = FALSE,
  show_row_dend = FALSE,
  show_heatmap_legend = TRUE,
  
  rect_gp = gpar(col = "black"),
  row_names_gp = gpar(cex = 0.9, fontface = "bold"),
  column_names_gp = gpar(cex = 0.9, fontface = "bold"),
  row_dend_width = unit(2, "cm"),
  column_dend_height = unit(2, "cm"),
  name = "Bootstrap\nsupport",
  column_title = sprintf("%s (Family) in %s", fb, unique(boot_ss_freq$plant)[pl]),
  column_names_rot = 45,
  row_title = "Basin bottoms",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    if(binary_mat[i, j] == 1){
      grid.text(
        "*",
        x = x,
        y = y,
        gp = gpar(fontsize = 14, col = "black", fontface = "bold")
      )
    }
  }
)
pdf(sprintf("%s/SSOTUcomp_heatmap_full_%s_%s.pdf",save.dir,unique(boot_ss_freq$plant)[pl],fb),
    width = 11, height = 3+nrow(ssc_f2)/4)
sscomp
dev.off()
#####
pl <- 3

if(all(colnames(boot_ss_freq) %in% cord)){
  ssc_f2 <- boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],cord]  
}else{
  ssc_f2 <- cbind(boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],cord],
                  boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],
                               setdiff(which(!colnames(boot_ss_freq) %in% cord),c(1:3))])
}

rownames(ssc_f2) <- boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],"new_ssid"]
#make full ss comp
fullss <- t(sapply(boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],"ssid"],id2bin,45))
dimnames(fullss) <- list(boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],"new_ssid"],
                         colnames(boot_ss_freq)[-c(1:3)])

if(ncol(fullss) == length(cord)){
  binary_mat <-   fullss[na.omit(match(rord,rownames(fullss))),cord]
}else{
  binary_mat <-   cbind(matrix(fullss[na.omit(match(rord,rownames(fullss))),cord],
                               nrow=length(na.omit(match(rord,rownames(fullss))))),
                        matrix( fullss[na.omit(match(rord,rownames(fullss))),which(!colnames(fullss) %in% cord)],
                                nrow=length(na.omit(match(rord,rownames(fullss)))))
  )
}


plnam <- c(plnam,unique(boot_ss_freq$plant)[pl])
ls_hmap[[pl]] <- list(bootsterap=ssc_f2[na.omit(match(rord,rownames(ssc_f2))),],
                      full_membership=binary_mat)

sscomp <- Heatmap(
  ssc_f2[na.omit(match(rord,rownames(ssc_f2))),],
  row_names_side = "left",
  column_names_side = "bottom",
  col = col_fun_green,
  show_column_dend = FALSE,
  show_row_dend = FALSE,
  show_heatmap_legend = TRUE,
  
  rect_gp = gpar(col = "black"),
  row_names_gp = gpar(cex = 0.9, fontface = "bold"),
  column_names_gp = gpar(cex = 0.9, fontface = "bold"),
  row_dend_width = unit(2, "cm"),
  column_dend_height = unit(2, "cm"),
  name = "Bootstrap\nsupport",
  column_title = sprintf("%s (Family) in %s", fb, unique(boot_ss_freq$plant)[pl]),
  column_names_rot = 45,
  row_title = "Basin bottoms",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    if(binary_mat[i, j] == 1){
      grid.text(
        "*",
        x = x,
        y = y,
        gp = gpar(fontsize = 14, col = "black", fontface = "bold")
      )
    }
  }
)
pdf(sprintf("%s/SSOTUcomp_heatmap_full_%s_%s.pdf",save.dir,unique(boot_ss_freq$plant)[pl],fb),
    width = 11, height = 3+nrow(ssc_f2)/4)
sscomp
dev.off()
#####
pl <- 4

if(all(colnames(boot_ss_freq) %in% cord)){
  ssc_f2 <- boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],cord]  
}else{
  ssc_f2 <- cbind(boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],cord],
                  boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],
                               setdiff(which(!colnames(boot_ss_freq) %in% cord),c(1:3))])
}

rownames(ssc_f2) <- boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],"new_ssid"]
#make full ss comp
fullss <- t(sapply(boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],"ssid"],id2bin,45))
dimnames(fullss) <- list(boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],"new_ssid"],
                         colnames(boot_ss_freq)[-c(1:3)])

if(ncol(fullss) == length(cord)){
  binary_mat <-   fullss[na.omit(match(rord,rownames(fullss))),cord]
}else{
  binary_mat <-   cbind(matrix(fullss[na.omit(match(rord,rownames(fullss))),cord],
                               nrow=length(na.omit(match(rord,rownames(fullss))))),
                        matrix( fullss[na.omit(match(rord,rownames(fullss))),which(!colnames(fullss) %in% cord)],
                                nrow=length(na.omit(match(rord,rownames(fullss)))))
  )
}


plnam <- c(plnam,unique(boot_ss_freq$plant)[pl])
ls_hmap[[pl]] <- list(bootsterap=ssc_f2[na.omit(match(rord,rownames(ssc_f2))),],
                      full_membership=binary_mat)

sscomp <- Heatmap(
  ssc_f2[na.omit(match(rord,rownames(ssc_f2))),],
  row_names_side = "left",
  column_names_side = "bottom",
  col = col_fun_green,
  show_column_dend = FALSE,
  show_row_dend = FALSE,
  show_heatmap_legend = TRUE,
  
  rect_gp = gpar(col = "black"),
  row_names_gp = gpar(cex = 0.9, fontface = "bold"),
  column_names_gp = gpar(cex = 0.9, fontface = "bold"),
  row_dend_width = unit(2, "cm"),
  column_dend_height = unit(2, "cm"),
  name = "Bootstrap\nsupport",
  column_title = sprintf("%s (Family) in %s", fb, unique(boot_ss_freq$plant)[pl]),
  column_names_rot = 45,
  row_title = "Basin bottoms",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    if(binary_mat[i, j] == 1){
      grid.text(
        "*",
        x = x,
        y = y,
        gp = gpar(fontsize = 14, col = "black", fontface = "bold")
      )
    }
  }
)
pdf(sprintf("%s/SSOTUcomp_heatmap_full_%s_%s.pdf",save.dir,unique(boot_ss_freq$plant)[pl],fb),
    width = 11, height = 3+nrow(ssc_f2)/4)
sscomp
dev.off()
#####
pl <- 5

if(all(colnames(boot_ss_freq) %in% cord)){
  ssc_f2 <- boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],cord]  
}else{
  ssc_f2 <- cbind(boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],cord],
                  boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],
                               setdiff(which(!colnames(boot_ss_freq) %in% cord),c(1:3))])
}

rownames(ssc_f2) <- boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],"new_ssid"]
#make full ss comp
fullss <- t(sapply(boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],"ssid"],id2bin,45))
dimnames(fullss) <- list(boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],"new_ssid"],
                         colnames(boot_ss_freq)[-c(1:3)])

if(ncol(fullss) == length(cord)){
  binary_mat <-   fullss[na.omit(match(rord,rownames(fullss))),cord]
}else{
  binary_mat <-   cbind(matrix(fullss[na.omit(match(rord,rownames(fullss))),cord],
                               nrow=length(na.omit(match(rord,rownames(fullss))))),
                        matrix( fullss[na.omit(match(rord,rownames(fullss))),which(!colnames(fullss) %in% cord)],
                                nrow=length(na.omit(match(rord,rownames(fullss)))))
  )
}


plnam <- c(plnam,unique(boot_ss_freq$plant)[pl])
ls_hmap[[pl]] <- list(bootsterap=ssc_f2[na.omit(match(rord,rownames(ssc_f2))),],
                      full_membership=binary_mat)

sscomp <- Heatmap(
  ssc_f2[na.omit(match(rord,rownames(ssc_f2))),],
  row_names_side = "left",
  column_names_side = "bottom",
  col = col_fun_green,
  show_column_dend = FALSE,
  show_row_dend = FALSE,
  show_heatmap_legend = TRUE,
  
  rect_gp = gpar(col = "black"),
  row_names_gp = gpar(cex = 0.9, fontface = "bold"),
  column_names_gp = gpar(cex = 0.9, fontface = "bold"),
  row_dend_width = unit(2, "cm"),
  column_dend_height = unit(2, "cm"),
  name = "Bootstrap\nsupport",
  column_title = sprintf("%s (Family) in %s", fb, unique(boot_ss_freq$plant)[pl]),
  column_names_rot = 45,
  row_title = "Basin bottoms",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    if(binary_mat[i, j] == 1){
      grid.text(
        "*",
        x = x,
        y = y,
        gp = gpar(fontsize = 14, col = "black", fontface = "bold")
      )
    }
  }
)
pdf(sprintf("%s/SSOTUcomp_heatmap_full_%s_%s.pdf",save.dir,unique(boot_ss_freq$plant)[pl],fb),
    width = 11, height = 3+nrow(ssc_f2)/4)
sscomp
dev.off()
#####
pl <- 6

if(all(colnames(boot_ss_freq) %in% cord)){
  ssc_f2 <- boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],cord]  
}else{
  ssc_f2 <- cbind(boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],cord],
                  boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],
                               setdiff(which(!colnames(boot_ss_freq) %in% cord),c(1:3))])
}

rownames(ssc_f2) <- boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],"new_ssid"]
#make full ss comp
fullss <- t(sapply(boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],"ssid"],id2bin,45))
dimnames(fullss) <- list(boot_ss_freq[boot_ss_freq$plant==unique(boot_ss_freq$plant)[pl],"new_ssid"],
                         colnames(boot_ss_freq)[-c(1:3)])

if(ncol(fullss) == length(cord)){
  binary_mat <-   fullss[na.omit(match(rord,rownames(fullss))),cord]
}else{
  binary_mat <-   cbind(matrix(fullss[na.omit(match(rord,rownames(fullss))),cord],
                               nrow=length(na.omit(match(rord,rownames(fullss))))),
                        matrix( fullss[na.omit(match(rord,rownames(fullss))),which(!colnames(fullss) %in% cord)],
                                nrow=length(na.omit(match(rord,rownames(fullss)))))
  )
}


plnam <- c(plnam,unique(boot_ss_freq$plant)[pl])
ls_hmap[[pl]] <- list(bootsterap=ssc_f2[na.omit(match(rord,rownames(ssc_f2))),],
                      full_membership=binary_mat)

sscomp <- Heatmap(
  ssc_f2[na.omit(match(rord,rownames(ssc_f2))),],
  row_names_side = "left",
  column_names_side = "bottom",
  col = col_fun_green,
  show_column_dend = FALSE,
  show_row_dend = FALSE,
  show_heatmap_legend = TRUE,
  
  rect_gp = gpar(col = "black"),
  row_names_gp = gpar(cex = 0.9, fontface = "bold"),
  column_names_gp = gpar(cex = 0.9, fontface = "bold"),
  row_dend_width = unit(2, "cm"),
  column_dend_height = unit(2, "cm"),
  name = "Bootstrap\nsupport",
  column_title = sprintf("%s (Family) in %s", fb, unique(boot_ss_freq$plant)[pl]),
  column_names_rot = 45,
  row_title = "Basin bottoms",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    if(binary_mat[i, j] == 1){
      grid.text(
        "*",
        x = x,
        y = y,
        gp = gpar(fontsize = 14, col = "black", fontface = "bold")
      )
    }
  }
)
pdf(sprintf("%s/SSOTUcomp_heatmap_full_%s_%s.pdf",save.dir,unique(boot_ss_freq$plant)[pl],fb),
    width = 11, height = 3+nrow(ssc_f2)/4)
sscomp
dev.off()

##
names(ls_hmap) <- plnam
saveRDS(ls_hmap, sprintf("%s/SSOTUcomp_heatmap_full_list.rds",save.dir))
