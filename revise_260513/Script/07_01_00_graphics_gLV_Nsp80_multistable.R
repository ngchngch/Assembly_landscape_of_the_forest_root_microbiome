# Set the working directory
setwd("/Users/ngch/Desktop/Now_using/data/sugaFB_ELA/revise_260513")

# Load required packages
library(tidyverse)
library(igraph)
library(ggplot2)
library(patchwork)
library(Rcpp)
library(RcppArmadillo)
library(vegan)
library(cluster)
library(parallel)

# Create the output directory
save.dir <- "Analysis/07_01_00_graphics_gLV_Nsp80_multistable"
dir.create(save.dir, showWarnings = FALSE, recursive = TRUE)

set.seed(123)

########################################################################
ldir <- list.dirs(dir_07_01);ldir <- ldir[grep("candidate",ldir)]
lf <- list.files(ldir,full.names = TRUE,recursive = TRUE); lf <- lf[grepl("dynamics_rank",lf)&grepl("rds",lf)]
nam <- sapply(strsplit(lf,"/"),"[",3)
for(i in 1:length(lf)){
  mat_obs <- readRDS(lf[i])
  
  # =============================================================================
  # 10. Visualization
  # =============================================================================
  Nsp <- 80
  obs_steps <- seq_len(nrow(mat_obs))+100
  
  ts_df <- as.data.frame(mat_obs) |>
    setNames(paste0("sp", seq_len(Nsp))) |>
    mutate(time = obs_steps) |>
    pivot_longer(
      -time,
      names_to = "species",
      values_to = "abundance"
    )
  
  # =============================================================================
  # 10-1. Visualize time series of all species
  # =============================================================================
  p1 <- ggplot(ts_df, aes(x = time, y = abundance, color = species)) +
    geom_line(alpha = 0.7, linewidth = 0.3) +
    
    scale_color_viridis_d(option = "plasma") +
    labs(
      x = "Time step",
      y = "Abundance",
      color = "Species"
    ) +
    theme_bw(base_size = 10) +
    theme(legend.position = "none",
          aspect.ratio =0.6,
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14))+
    scale_x_continuous(
      limits = c(101, 3000),
      breaks = c(100, 1000,2000,3000),
      expand = expansion(mult = 0.02)
    ) +
    coord_cartesian(xlim = c(min(ts_df$time), max(ts_df$time)))
  
  p1
  
  ggsave(
    filename = sprintf("%s/dynamics_all_%s.pdf", save.dir,nam[i]),
    plot = p1,
    height = 6,
    width = 8
  )
  
  # =============================================================================
  # 10-6. Calculate PCoA
  # =============================================================================
  bc_dist_matrix <- vegan::vegdist(mat_obs, method = "bray")
  
  pcoa_result <- cmdscale(bc_dist_matrix, k = 2, eig = TRUE)
  
  eigenvalues <- pcoa_result$eig
  eigenvalues_positive <- pmax(eigenvalues, 0)
  var_explained <- eigenvalues_positive / sum(eigenvalues_positive) * 100
  
  pcoa_df <- data.frame(
    PC1 = pcoa_result$points[, 1],
    PC2 = pcoa_result$points[, 2],
    time = seq_len(nrow(mat_obs))
  )
  
  cat(sprintf("PCoA axis 1 explained variance: %.1f%%\n", var_explained[1]))
  cat(sprintf("PCoA axis 2 explained variance: %.1f%%\n", var_explained[2]))
  
  p_pcoa <- ggplot(pcoa_df, aes(x = PC1, y = PC2, colour = time)) +
    geom_path(linewidth = 0.4, alpha = 0.7) +
    geom_point(size = 0.8, alpha = 0.6) +
    scale_colour_viridis_c(
      option = "plasma",
      name = "Time step"
    ) +
    labs(
      x = sprintf("PCoA 1 (%.1f%%)", var_explained[1]),
      y = sprintf("PCoA 2 (%.1f%%)", var_explained[2])
    ) +
    theme_bw(base_size = 12) +
    theme(
      aspect.ratio = 1,
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    )
  
  print(p_pcoa)
  
  ggsave(
    filename = sprintf("%s/pcoa_bray_curtis_%s.pdf", save.dir,nam[i]),
    plot = p_pcoa,
    width = 6,
    height = 5
  )
  
}
