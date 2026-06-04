setwd("/Users/ngch/Desktop/Now_using/data/sugaFB_ELA/revise_260513")
source("packages/fun_gLV_simulation_revise_260514.R")

# %% loading library and functions
library(tidyverse)
library(ggplot2)
library(doParallel)
library(igraph)
library(deSolve)
library(vegan)
library(cluster)
library(patchwork)
library(rELA)
library(philentropy)

save.dir <- "Analysis/06_03_gLV_cluster_numbers_rep15_260522"
dir.create(save.dir)

dir_06_01 <- "06_01_gLV_samplesize_validation_260108" 

Alist <- list.files(dir_06_01,pattern="Amatrix",full.names=TRUE,recursive=TRUE)

rlist <- list.files(dir_06_01,pattern="rparameter",full.names=TRUE,recursive=TRUE)

n.core <- 8
time <- 100
rep <- 100
total_rep <- 15
set.seed(123454567)
seeds <- sapply(strsplit(Alist,"/"),"[",2)
nsp <- sapply(strsplit(Alist,"_"),function(x){as.numeric(gsub("Nsp|.rds","",x[[length(x)]]))})

cluster_summary <- NULL
#,"seed2345","seed3456"
for(tr in 1:total_rep){
for(seed in c("seed1234","seed2345","seed3456")){
  for(sp in c(20,30,40,50)){
    cat("\r","rep:",tr, seed, sp)
    
    A <- readRDS(Alist[seeds==seed & nsp==sp])
    r_set <- readRDS(rlist[seeds==seed & nsp==sp])
    
    run_res <- mclapply(1:rep,function(x){
      mat = gen_dyn(time, dt=1000, A[[1]], r_set[[1]], runif(sp), type="gLV") #ode(y = N0(N), times = seq(1,500), func = lv_fun, parms = list(A=A, r=runif(N)))[,-1]
      end_points <- colMeans(mat[95:100,])  
      return(end_points)
    },mc.cores = n.core) %>% do.call(rbind,.)
   
    ##binarization  
    mat0 <- run_res
    d_bray <- vegdist(mat0,method="bray")
    
    qt <- quantile(mat0, probs = c(0.1,0.9))
    bin_th <- seq(qt[1],qt[2], length.out = 30)
    
    cr <- c()
    for(i in 1:length(bin_th)){
      cat("\n", i, "/", length(bin_th))
      bth <- bin_th[i]
      mat2 <- mat0
      mat2[mat2<=bth] <- 0
      mat2[mat2>0] <- 1
      d_jac <- vegdist(mat2,method="jaccard")
      cr[i] <- cor(d_bray, d_jac,method = "kendall")
    }
    rm(i)
    
    max_cr <- which.max(cr)
    best_bin_th <- bin_th[max_cr[order(max_cr)[1]]]
    
    mat2 <- mat0
    mat2[mat2<=best_bin_th] <- 0
    mat2[mat2>0] <- 1
    
    #clustering ----------------------------------------------------------
    
    ## Jaccard distance
    d_jac <- vegdist(mat2, method = "jaccard")
    
    ## ward.D2 clustering
    hc <- hclust(d_jac, method = "average")
    
    ## silhouette analysis
    max_k <- min(10, nrow(mat2) - 1)  # upper limit of cluster number
    
    sil_width <- rep(NA, max_k)
    
    for(k in 2:max_k){
      
      cl <- cutree(hc, k = k)
      
      sil <- silhouette(cl, d_jac)
      
      sil_width[k] <- mean(sil[, 3])
    }
    
    ## optimal cluster number
    best_k <- which.max(sil_width)
    best_sil <- sil_width[best_k]
    
    ## final cluster assignment
    final_cluster <- cutree(hc, k = best_k)
    
    ## save clustering results
    
    
    ## record summary
    cluster_summary <- rbind(cluster_summary,
                             data.frame(
                               rep=tr,
      seed = seed,
      Nsp = sp,
      sample_size = rep,
      optimal_cluster_n = best_k,
      mean_silhouette = best_sil,
      bin_threshold = best_bin_th
    ))
    
  }
}
}
saveRDS(cluster_summary, file = sprintf("%s/cluster_summary.rds",save.dir))

#3456=20, 2345=30,1234=40,
cluster_summary %>%
  group_by(Nsp,seed) %>%
  summarise(mean_cluster_n = median(optimal_cluster_n))

g <- ggplot(cluster_summary,
            aes(x = as.factor(Nsp),
                y = optimal_cluster_n)) +
  geom_boxplot(
    outlier.color = NA,
    color = "red",
    linewidth = 0.7
  ) +
  geom_jitter(
    width = 0.2,
    height = 0,
    color = "red",
    alpha = 0.5,
    size = 1.5
  ) +
  facet_wrap(~seed, ncol = 1, scales = "free") +
  theme_bw() +
  labs(
    x = "Number of species",
    y = "Optimal number of equilibrium-state clusters"
  ) +
  theme(
    aspect.ratio = 1,
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    strip.text = element_text(size = 12),
    panel.grid.minor = element_blank()
  )

ggsave(
  sprintf("%s/gLV_equibriums.pdf", save.dir),
  plot = g,
  width = 3,
  height = 6
)

#
dir_06_02 <- "Basedata/06_02_gLV_graphics_samplesize_validation"

sm <- readRDS(sprintf("%s/gLV_sample_size_validation_summary.rds",dir_06_02))
## reference median number of clusters for each condition
cluster_median <- cluster_summary %>%
  group_by(Nsp, seed) %>%
  summarise(
    median_cluster_n = median(optimal_cluster_n),
    .groups = "drop"
  )
cluster_median$seed <- substr(cluster_median$seed, 1,5)
g <- ggplot(sm,
            aes(x = as.factor(nsample),
                y = nbasin)) +
  
  ## add median reference line
  geom_hline(
    data = cluster_median,
    aes(yintercept = median_cluster_n),
    color = "red",
    linewidth = 0.8
  ) +
  geom_boxplot(outlier.color = NA) +
  
  geom_jitter(
    width = 0.2,
    height = 0,
    alpha = 0.5
  ) +
  facet_grid(seed ~ Nsp, scales = "free") +
  
  theme_bw() +
  
  labs(
    x = "Number of samples used in energy landscape analysis",
    y = "Number of basins identified in ELA"
  ) +
  
  theme(
    aspect.ratio = 1,
    axis.text.x = element_text(size = 10),
    axis.title = element_text(size = 12),
    axis.text.y = element_text(size = 10),
    strip.text = element_text(size = 12),
    panel.grid.minor = element_blank()
  )

ggsave(
  sprintf("%s/gLV_sample_size_nbasin.pdf", save.dir),
  plot = g,
  width = 16,
  height = 8
)
