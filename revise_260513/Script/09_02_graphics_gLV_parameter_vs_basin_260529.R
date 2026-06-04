library(tidyverse)
library(ggpubr)
library(FSA)
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

save.dir <- "Analysis/09_02_graphics_gLV_parameter_vs_basin"
dir.create(save.dir)
############################################################
# load summary
############################################################
dir_09_01 <- "Output_spacon/09_01_gLV_parameter_vs_basin"


cluster_summary <- read.csv(
  sprintf("%s/cluster_summary.csv", dir_09_01)
)

############################################################
# factor化
############################################################

cluster_summary$connectance <-
  factor(cluster_summary$connectance)

############################################################
# 1. connectance vs n_basin
############################################################

p_conn <- ggplot(
  cluster_summary,
  aes(x = connectance, y = n_basin,color=connectance)
) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(
    width = 0.15,
    alpha = 0.6,
    size = 2
  ) +
  theme_bw(base_size = 14) +
  labs(
    x = "Connectance",
    y = "Number of basins"
  )

ggsave(
  sprintf("%s/connectance_vs_basin.pdf", save.dir),
  p_conn,
  width = 5,
  height = 4
)

###########################################################
## statistics : connectance
############################################################
pairwise_wilcox_conn <- pairwise.wilcox.test(  x = cluster_summary$n_basin,
                                               g = cluster_summary$connectance,
                                               p.adjust.method = "fdr")
print(pairwise_wilcox_conn)

############################################################
# 2. mean interaction strength vs n_basin
############################################################

p_strength <- ggplot(
  cluster_summary,
  aes(
    x = mean_interaction_strength,
    y = n_basin
  )
) +
  geom_point(aes(color=connectance),alpha = 0.7) +
  theme_bw(base_size = 14) +
  theme(aspect.ratio = 1)+
  labs(
    x = "Mean interaction strength",
    y = "Number of basins"
  )

ggsave(
  sprintf(
    "%s/mean_strength_vs_basin.pdf",
    save.dir
  ),
  p_strength,
  width = 5,
  height = 4
)

############################################################
# statistics : interaction strength
############################################################

cor_strength <- cor.test(
  cluster_summary$mean_interaction_strength,
  cluster_summary$n_basin,
  method = "kendal"
)

print(cor_strength)

############################################################
# save statistics
############################################################

capture.output(
  kw_conn,
  file = sprintf(
    "%s/stat_connectance_kw.txt",
    save.dir
  )
)

capture.output(
  cor_strength,
  file = sprintf(
    "%s/stat_strength_spearman.txt",
    save.dir
  )
)
