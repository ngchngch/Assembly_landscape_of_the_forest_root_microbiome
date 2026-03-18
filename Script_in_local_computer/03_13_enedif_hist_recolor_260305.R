
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
library("ggrepel")

########################################################################
#read original functions
#source("packages/01_1_function.R")


####
########################################################################
dir_03_10 <- "Output/03_10_graphics_states_flow_flow_spl_250508"
ELA_prep_dir <- "Output_supercomputer/02_01_ELA_prep_abundance_threshold"
dir_02_05 <- "Output_supercomputer/02_05_summary_taxa_select"
dir_02_06 <- "Output_supercomputer/02_06_ELA"
dir_02_06_02_02 <- "Output_supercomputer/02_06_02_02_basincomposition_freq_summary"

#setwd("/Users/ngch/Desktop/Now_using/data/sugaFB_ELA/revise_251212/output_spacom/result_for_revise/analysis_in_local_computer")
save.dir <- "Output/03_13_enedif_hist_recolor_260305"
dir.create(save.dir)

########################################################################
# dir_03_10 <- "/Users/ngch/Desktop/Now_using/data/sugaFB_ELA/for_Github/Output/03_10_graphics_states_flow_flow_Spl_250508"
# ELA_prep_dir <- "/Users/ngch/Desktop/Now_using/data/sugaFB_ELA/for_Github/Output_supercomputer/02_01_ELA_prep_abundance_threshold"
# dir_02_05 <- "/Users/ngch/Desktop/Now_using/data/sugaFB_ELA/for_Github/Output_supercomputer/02_05_summary_taxa_select"
# dir_02_06 <- "/Users/ngch/Desktop/Now_using/data/sugaFB_ELA/for_Github/Output_supercomputer/02_06_ELA"
# dir_02_06_03 <- "/Users/ngch/Desktop/Now_using/data/sugaFB_ELA/revise_251212/output_spacom/result_for_revise/02_06_03_summarize_ELA_bootstrap_260304"

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

basin_type <- cbind(rest=c("F_B1","F_B2","F_B3","F_B6","F_B4","F_B5","F_B7"),
                    type=c("EcM","EcM","EcM","Endo","AM","AM","AM"))
bres <- readRDS(sprintf("%s/energy_barrier_boot_%s.rds",dir_02_06_03,fb))
bres$rename_SS2 <- factor(bres$rename_SS,
                               levels=basin_type[,1])
bres$plant2 <- factor(sprintf("*%s*",bres$plant),
                      levels=c("*Pinus*","*Larix*","*Betula*",
                               "*Populus*","*Acer*","*Juglans*"))
bres_type <- readRDS(sprintf("%s/energy_barrier_boot_amongBT_%s.rds",dir_02_06_03,fb))


g <- ggplot(bres[!is.infinite(bres$energy_bound),],aes(x=rename_SS2,y=energy_bound))+
  geom_hline(yintercept = 0,linetype="dashed")+
  geom_violin(aes(fill=rename_SS2),show.legend = FALSE)+
  geom_point(alpha=0.05)+
  facet_wrap(~plant2,scales = "free")+
  labs(y="Energy differences (basin botoms vs the minimum boundaries)",x="Basin bottoms")+
  theme_bw()+
  theme(aspect.ratio = 1,
        strip.text = element_markdown(size=16),
        axis.text.x = element_text(size=12,angle = 45, hjust=1),
        axis.text.y=element_text(size=12),
        axis.title = element_text(size=14))+
  scale_fill_manual(values = sscolor)

ggsave(sprintf("%s/energy_diff_boundary_full_%s.pdf",save.dir,fb),g,width=12,height=8)

g <- ggplot(bres_type[!is.infinite(bres_type$energy_bound),],
            aes(x=rename_SS2,y=energy_bound))+
  geom_hline(yintercept = 0,linetype="dashed")+
  geom_violin(aes(fill=rename_SS2),show.legend = FALSE)+
  geom_point(alpha=0.05)+
  facet_grid(to~sprintf("*%s*",plant),scales = "free")+
  labs(y="Energy differences\n(basin botoms vs the minimum boundaries)",x="Basin bottoms")+
  theme_bw()+
  theme(aspect.ratio = 1,
        strip.text = element_markdown(size=16),
        axis.text.x = element_text(size=12,angle = 45, hjust=1),
        axis.text.y=element_text(size=12),
        axis.title = element_text(size=14))+
  scale_fill_manual(values = sscolor)

g
ggsave(sprintf("%s/energy_diff_amongBT_boundary_full_%s.pdf",save.dir,fb),g,width=12,height=8)

####################
fb <- "Prokaryota"
basin_type <- cbind(rest=c("P_B1","P_B2","P_B3","P_B4","P_B5","P_B6","P_B7","P_B8",
                           "P_B9","P_B10","P_B11","P_B12","P_B13","P_B14","P_B15",
                           "P_B16","P_B17","P_B18","P_B19","P_B20","P_B21","P_B22","P_B23","P_B24"),
                    type=c(rep("Group III",8),
                           rep("Group II",7),
                           rep("Group I",9)
                    ))
rownames(basin_type) <- basin_type[,1]

bres <- readRDS(sprintf("%s/energy_barrier_boot_%s.rds",dir_02_06_03,fb))
bres$rename_SS2 <- factor(bres$rename_SS,
                          levels=basin_type[,1])
bres$plant2 <- factor(sprintf("*%s*",bres$plant),
                      levels=c("*Pinus*","*Larix*","*Betula*",
                               "*Populus*","*Acer*","*Juglans*"))
bres_type <- readRDS(sprintf("%s/energy_barrier_boot_amongBT_%s.rds",dir_02_06_03,fb))


g <- ggplot(bres[!is.infinite(bres$energy_bound),],aes(x=rename_SS2,y=energy_bound))+
  geom_hline(yintercept = 0,linetype="dashed")+
  geom_violin(aes(fill=rename_SS2),show.legend = FALSE)+
  geom_jitter(aes(color=rename_SS2),alpha=0.03,width = 0.2,height = 0,show.legend = FALSE)+
  facet_wrap(~plant2,scales = "free")+
  labs(y="Energy differences (basin botoms vs the minimum boundaries)",x="Basin bottoms")+
  theme_bw()+
  theme(aspect.ratio = 1,
        strip.text = element_markdown(size=16),
        axis.text.x = element_text(size=12,angle = 45, hjust=1),
        axis.text.y=element_text(size=12),
        axis.title = element_text(size=14))+
  scale_fill_manual(values = sscolor)+
  scale_color_manual(values = sscolor)

ggsave(sprintf("%s/energy_diff_boundary_full_%s.pdf",save.dir,fb),g,width=12,height=8)

g <- ggplot(bres_type[!is.infinite(bres_type$energy_bound),],
            aes(x=rename_SS2,y=energy_bound))+
  geom_hline(yintercept = 0,linetype="dashed")+
  geom_violin(aes(fill=rename_SS2),show.legend = FALSE)+
  geom_jitter(aes(color=rename_SS2),alpha=0.03,width = 0.2,height = 0,show.legend = FALSE)+
  facet_grid(to~sprintf("*%s*",plant),scales = "free")+
  labs(y="Energy differences\n(basin botoms vs the minimum boundaries)",x="Basin bottoms")+
  theme_bw()+
  theme(aspect.ratio = 1,
        strip.text = element_markdown(size=16),
        axis.text.x = element_text(size=11,angle = 45, hjust=1),
        axis.text.y=element_text(size=11),
        axis.title = element_text(size=14))+
  scale_fill_manual(values = sscolor)+
  scale_color_manual(values = sscolor)

g
ggsave(sprintf("%s/energy_diff_amongBT_boundary_full_%s.pdf",save.dir,fb),g,width=14,height=8)
