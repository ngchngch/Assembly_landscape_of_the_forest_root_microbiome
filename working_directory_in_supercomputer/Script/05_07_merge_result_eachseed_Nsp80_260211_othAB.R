library(ggpubr)
library(ggplot2)
library(dplyr)



#ELA_prep_dir='02_01_ELA_prep_abundance_threshold'
#n.core=8
#dir_02_03='02_03_summarize_occ_taxa_th'





#########################################################################
#save.dir <- "05_07_merge_result_eachseed_Nsp80_260128"
dir.create(save.dir)

########################################################################
#read original functions
source("packages/01_1_function.R")
########################################################################
sc <- readRDS(sprintf("%s/top3_seed_across.rds",dir_05_01_00))

res_or <- readRDS(sprintf("%s/ELA_withRA_summary_allseed.rds",dir_05_04))
res_rand <- readRDS(sprintf("%s/summarized_ELA_withRA_4step_Nsp80.rds",dir_05_05))
basin_freq <- do.call(rbind,lapply(list.files(dir_05_06,"ELA_evenness_wiwoi_sp",recursive=TRUE,full.names=TRUE),readRDS))

res_or$ID <- res_or$id
res_or2 <- res_or[,-c(2,3)]

basin_freq$ID <- basin_freq$id
basin_freq$connectance <- as.numeric(gsub("conn","",basin_freq$conn))

res1 <- merge(res_rand[,c("ID","seed","connectance","ra","z_dland","z_deven")],
             res_or2,
             by=c("ID","seed","connectance","ra"))

res2 <- merge(res1,
              basin_freq[,c("ID","seed","connectance","d_nbasin","d_shannon")],
              by=c("ID","seed","connectance"))

saveRDS(res2,
        sprintf("%s/summarized_ELA_withRA_basinfreq_NSP80.rds",save.dir))

for(perc in c("perc15","perc35","perc65","perc85")){
  g <- ggplot(na.omit(res2[res2$ra==perc,]),
              aes(x=z_dland,y=degree_out_abs))+
    stat_cor(
      method = "spearman",
      label.x.npc = "left",
      label.y.npc = "top",
      cor.coef.name = "rho"
    ) +
    geom_point() +
    labs(x="*delta*topography",
         y="weighted out degree")+
    facet_wrap(seed~connectance,scales="free") +
    theme_bw()+
    theme(aspect.ratio = 1,
          text = element_text(size=14),
    )
  
  ggsave(sprintf("%s/zdL_wdegree_out_abs_NSP80_%s.pdf",save.dir,perc),
         plot=g,
         width=12,
         height=8)

g <- ggplot(na.omit(res2[res2$ra==perc,]),
              aes(x=z_dland,y=mAB*degree_out_abs,color=mAB))+
    stat_cor(
      method = "spearman",
      label.x.npc = "left",
      label.y.npc = "top",
      cor.coef.name = "rho"
    ) +
    geom_point() +
    labs(x="*delta*topography",
         y="weighted out degree*mAB")+
    facet_wrap(seed~connectance,scales="free") +
    theme_bw()+
    theme(aspect.ratio = 1,
          text = element_text(size=14),
    )
  
  ggsave(sprintf("%s/zdL_mABxwdegree_out_abs_NSP80_%s.pdf",save.dir,perc),
         plot=g,
         width=12,
         height=8)

  g <- ggplot(na.omit(res2[res2$ra==perc,]),
              aes(x=z_dland,y=mBC_diffstart))+
    stat_cor(
      method = "spearman",
      label.x.npc = "left",
      label.y.npc = "top",
      cor.coef.name = "rho"
    ) +
    geom_point() +
    labs(x="*delta*topography",
         y="Mean community distance\n(Bray-Crutis distance with/without a species)")+
    facet_wrap(seed~connectance,scales="free") +
    theme_bw()+
    theme(aspect.ratio = 1,
          text = element_text(size=14),
    )
  
  ggsave(sprintf("%s/zdL_mBCdiff_NSP80_%s.pdf",save.dir,perc),
         plot=g,
         width=12,
         height=8)
  
  g <- ggplot(na.omit(res2[res2$ra==perc,]),
              aes(x=z_dland,y=mJC_diffstart))+
    stat_cor(
      method = "spearman",
      label.x.npc = "left",
      label.y.npc = "top",
      cor.coef.name = "rho"
    ) +
    geom_point() +
    labs(x="*delta*topography",
         y="Mean community distance\n(Jaccard distance with/without a species)")+
    facet_wrap(seed~connectance,scales="free") +
    theme_bw()+
    theme(aspect.ratio = 1,
          text = element_text(size=14),
    )
  
  ggsave(sprintf("%s/zdL_mJCdiff_NSP80_%s.pdf",save.dir,perc),
         plot=g,
         width=12,
         height=8)
  
  g <- ggplot(na.omit(res2[res2$ra==perc,]),
              aes(x=z_dland,y=mBC_diffstart))+
    stat_cor(
      method = "spearman",
      label.x.npc = "left",
      label.y.npc = "top",
      cor.coef.name = "rho"
    ) +
    geom_point() +
    labs(x="*delta*topography",
         y="Mean community distance\n(Bray-Crutis distance with/without a species)")+
    facet_wrap(seed~connectance,scales="free") +
    theme_bw()+
    theme(aspect.ratio = 1,
          text = element_text(size=14),
    )
  
  ggsave(sprintf("%s/zdL_mBC_ds_l2r_NSP80_%s.pdf",save.dir,perc),
         plot=g,
         width=12,
         height=8)
  
  
  
  g <- ggplot(na.omit(res2[res2$ra==perc,]),
              aes(x=z_deven,y=shannon_cls_bray))+
    stat_cor(
      method = "spearman",
      label.x.npc = "left",
      label.y.npc = "top",
      cor.coef.name = "rho"
    ) +
    geom_point() +
    facet_wrap(seed~connectance,scales="free") +
    theme_bw()+
    theme(aspect.ratio = 1,
          text = element_text(size=14),
    )
  
  ggsave(sprintf("%s/zdE_shannon_cls_BC_NSP80_%s.pdf",save.dir,perc),
         plot=g,
         width=12,
         height=8)
  
  g <- ggplot(na.omit(res2[res2$ra==perc,]),
              aes(x=z_deven,y=shannon_cls_jac))+
    stat_cor(
      method = "spearman",
      label.x.npc = "left",
      label.y.npc = "top",
      cor.coef.name = "rho"
    ) +
    geom_point() +
    facet_wrap(seed~connectance,scales="free") +
    theme_bw()+
    theme(aspect.ratio = 1,
          text = element_text(size=14),
    )
  
  ggsave(sprintf("%s/zdE_shannon_cls_JC_NSP80_%s.pdf",save.dir,perc),
         plot=g,
         width=12,
         height=8)
      
  g <- ggplot(na.omit(res2[res2$ra==perc,]),
              aes(x=z_deven,y=d_shannon))+
    stat_cor(
      method = "spearman",
      label.x.npc = "left",
      label.y.npc = "top",
      cor.coef.name = "rho"
    ) +
    geom_point() +
    facet_wrap(seed~connectance,scales="free") +
    theme_bw()+
    theme(aspect.ratio = 1,
          text = element_text(size=14),
    )
  
  ggsave(sprintf("%s/zdE_delta_basin_shannon_NSP80_%s.pdf",save.dir,perc),
         plot=g,
         width=12,
         height=8)
  
  for(r in 1:15){
    seed <- sc[r,2]
    conn0 <- sc[r,1]
    
    png(sprintf("%s/heatmap_ELA_seed%s_conn%s_NSP80_%s.png",save.dir,seed,conn0,perc),
        width=1800,
        height=1800)
    pairs(res2[res2$seed == seed & res2$connectance==conn0 & res2$ra==perc,
               c("z_dland","z_deven","d_land","d_even","mBC_diffstart","mJC_diffstart",
                 "mAB","r","degree_out","degree_out_abs","degree_in",
                 "betweenness",   "closenes",   "pagerank",
                 "delta_lambda_max","delta_modularity","l2r_id_BC","l2r_id_JC")])
    dev.off()
  }
}
