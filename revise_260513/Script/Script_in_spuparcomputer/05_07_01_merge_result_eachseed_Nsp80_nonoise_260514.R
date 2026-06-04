library(ggpubr)
library(ggplot2)
library(dplyr)



#ELA_prep_dir='02_01_ELA_prep_abundance_threshold'
#n.core=8
#dir_02_03='02_03_summarize_occ_taxa_th'





#########################################################################
#save.dir <- "05_07_01_merge_result_eachseed_Nsp80_nonoise_260514"
dir.create(save.dir)

########################################################################
#read original functions
source("packages/01_1_function.R")
########################################################################
sc <- readRDS(sprintf("%s/top3_seed_across.rds",dir_05_01_00))

res_or <- readRDS(sprintf("%s/ELA_withRA_summary_allseed.rds",dir_05_04))
res_rand <- readRDS(sprintf("%s/summarized_ELA_withRA_4step_Nsp80.rds",dir_05_05))

res_or$ID <- res_or$id
res_or2 <- res_or[,-c(2,3)]

res1 <- merge(res_rand[,c("ID","seed","connectance","ra","z_dland","z_deven")],
             res_or2,
             by=c("ID","seed","connectance","ra"))

res2 <- res1

saveRDS(res2,
        sprintf("%s/summarized_ELA_withRA_basinfreq_NSP80.rds",save.dir))

perc = "perc50"
  
  cr_test <- NULL
  cr_test_jc <- NULL
  for(nsc in 1:nrow(sc)){
     
      cr <- cor.test(na.omit(res2[res2$ra==perc & res2$seed==sc[nsc,2] & res2$connectance==sc[nsc,"connect"],"z_dland"]),
                             na.omit(res2[res2$ra==perc & res2$seed==sc[nsc,2] & res2$connectance==sc[nsc,"connect"],"mBC_diffstart"]),
                     method="kendall")
      
      cr_test <- rbind(cr_test,data.frame(seed=sc[nsc,2],conn=sc[nsc,"connect"],tau=cr$estimate,raw_p=cr$p.value))
      cr_jc <- cor.test(na.omit(res2[res2$ra==perc & res2$seed==sc[nsc,2] & res2$connectance==sc[nsc,"connect"],"z_dland"]),
                     na.omit(res2[res2$ra==perc & res2$seed==sc[nsc,2] & res2$connectance==sc[nsc,"connect"],"mJC_diffstart"]),
                     method="kendall")
      
      cr_test_jc <- rbind(cr_test_jc,data.frame(seed=sc[nsc,2],conn=sc[nsc,"connect"],tau=cr_jc$estimate,raw_p=cr_jc$p.value))
  }
  
  cr_test$p_BH <- p.adjust(cr_test$raw_p,method="BH")
  cr_test_jc$p_BH <- p.adjust(cr_test_jc$raw_p,method="BH")
  saveRDS(list(Bray=cr_test,Jaccard=cr_test_jc),
          sprintf("%s/cor_test_zdL_mBCdiff_NSP80_%s.rds",save.dir,perc))
  
  

  g <- ggplot(na.omit(res2[res2$ra==perc,]),
              aes(x=z_dland,y=mBC_diffstart))+
    stat_cor(
      method = "kendall",
      label.x.npc = "left",
      label.y.npc = "top",
      cor.coef.name = "tau"
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
      method = "kendall",
      label.x.npc = "left",
      label.y.npc = "top",
      cor.coef.name = "tau"
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
