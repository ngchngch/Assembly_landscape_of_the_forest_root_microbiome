library(ggpubr)
library(ggplot2)
library(dplyr)



#ELA_prep_dir='02_01_ELA_prep_abundance_threshold'
#n.core=8
#dir_02_03='02_03_summarize_occ_taxa_th'
#dir_07_05 <- "/Users/ngch/Desktop/Now_using/data/sugaFB_ELA/revise_260513/07_05_Zconv_ELA_withRA_4step_Nsp80_multistable_260519"




#########################################################################
#save.dir <- "07_06_merge_result_eachseed_Nsp80_multistable"
dir.create(save.dir)

########################################################################
#read original functions
source("packages/01_1_function.R")
########################################################################

keylist <- list.files(dir_07_02,pattern="keystoneness_summary",full.names = TRUE,recursive = TRUE)
  
tops <- strsplit(keylist,"/") %>%
  lapply(function(x)x[length(x)-1]) %>%
  unlist() %>%
  strsplit("_") %>%
  lapply(function(x) as.numeric(gsub("top","",x[1]))) %>%
  unlist()


keystone_summary <- lapply(1:3,function(x){cbind(top=tops[x],readRDS(keylist[x]))}) %>%
  do.call(rbind,.)
res_rand <- readRDS(sprintf("%s/summarized_ELA_withRA_4step_Nsp80.rds",dir_07_05))

keystone_summary$ID <- keystone_summary$id

res1 <- merge(res_rand[,c("ID","ra","top","z_dland","d_land",
                          "z_deven","d_even")],
              keystone_summary,
             by=c("ID","top"))

res2 <- res1[which(!is.na(res1$z_dland)),]

perc = "perc50"
  
  cr_test <- NULL
  cr_test_jc <- NULL
  for(nsc in 1:length(tops)){
     
      cr <- cor.test(na.omit(res2[res2$ra==perc & res2$top==nsc,"z_dland"]),
                             na.omit(res2[res2$ra==perc & res2$top==nsc,"mBC_diffstart"]),
                     method="kendall")
      
      cr_test <- rbind(cr_test,data.frame(top=nsc,tau=cr$estimate,raw_p=cr$p.value))
      cr_jc <- cor.test(na.omit(res2[res2$ra==perc & res2$top==nsc,"z_dland"]),
                     na.omit(res2[res2$ra==perc & res2$top==nsc,"mJC_diffstart"]),
                     method="kendall")
      
      cr_test_jc <- rbind(cr_test_jc,data.frame(top=nsc,tau=cr_jc$estimate,raw_p=cr_jc$p.value))
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
    facet_wrap(~top,scales="free") +
    theme_bw()+
    theme(aspect.ratio = 1,
          text = element_text(size=14),
    )
  
  ggsave(sprintf("%s/zdL_mBCdiff_NSP80_%s.pdf",save.dir,perc),
         plot=g,
         width=10,
         height=5)
  
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
    facet_wrap(~top,scales="free") +
    theme_bw()+
    theme(aspect.ratio = 1,
          text = element_text(size=14),
    )
  
  ggsave(sprintf("%s/zdL_mJCdiff_NSP80_%s.pdf",save.dir,perc),
         plot=g,
         width=10,
         height=5)
