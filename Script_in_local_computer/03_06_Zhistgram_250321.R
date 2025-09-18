
##########################################################################



set.seed(1234)

library(renv)
library(ggstar)
library(ggplot2)
library(ggtext)
library(patchwork)

#detectCores()
n.core <- 8 # not threads!!

dir_03_02_02 <- "Output/03_02_02_Zconv_fixP"
ELA_prep_dir <- "Output_supercomputer/02_01_ELA_prep_abundance_threshold"

#########################################################################
save.dir <- "Output/03_06_Zhistgram_250321"
dir.create(save.dir)

info <- readRDS(sprintf("%s/comp_sample_info_plant2.rds",ELA_prep_dir))
rmat <- readRDS(sprintf("%s/NoCLR_ExpVar_clrRA_target_taxa.rds",ELA_prep_dir))
ramat <- cbind(rmat$Fungi,rmat$Prokaryota)

ocmat <- list(Fungi=readRDS(sprintf("%s/ocmat_remove_M2SD_Fungi_Family.rds",ELA_prep_dir)),
              Prokaryote=readRDS(sprintf("%s/ocmat_remove_M2SD_Prokaryota_Family.rds",ELA_prep_dir)))

tx_f <- readRDS("../../Base_data/Fungi/taxa_list_mod.rds")
tx_b <- as.data.frame(readRDS("../../Base_data/OTU_basedata_set/taxonomy_list_dada2.rds"))
tx_m <- rbind(tx_f[,c("Order","Genus")],tx_b[,c("Order","Genus")])



zval_f_org <- readRDS(sprintf("%s/zvalue_Fungi.rds",dir_03_02_02))
zval_p_org <- readRDS(sprintf("%s/zvalue_Prokaryote.rds",dir_03_02_02))

zval_f_org$plant2 <- factor(sprintf("*%s*",zval_f_org$plant),
                               levels=c("*Pinus*","*Larix*","*Betula*",
                                        "*Populus*","*Acer*","*Juglans*"))
zval_p_org$plant2 <- factor(sprintf("*%s*",zval_p_org$plant),
                               levels=c("*Pinus*","*Larix*","*Betula*",
                                        "*Populus*","*Acer*","*Juglans*"))

pland_BH <- p.adjust(c(zval_f_org$p_land,
                       zval_p_org$p_land),method = "BH")

zval_f_org$p_land_BH <- pland_BH[1:nrow(zval_f_org)]

zval_p_org$p_land_BH <- pland_BH[nrow(zval_f_org)+1:nrow(zval_p_org)]

peven_BH <- p.adjust(c(zval_f_org$p_even,
                       zval_p_org$p_even),method = "BH")

zval_f_org$p_even_BH <- peven_BH[1:nrow(zval_f_org)]

zval_p_org$p_even_BH <- peven_BH[nrow(zval_f_org)+1:nrow(zval_p_org)]


zval_f_org$landscape <- ifelse(zval_f_org$p_land_BH<0.05,
                        "*P* (FDR) < 0.05",
                        "*N.S.*")
zval_f_org$landscape <-factor(zval_f_org$landscape,
                              levels = c("*P* (FDR) < 0.05","*N.S.*"))
zval_p_org$landscape <- ifelse(zval_p_org$p_land_BH<0.05,
                               "*P* (FDR) < 0.05",
                               "*N.S.*")
zval_p_org$landscape <-factor(zval_p_org$landscape,
                              levels = c("*P* (FDR) < 0.05","*N.S.*"))

zval_f_org$evenness <- ifelse(zval_f_org$p_even_BH < 0.025,
                              "*P* (FDR) < 0.05",
                              "*N.S.*")
zval_f_org$evenness <-factor(zval_f_org$evenness,
                             levels = c("*P* (FDR) < 0.05","*N.S.*"))

zval_p_org$evenness <- ifelse(zval_p_org$p_even_BH < 0.025,
                              "*P* (FDR) < 0.05",
                              "*N.S.*")

zval_p_org$evenness <-factor(zval_p_org$evenness,
                             levels = c("*P* (FDR) < 0.05","*N.S.*"))

zval_f_org$ra2 <- ifelse(zval_f_org$ra=="perc25",
                         "25%",
                         ifelse(zval_f_org$ra=="median",
                                "50%",
                                "75%"))

zval_p_org$ra2 <- ifelse(zval_p_org$ra=="perc25",
                         "25%",
                         ifelse(zval_p_org$ra=="median",
                                "50%",
                                "75%"))

for(ra in c("25%","50%","75%")){
  zvalf <- zval_f_org[zval_f_org$ra2==ra,]
  zvalf$z <- "\u0394*S*"
  zvalf$z2 <- "\u0394*H*"
  g <- ggplot(zvalf,
              aes(x=z_land))+
    geom_histogram(aes(fill=landscape,color=landscape),linewidth=0.3)+
    facet_grid(plant2~z,scales = "free_y")+
    labs(x="Z-score",
         y="Frequency",
         fill="Stability landscape\nre-organization",
         color="Stability landscape\nre-organization")+
    theme_bw()+
    theme(aspect.ratio = 1,
          axis.title = element_markdown(size = 13),
          axis.text = element_text(size = 11),
          legend.position="none",
          strip.text.x = element_markdown(size=13),
          strip.text.y = element_blank(),
          strip.background.y = element_blank(),
          legend.text = element_markdown(size=11))+
    scale_fill_manual(values=c("red","gray80"))+
    scale_color_manual(values=c("red","black"))
  g  
  
  
  g2 <- ggplot(zvalf,
               aes(x=z_even))+
    geom_histogram(aes(fill=evenness,color=evenness),linewidth=0.3)+
    facet_grid(plant2~z2,scales = "free_y")+
    labs(x="Z-score",
         #y="Frequency",
         fill="Stochasticity increase/decrease",
         color="Stochasticity increase/decrease")+
    theme_bw()+
    theme(aspect.ratio = 1,
          axis.title.x = element_markdown(size = 13),
          axis.title.y = element_blank(),
          axis.text = element_text(size = 11),
          legend.title = element_blank(),
          strip.text.y = element_markdown(size = 13),
          strip.text.x = element_markdown(size=13),
          legend.text = element_markdown(size=11))+
    scale_fill_manual(values=c("red","gray80"))+
    scale_color_manual(values=c("red","black"))
    
  
  g_merge <-  g+g2+plot_layout(guides="collect")
  
  ggsave(sprintf("%s/Zhist_Fungi_%s.pdf",
                 save.dir,gsub("%","",ra)),
         g_merge,
         width=5,height=8,
         device = cairo_pdf)
  
  
  zvalf <- zval_p_org[zval_p_org$ra2==ra,]
  zvalf$z <- "\u0394*S*"
  zvalf$z2 <- "\u0394*H*"
  g <- ggplot(zvalf,
              aes(x=z_land))+
    geom_histogram(aes(fill=landscape,color=landscape),linewidth=0.3)+
    facet_grid(plant2~z,scales = "free_y")+
    labs(x="Z-score",
         y="Frequency",
         fill="Stability landscape\nre-organization",
         color="Stability landscape\nre-organization")+
    theme_bw()+
    theme(aspect.ratio = 1,
          axis.title = element_markdown(size = 13),
          axis.text = element_text(size = 11),
          legend.position="none",
          strip.text.x = element_markdown(size=13),
          strip.text.y = element_blank(),
          strip.background.y = element_blank(),
          legend.text = element_markdown(size=11))+
    scale_fill_manual(values=c("red","gray80"))+
    scale_color_manual(values=c("red","black"))
  g  
  
  
  g2 <- ggplot(zvalf,
               aes(x=z_even))+
    geom_histogram(aes(fill=evenness,color=evenness),linewidth=0.3)+
    facet_grid(plant2~z2,scales = "free_y")+
    labs(x="Z-score",
         #y="Frequency",
         fill="Stochasticity increase/decrease",
         color="Stochasticity increase/decrease")+
    theme_bw()+
    theme(aspect.ratio = 1,
          axis.title.x = element_markdown(size = 13),
          axis.title.y = element_blank(),
          axis.text = element_text(size = 11),
          legend.title = element_blank(),
          strip.text.y = element_markdown(size = 13),
          strip.text.x = element_markdown(size=13),
          legend.text = element_markdown(size=11))+
    scale_fill_manual(values=c("red","gray80"))+
    scale_color_manual(values=c("red","black"))
  
  
  g_merge <-  g+g2+plot_layout(guides="collect")
  
  ggsave(sprintf("%s/Zhist_Prokaryotes_%s.pdf",
                 save.dir,gsub("%","",ra)),
         g_merge,
         width=5,height=8,
         device = cairo_pdf)
  
}
  