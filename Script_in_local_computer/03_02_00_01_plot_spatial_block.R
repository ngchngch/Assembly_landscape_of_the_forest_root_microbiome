setwd("/Users/ngch/Desktop/Now_using/data/sugaFB_ELA/revise_251212/output_spacom/result_for_revise/analysis_in_local_computer")
library(ggplot2)
library(ggspatial)
#################
dir_03_02_00 <- "Output_supercomputer/03_02_00_prep_spatial_block"

##################
save.dir <- "Output/03_02_00_01_plot_spatial_block"
dir.create(save.dir, showWarnings = FALSE)

info <- readRDS(sprintf("%s/comp_sample_info_plant2_spblock.rds",dir_03_02_00))
#info <- readRDS("/Users/ngch/Desktop/Now_using/data/sugaFB_ELA/revise_251212/output_spacom/result_for_revise/03_02_00_prep_spatial_block/comp_sample_info_plant2_spblock.rds")
  
 g <- ggplot(info, aes(x = y_m, y = x_m,color=as.factor(sp_block))) +
    geom_point(size = 3) +
    
    # 方位（北矢印）
    annotation_north_arrow(
      location = "tr",          # "tl", "tr", "bl", "br"
      which_north = "true",
      style = north_arrow_fancy_orienteering,
      height = unit(1.2, "cm"),
      width  = unit(1.2, "cm"),
      pad_x = unit(0.3, "cm"),
      pad_y = unit(0.3, "cm")
    ) +
    
    # 縮尺
    annotation_scale(
      location = "br",
      width_hint = 0.3,
      pad_x = unit(0.3, "cm"),
      pad_y = unit(0.3, "cm")
    ) +
    
    coord_equal() +
    labs(color="Spatial block")+
    theme_bw()+
    theme(axis.text=element_blank(),
          aspect.ratio = 1,
          axis.title=element_blank())+
    scale_color_manual(values=palette(value = "Okabe-Ito"))
  
 ggsave(sprintf("%s/spatial_block.pdf",save.dir),width=6,height=5)
 ggsave(sprintf("%s/spatial_block.png",save.dir),width=6,height=5,dpi=300)
 
 