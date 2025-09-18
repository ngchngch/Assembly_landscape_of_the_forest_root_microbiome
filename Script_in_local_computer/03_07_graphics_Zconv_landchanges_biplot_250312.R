
##########################################################################
#install.packages("/Volumes/8TBHDD_NGCH/sugadaira_bacteria_2023/240801_SSchange_randamize/packages/rELA.v0.51.tar")


set.seed(1234)

library(renv)
library(ggstar)
library(ggplot2)
library(ggtext)
library(ggrepel)
library(patchwork)

#detectCores()
n.core <- 8 # not threads!!

dir_03_02_02 <- "Output/03_02_02_Zconv_fixP"
source("packages/01_1_function.R")
#########################################################################
save.dir <- "Output/03_07_graphics_Zconv_landchanges_each_biplot_250312"
dir.create(save.dir)

tx_f <- readRDS("Base_data/Fungi/taxa_list_mod.rds")
tx_p <- as.data.frame(readRDS("Base_data/OTU_basedata_set/taxonomy_list_dada2.rds"))
tx_m <- rbind(tx_f[,c("Order","Genus")],tx_p[,c("Order","Genus")])

zval_f_org <- readRDS(sprintf("%s/zvalue_Fungi.rds",dir_03_02_02))
zval_p_org <- readRDS(sprintf("%s/zvalue_Prokaryote.rds",dir_03_02_02))

pland_BH <- p.adjust(c(zval_f_org$p_land,zval_p_org$p_land),method = "BH")
zval_f_org$p_land_BH <- pland_BH[1:length(zval_f_org$p_land)]
zval_p_org$p_land_BH <- pland_BH[(length(zval_f_org$p_land)+1):length(pland_BH)]

peven_BH <- p.adjust(c(zval_f_org$p_even,zval_p_org$p_even),method = "BH")
zval_f_org$p_even_BH <- peven_BH[1:length(zval_f_org$p_even)]
zval_p_org$p_even_BH <- peven_BH[(length(zval_f_org$p_even)+1):length(peven_BH)]

zval_f_org$target <- ifelse(zval_f_org$Taxa %in% tx_f[,"Genus"],
                            paste("Fng",sprintf("*%s*",zval_f_org$Taxa),sep="_"),
                            paste("Prok",sprintf("*%s*",zval_f_org$Taxa),sep="_"))

zval_f_org$target_k <- ifelse(zval_f_org$Taxa %in% tx_f[,"Genus"],
                              "Fungi",
                              "Prokaryote")

zval_f_org$target_p <- ifelse(zval_f_org$Taxa %in% tx_f[,"Genus"],
                              paste("Fng",tx_f[,"Phylum"][match(zval_f_org$Taxa,tx_f[,"Genus"])],sep="_"),
                              paste("Prok",tx_p[,"Phylum"][match(zval_f_org$Taxa,tx_p[,"Genus"])],sep="_"))

zval_f_org$target_c <- ifelse(zval_f_org$Taxa %in% tx_f[,"Genus"],
                              paste("Fng",tx_f[,"Class"][match(zval_f_org$Taxa,tx_f[,"Genus"])],sep="_"),
                              paste("Prok",tx_p[,"Class"][match(zval_f_org$Taxa,tx_p[,"Genus"])],sep="_"))

zval_f_org$guild <- tx_f[,"Guild"][match(zval_f_org$Taxa,tx_f[,"Genus"])]

write.csv(zval_f_org[which(zval_f_org$guild=="EcMF"),],
          file = sprintf("%s/zvalue_EcMF_fungi.csv",save.dir),
          row.names = FALSE, quote = FALSE)

zval_p_org$target <- ifelse(zval_p_org$Taxa %in% tx_f[,"Genus"],
                            paste("Fng",sprintf("*%s*",zval_p_org$Taxa),sep="_"),
                            paste("Prok",sprintf("*%s*",zval_p_org$Taxa),sep="_"))

zval_p_org$target_k <- ifelse(zval_p_org$Taxa %in% tx_f[,"Genus"],
                              "Fungi",
                              "Prokaryote")

zval_p_org$target_p <- ifelse(zval_p_org$Taxa %in% tx_f[,"Genus"],
                              paste("Fng",tx_f[,"Phylum"][match(zval_p_org$Taxa,tx_f[,"Genus"])],sep="_"),
                              paste("Prok",tx_p[,"Phylum"][match(zval_p_org$Taxa,tx_p[,"Genus"])],sep="_"))

zval_p_org$target_c <- ifelse(zval_p_org$Taxa %in% tx_f[,"Genus"],
                              paste("Fng",tx_f[,"Class"][match(zval_p_org$Taxa,tx_f[,"Genus"])],sep="_"),
                              paste("Prok",tx_p[,"Class"][match(zval_p_org$Taxa,tx_p[,"Genus"])],sep="_"))

zval_p_org$guild <- tx_f[,"Guild"][match(zval_p_org$Taxa,tx_f[,"Genus"])]

write.csv(zval_f_org[which(zval_f_org$guild=="EcMF"),],
          file = sprintf("%s/zvalue_EcMF_prok.csv",save.dir),
          row.names = FALSE, quote = FALSE)

# Prok_Proteobacteria       Prok_Chloroflexi   Prok_Acidobacteriota  Prok_Actinobacteriota 
# "#00005B"                    "#E8FFFF"        "#AB005D"              "#FF8DD7" 
# Prok_Verrucomicrobiota      Fng_Basidiomycota       Prok_Myxococcota   Prok_Planctomycetota 
# "#008203"              "#242424"              "#00D087"              "#E425A6" 
# Fng_Ascomycota    Prok_Armatimonadota        Prok_Firmicutes  Prok_Bdellovibrionota 
# "#FFFF00"              "#0000E1"              "#049CC6"              "#3E5025" 
# Prok_Bacteroidota   Prok_Gemmatimonadota  Fng_Mortierellomycota      Prok_Nitrospirota 
# "#FAAD3B"              "#AD5700"              "#31A800"              "#4561FE" 
# Fng_Kickxellomycota      Fng_Rozellomycota     Prok_Spirochaetota   Prok_Patescibacteria 
# "#F06A00"              "#18E6FF"              "#69323D"              "#00657A" 

col_phyl <- c("#2DDBC4","#E8FFFF","#AB005D","#FF8DD7",
              "#FC2200","#F06A00","#5A7966","#0D41F5",
              "#FFFF00","#008203","#049CC6","#FFFFD6",
              "#F546FF","#AD5700","#C393EB","#6776B3",
              "#00005B","#18E6FF","#69323D","#B8A55F")
names(col_phyl) <- unique(c(zval_f_org$target_p,zval_p_org$target_p))

zval_l <- rbind(cbind(fb="Fungi",zval_f_org),
                cbind(fb="Prokaryote",zval_p_org))
#coloring taxa selection

color_taxa <- c()

for(l in 1:length(unique(zval_l$ra))){#l <- 1
  zval1 <- na.omit(zval_l[zval_l$ra == unique(zval_l$ra)[l],])
  for(pl in 1:length(unique(zval1$plant))){
    for(f2 in 1:length(unique(zval1$target_k))){
      for(f in 1:length(unique(zval1$fb))){#pl <- 1;f <- 1
        zval2 <- zval1[which(zval1$fb == unique(zval1$fb)[f] &
                               zval1$plant == unique(zval1$plant)[pl] &
                               zval1$target_k == unique(zval1$target_k)[f2]),]
        
        top <- zval2[order(zval2$z_land,decreasing = TRUE)[1:5],"target"]
        sig <- zval2[zval2$p_land_BH<0.05,"target"]
        
        color_taxa <- c(color_taxa,
                        intersect(top,sig))
      }
    }
    
  }  
}

pick_taxa <- c(unique(color_taxa),"Fng_Others","Prok_Others")

taxcol_genus <- cbind(c(pick_taxa[order(pick_taxa)]),
      c(color[1:sum(grepl("Fng_",unique(color_taxa)))],
        "grey40",
        color[sum(grepl("Fng_",unique(color_taxa)))+1:sum(grepl("Prok_",unique(color_taxa)))],
        "grey80"))

rownames(taxcol_genus) <- taxcol_genus[,1]

saveRDS(taxcol_genus,sprintf("%s/taxcol_genus.rds",save.dir))

le <- zval_l
  lef <- le[!is.na(le$p_land),]
  lef$z_land[is.na(lef$z_land)] <- 0
  lef$z_even[is.na(lef$z_even)] <- 0
  
  lef$landscape <- ifelse(lef$p_land_BH<0.05,
                          ifelse(lef$p_even_BH<0.025,
                                 ifelse(lef$z_even>0,
                                        "with \u0394*H* increase",
                                        "with \u0394*H* reduction"),
                                 "without \u0394*H* change"),
                          "No re-organization")
  
  
  
  lef$target2 <- ifelse(lef$target %in% rownames(taxcol_genus),
                        lef$target,
                        ifelse(grepl("Fng_",lef$target),"Fng_Others",
                               "Prok_Others"))
  
  lef$plant2 <- sprintf("*%s*",lef$plant)
  
  lef$tag_ord <- factor(lef$target2,
                        levels = rownames(taxcol_genus))
  
  #top 2 selct for plot annotation
  
  
  #topn <- 3
  lef$Top <- FALSE
  topmat <- NULL
  topmat_e <- NULL
  
  for(i in 1:length(unique(lef$fb))){#i<- 1
    for(j in 1:length(unique(lef$plant))){#j<-1
      for(ra in 1:length(unique(lef$ra))){#ra <- 1
        #top2 landscape change
        lef2 <- lef[which(lef$fb==unique(lef$fb)[i]&lef$plant==unique(lef$plant)[j]&lef$ra==unique(lef$ra)[ra]&lef$p_land_BH<0.05),]
        
        # lef2 <- NULL
        # for(k in 1:length(unique(lef$Taxa))){#k <- 1
        #   lef1 <- sel_lef[which(sel_lef$Taxa==unique(sel_lef$Taxa)[k]),]
        #   lef2 <- rbind(lef2,lef1[which.max(lef1$z_land),])
        # }
        # 
        th <- 2
        
        if(nrow(lef2)>th){
          sel_tag <-   lef2[order(lef2$z_land,decreasing = TRUE)[1:th],]
          
          topmat <-   rbind(topmat,lef2[order(lef2$z_land,decreasing = TRUE)[1:th],])
          
        }else{
          sel_tag <-   lef2
          
          topmat <-   rbind(topmat,sel_tag)
        }
        
        for(l in 1:nrow(sel_tag)){
          lef[which(lef$fb==sel_tag$fb[l]&
                      lef$plant==sel_tag$plant[l]&
                      lef$target%in%sel_tag$target[l]&
                      lef$ra==sel_tag$ra[l]),"Top"] <- TRUE
          
        }  
        
        
        #top2 evenness increase
        sel_lef2 <- lef2[lef2$p_even_BH<0.025,]
        pos_lef <- sel_lef2[sel_lef2$z_even>0,]
        neg_lef <- sel_lef2[sel_lef2$z_even<0,]
       
        if(nrow(pos_lef)>th & nrow(neg_lef)>th){
          sel_tag2 <- rbind(pos_lef[order(pos_lef$z_even,decreasing = TRUE)[1:th],],
                            neg_lef[order(neg_lef$z_even)[1:th],])
          topmat_e <- rbind(topmat_e,sel_tag2)
        }else{
          if(nrow(pos_lef)>th){
            sel_tag2 <- rbind(pos_lef[order(pos_lef$z_even,decreasing = TRUE)[1:th],],
                              neg_lef)
            topmat_e <- rbind(topmat_e,sel_tag2)
          }else{
            if(nrow(neg_lef)>th){
              sel_tag2 <- rbind(pos_lef,
                                neg_lef[order(neg_lef$z_even)[1:th],])
              topmat_e <- rbind(topmat_e,sel_tag2)
            }else{
              sel_tag2 <- rbind(pos_lef,
                                neg_lef)
              topmat_e <- rbind(topmat_e,sel_tag2)
          }
        }}
          
         #topmat_e_upper <- rbind(topmat_e_upper,lef2[order(lef2$z_even,decreasing = TRUE)[1:th],])
        
        # }else{
        #   sel_tag2 <-lef2
        #   topmat_e_upper <- rbind(topmat_e_upper,lef2)
        #   
        # }
        
        for(l in 1:nrow(sel_tag2)){
          lef[which(lef$fb==sel_tag2$fb[l]&
                      lef$plant==sel_tag2$plant[l]&
                      lef$target%in%sel_tag2$target[l]&
                      lef$ra==sel_tag2$ra[l]),"Top"] <- TRUE
          
        }  
      }}}
  
  
  
  lef$target_label <- ifelse(lef$Taxa %in% tx_f[,"Genus"],
                             paste("Fng",sprintf("~italic('%s')",lef$Taxa),sep="_"),
                             paste("Prok",sprintf("~italic('%s')",lef$Taxa),sep="_"))
  
  
  sum(lef$Top)
  write.csv(topmat,sprintf("%s/delta_landscape_topmat.csv",save.dir))
  write.csv(topmat_e,sprintf("%s/delta_landscape_topmat_e.csv",save.dir))
  saveRDS(unique(c(topmat$Taxa,topmat_e$Taxa)),
          sprintf("%s/Taxa_for_flow_diagram.rds",save.dir))
  ##############
  ###################
  lef <- lef[order(match(lef$ra, c("perc25","median","perc75"))),]
  fb <- "Fungi"
  
  
  for(ra in c("perc25","median","perc75")){
    g1 <- list(NULL)
    a <- 0 
  for(plants in c("Acer","Betula","Pinus","Juglans","Populus","Larix")){
    a <- a+1
      g1[[a]] <- ggplot(lef[which(lef$fb==fb & lef$plant %in% plants & lef$ra %in% ra),],
                   aes(x=z_land,y=z_even))+
        geom_vline(xintercept = 0,linewidth=0.2,color="gray60")+
        geom_hline(yintercept = 0,linewidth=0.2,color="gray60")+
       # geom_path(aes(group=Taxa),color="gray",
        #          arrow = arrow(length = unit(0.2, "cm")),
         #         linewidth=0.5)+
        geom_star(data=function(x){x[x$landscape=="No re-organization",]},
                  aes(size=landscape,starshape=landscape,fill=target_p))+
        geom_star(data=function(x){x[x$landscape!="No re-organization",]},
                  aes(size=landscape,starshape=landscape,fill=target_p))+
        labs(x="Z-standardized \u0394*S* \n( Fungal community )",
             y="Z-standardized \u0394*H* \n( Fungal community )",
             fill="Genus",
             size="Energy landscape re-organization",
             starshape="Energy landscape re-organization")+
        theme_bw()+
        facet_wrap(~plant2,nrow=1,scales="free")+
        scale_fill_manual(values = col_phyl)+
        scale_size_manual(values = c("No re-organization"=2,
                                     "without \u0394*H* change"=4,
                                     "with \u0394*H* reduction"=4,
                                     "with \u0394*H* increase"=4))+
        scale_starshape_manual(values = c("No re-organization"=15,
                                          "without \u0394*H* change"=28,
                                          "with \u0394*H* reduction"=23,
                                          "with \u0394*H* increase"=11))+
        theme(aspect.ratio = 1,
              legend.position = "none",
              strip.text = element_markdown(size = 15),
              axis.text = element_text(size = 10),
              axis.title = element_markdown(size = 12),
              legend.text = element_markdown(size=9))+
        guides(fill=guide_legend(override.aes = list(starshape=15,size=3),
                                 ncol=2))
      
      
      
       ggsave(sprintf("%s/Z_land_vs_even_change_%s_nolab_%s_%s.pdf",save.dir,fb,plants,ra),g1[[a]],width=4,height=4,
              device = cairo_pdf)
      # 
      ggsave(sprintf("%s/Z_land_vs_even_change_%s_%s_%s.pdf",save.dir,fb,plants,ra),g1[[a]]+
               geom_text_repel(data=function(x){x[x$Top,]},
                               min.segment.length = 0.1,
                               aes(label=target_label),size=3,
                               parse = TRUE),width=6,height=6,
             device = cairo_pdf)
      
      if(a==1){
        ggsave(plot=g_legend(g1[[a]]+theme(legend.position="bottom",legend.direction = "vertical")+
                               guides(fill=guide_legend(override.aes = list(starshape=15,size=3),
                                                        nrow=4))),
               sprintf("%s/legend_Z_land_vs_even_change_%s_%s_%s.pdf",save.dir,fb,plants,ra),width=16,height=4,
               device = cairo_pdf)
        
      }
      
  }
    mg <- g1[[1]] +
      geom_text_repel(data=function(x){x[x$Top,]},
                      min.segment.length = 0.1,
                      aes(label=target_label),size=3,
                      parse = TRUE)+
      g1[[2]] +
      geom_text_repel(data=function(x){x[x$Top,]},
                      min.segment.length = 0.1,
                      aes(label=target_label),size=3,
                      parse = TRUE)+
      g1[[3]] +
      geom_text_repel(data=function(x){x[x$Top,]},
                      min.segment.length = 0.1,
                      aes(label=target_label),size=3,
                      parse = TRUE)+
      g1[[4]] +
      geom_text_repel(data=function(x){x[x$Top,]},
                      min.segment.length = 0.1,
                      aes(label=target_label),size=3,
                      parse = TRUE)+
      g1[[5]] +
      geom_text_repel(data=function(x){x[x$Top,]},
                      min.segment.length = 0.1,
                      aes(label=target_label),size=3,
                      parse = TRUE)+
      g1[[6]] +
      geom_text_repel(data=function(x){x[x$Top,]},
                      min.segment.length = 0.1,
                      aes(label=target_label),size=3,
                      parse = TRUE)+
      plot_layout(ncol=3)
    ggsave(sprintf("%s/Z_land_vs_even_change_%s_%s.pdf",save.dir,fb,ra),mg,width=12,height=8,
           device = cairo_pdf)
      
    mg <- g1[[1]] +theme(axis.title.x = element_blank())+
      g1[[2]] +theme(axis.title = element_blank())+
      g1[[3]] +theme(axis.title = element_blank())+
      g1[[4]] +theme(axis.title = element_blank())+
      g1[[5]] +theme(axis.title.y = element_blank())+
      g1[[6]] +theme(axis.title = element_blank())+
      plot_layout(ncol=3)
    ggsave(sprintf("%s/Z_land_vs_even_change_nolab_%s_%s.pdf",save.dir,fb,ra),mg,width=10,height=7,
           device = cairo_pdf)
  }
  
  
  fb <- "Prokaryote"
  for(ra in c("perc25","median","perc75")){#ra <- "median"
    g1 <- list(NULL)
    a <- 0 
    for(plants in c("Acer","Betula","Pinus","Juglans","Populus","Larix")){#plants <- "Acer"
      lef2 <- lef[which(lef$fb==fb & lef$plant %in% plants  & lef$ra %in% ra),]
      a <- a+1
      g1[[a]] <- ggplot(lef[which(lef$fb==fb & lef$plant %in% plants & lef$ra %in% ra),],
                        aes(x=z_land,y=z_even))+
        geom_vline(xintercept = 0,linewidth=0.2,color="gray60")+
        geom_hline(yintercept = 0,linewidth=0.2,color="gray60")+
        # geom_path(aes(group=Taxa),color="gray",
        #          arrow = arrow(length = unit(0.2, "cm")),
        #         linewidth=0.5)+
        geom_star(data=function(x){x[x$landscape=="No re-organization",]},
                  aes(size=landscape,starshape=landscape,fill=target_p))+
        geom_star(data=function(x){x[x$landscape!="No re-organization",]},
                  aes(size=landscape,starshape=landscape,fill=target_p))+
        labs(x="Z-standardized \u0394*S* \n( Prokaryotic community )",
             y="Z-standardized \u0394*H* \n( Prokaryotic community )",
             fill="Genus",
             size="Energy landscape re-organization",
             starshape="Energy landscape re-organization")+
        theme_bw()+
        facet_wrap(~plant2,nrow=1,scales="free")+
        scale_fill_manual(values = col_phyl)+
        scale_size_manual(values = c("No re-organization"=2,
                                     "without \u0394*H* change"=4,
                                     "with \u0394*H* reduction"=4,
                                     "with \u0394*H* increase"=4))+
        scale_starshape_manual(values = c("No re-organization"=15,
                                          "without \u0394*H* change"=28,
                                          "with \u0394*H* reduction"=23,
                                          "with \u0394*H* increase"=11))+
        theme(aspect.ratio = 1,
              legend.position = "none",
              strip.text = element_markdown(size = 15),
              axis.text = element_text(size = 10),
              axis.title = element_markdown(size = 12),
              legend.text = element_markdown(size=9))+
        guides(fill=guide_legend(override.aes = list(starshape=15,size=3),
                                 ncol=2))
      
      
      
      ggsave(sprintf("%s/Z_land_vs_even_change_%s_nolab_%s_%s.pdf",save.dir,fb,plants,ra),g1[[a]],width=4,height=4,
             device = cairo_pdf)
      # 
      ggsave(sprintf("%s/Z_land_vs_even_change_%s_%s_%s.pdf",save.dir,fb,plants,ra),g1[[a]]+
               geom_text_repel(data=function(x){x[x$Top,]},
                               min.segment.length = 0.1,
                               aes(label=target_label),size=3,
                               parse = TRUE),width=6,height=6,
             device = cairo_pdf)
      
      if(a==1){
        ggsave(plot=g_legend(g1[[a]]+theme(legend.position="bottom")+
                               guides(fill=guide_legend(override.aes = list(starshape=15,size=3),
                                                        nrow=4))),
               sprintf("%s/legend_Z_land_vs_even_change_%s_%s_%s.pdf",save.dir,fb,plants,ra),width=16,height=4,
               device = cairo_pdf)
        
      } 
    }
    mg <- g1[[1]] +
      geom_text_repel(data=function(x){x[x$Top,]},
                      min.segment.length = 0.1,
                      aes(label=target_label),size=3,
                      parse = TRUE)+
      g1[[2]] +
      geom_text_repel(data=function(x){x[x$Top,]},
                      min.segment.length = 0.1,
                      aes(label=target_label),size=3,
                      parse = TRUE)+
      g1[[3]] +
      geom_text_repel(data=function(x){x[x$Top,]},
                      min.segment.length = 0.1,
                      aes(label=target_label),size=3,
                      parse = TRUE)+
      g1[[4]] +
      geom_text_repel(data=function(x){x[x$Top,]},
                      min.segment.length = 0.1,
                      aes(label=target_label),size=3,
                      parse = TRUE)+
      g1[[5]] +
      geom_text_repel(data=function(x){x[x$Top,]},
                      min.segment.length = 0.1,
                      aes(label=target_label),size=3,
                      parse = TRUE)+
      g1[[6]] +
      geom_text_repel(data=function(x){x[x$Top,]},
                      min.segment.length = 0.1,
                      aes(label=target_label),size=3,
                      parse = TRUE)+
      plot_layout(ncol=3)
    ggsave(sprintf("%s/Z_land_vs_even_change_%s_%s.pdf",save.dir,fb,ra),mg,width=12,height=8,
           device = cairo_pdf)
    
    mg <- g1[[1]] +theme(axis.title.x = element_blank())+coord_cartesian(ylim = c(-4.7,3))+
      g1[[2]] +theme(axis.title = element_blank())+coord_cartesian(ylim = c(-5.6,3))+
      g1[[3]] +theme(axis.title = element_blank())+coord_cartesian(xlim = c(NA,11))+
      g1[[4]] +theme(axis.title = element_blank())+coord_cartesian(xlim = c(NA,15.5))+
      g1[[5]] +theme(axis.title.y = element_blank())+coord_cartesian(ylim = c(-5.3,NA),xlim = c(NA,13))+
      g1[[6]] +theme(axis.title = element_blank())+coord_cartesian(ylim = c(-4.8,NA),xlim = c(NA,13.5))+
      plot_layout(ncol=3)
    ggsave(sprintf("%s/Z_land_vs_even_change_nolab_%s_%s.pdf",save.dir,fb,ra),mg,width=10,height=7,
           device = cairo_pdf)
  }
  
  