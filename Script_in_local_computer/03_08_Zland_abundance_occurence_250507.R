
###########################################
##########################################################################



source("packages/01_1_function.R")
set.seed(1234)

library(renv)
library(ggstar)
library(ggplot2)
library(ggtext)
library(patchwork)

library(ggrepel)

#detectCores()
n.core <- 8 # not threads!!

dir_03_02_02 <- "Output/03_02_02_Zconv_fixP"
dir_03_07 <- "Output/03_07_graphics_Zconv_landchanges_each_biplot_250312"
ELA_prep_dir <- "Output_supercomputer/02_01_ELA_prep_abundance_threshold"

#########################################################################
save.dir <- "Output/03_08_Zland_abundance_occurence_250623"
dir.create(save.dir)

info <- readRDS(sprintf("%s/comp_sample_info_plant2.rds",ELA_prep_dir))
rmat <- readRDS(sprintf("%s/NoCLR_ExpVar_clrRA_target_taxa.rds",ELA_prep_dir))
ramat <- cbind(rmat$Fungi,rmat$Prokaryota)

ocmat <- list(Fungi=readRDS(sprintf("%s/ocmat_remove_M2SD_Fungi_Family.rds",ELA_prep_dir)),
              Prokaryote=readRDS(sprintf("%s/ocmat_remove_M2SD_Prokaryota_Family.rds",ELA_prep_dir)))

tx_f <- readRDS("../../Base_data/Fungi/taxa_list_mod.rds")
tx_p <- as.data.frame(readRDS("../../Base_data/OTU_basedata_set/taxonomy_list_dada2.rds"))
tx_m <- rbind(tx_f[,c("Order","Genus")],tx_b[,c("Order","Genus")])



zval_f_org <- readRDS(sprintf("%s/zvalue_Fungi.rds",dir_03_02_02))
zval_p_org <- readRDS(sprintf("%s/zvalue_Prokaryote.rds",dir_03_02_02))


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


pland_BH <- p.adjust(c(zval_f_org$p_land,zval_p_org$p_land),method = "BH")
zval_f_org$p_land_BH <- pland_BH[1:length(zval_f_org$p_land)]
zval_p_org$p_land_BH <- pland_BH[(length(zval_f_org$p_land)+1):length(pland_BH)]

peven_BH <- p.adjust(c(zval_f_org$p_even,zval_p_org$p_even),method = "BH")
zval_f_org$p_even_BH <- peven_BH[1:length(zval_f_org$p_even)]
zval_p_org$p_even_BH <- peven_BH[(length(zval_f_org$p_even)+1):length(peven_BH)]

#reset inf values
zval_l <- list(perc25=rbind(cbind(fb="Fungi",zval_f_org[zval_f_org$ra=="perc25",]),
                            cbind(fb="Prokaryote",zval_p_org[zval_p_org$ra=="perc25",])),
               median=rbind(cbind(fb="Fungi",zval_f_org[zval_f_org$ra=="median",]),
                            cbind(fb="Prokaryote",zval_p_org[zval_p_org$ra=="median",])),
               perc75=rbind(cbind(fb="Fungi",zval_f_org[zval_f_org$ra=="perc75",]),
                            cbind(fb="Prokaryote",zval_p_org[zval_p_org$ra=="perc75",])))


#coloring taxa selection
taxcol_genus <- readRDS(sprintf("%s/taxcol_genus.rds",dir_03_07))


for(i in 1:length(zval_l)){#i <- 2
  ra <- names(zval_l)[i]
  le <- zval_l[[i]]
  lef <- le[!is.na(le$p_land),]
  lef$z_land[is.na(lef$z_land)] <- 0
  lef$z_even[is.na(lef$z_even)] <- 0
  
  lef$landscape <- ifelse(lef$p_land_BH<0.05,
                          ifelse(lef$p_even_BH<0.025,
                                 ifelse(lef$z_even>0,
                                        "with \u0394*evenness* increase",
                                        "with \u0394*evenness* reduction"),
                                 "without \u0394*evenness* change"),
                          "No re-organization")
  
  
  
  lef$target2 <- ifelse(lef$target %in% rownames(taxcol_genus),
                        lef$target,
                        ifelse(grepl("Fng_",lef$target),"Fng_Others",
                               "Prok_Others"))
  
  lef$plant2 <- factor(sprintf("*%s*",lef$plant),
                       sprintf("*%s*",c("Pinus","Larix","Betula","Populus","Acer","Juglans")))
  
  #top 2 selct for plot annotation
  
  #topn <- 3
  lef$Top <- FALSE
  topmat <- NULL
  topmat_e_upper <- NULL
  topmat_e_lower <- NULL
  
  for(k in 1:length(unique(lef$fb))){#i<- 1
    for(j in 1:length(unique(lef$plant))){#j<-1
      #top2 landscape change
      sel_lef <- lef[which(lef$fb==unique(lef$fb)[k]&lef$plant==unique(lef$plant)[j]&lef$p_land_BH<0.05),]
      
      if(length(sel_lef$z_land)>1){
        th <- quantile(sel_lef$z_land,0.95)
        
        sel_tag <-   sel_lef[sel_lef$z_land>=th,"target"]
        
        topmat <-   rbind(topmat,sel_lef[sel_lef$z_land>=th,])
      }else{
        sel_tag <-   sel_lef[,"target"]
        topmat <-   rbind(topmat,sel_lef)
      }
      lef[which(lef$fb==unique(lef$fb)[k]&lef$plant==unique(lef$plant)[j]&lef$target%in%sel_tag),"Top"] <- TRUE
      
      
    
  }}
  
  
  lef$abundance <- NA
  lef$occurence <- NA
  for(a in 1:length(unique(lef$target))){
    for(b in 1:length(unique(lef$plant))){#a <- 115;b<- 1
      
      rd_sel <- ramat[which(rownames(ramat) %in% info[info$plant==unique(lef$plant)[b],"ID"]),
                      which(colnames(ramat) %in% unique(lef$Taxa)[a])]
      
      lef[which(lef$target==unique(lef$target)[a]&lef$plant==unique(lef$plant)[b]),"abundance"] <- mean(rd_sel)
      lef[which(lef$target==unique(lef$target)[a]&lef$plant==unique(lef$plant)[b]),"occurence"] <- sum(rd_sel>0)/length(rd_sel)
    }
  }
  
  
  lef$target_label <- ifelse(lef$target %in% tx_f[,"Genus"],
                             sprintf("~italic('%s')",lef$Taxa),
                             sprintf("~italic('%s')",lef$Taxa))
  
 write.csv(lef,sprintf("%s/landscape_%s.csv",save.dir,ra))
  #####################################
########################################
  
  plants <- c("Acer","Betula","Pinus","Populus","Larix","Juglans")
  
  
  fb <- "Fungi"
  
  glef <- na.omit(lef[which(lef$fb==fb&lef$plant%in%plants),])
  
  cr_ab <- c()
  cr_ab2 <- c()
  cr_oc <- c()
  p_ab <- c()
  p_ab2 <- c()
  p_oc <- c()
  
 for(pl in 1:length(plants)){#pl <- 1
   crab <-  cor.test(glef[which(glef$plant==plants[pl]&glef$target_k=="Fungi"),
                          "abundance"],
                     glef[which(glef$plant==plants[pl]&glef$target_k=="Fungi"),
                          "z_land"],method="kendall")
   crab2 <-  cor.test(glef[which(glef$plant==plants[pl]&glef$target_k=="Prokaryote"),
                          "abundance"],
                     glef[which(glef$plant==plants[pl]&glef$target_k=="Prokaryote"),
                          "z_land"],method="kendall")
   
   
  cr_ab[pl] <- crab$estimate
   p_ab[pl] <- crab$p.value
   
   cr_ab2[pl] <- crab2$estimate
   p_ab2[pl] <- crab2$p.value
   
   croc <-  cor.test(glef[which(glef$plant==plants[pl]),"occurence"],
                     glef[which(glef$plant==plants[pl]),"z_land"],method="kendall")
   cr_oc[pl] <- croc$estimate
   p_oc[pl] <- croc$p.value
 }
 
  write.csv(cbind(plant=plants,
                  cr_ab_f=cr_ab,
                  p_ab_f=p.adjust(p_ab,method="BH"),
                  cr_ab_p=cr_ab2,
                  p_ab_p=p.adjust(p_ab2,method="BH"),
                  cr_oc=cr_oc,
                  p_oc=p.adjust(p_oc,method="BH")),
            sprintf("%s/cor_test_%s_%s.csv",save.dir,ra,fb))
  
  g1 <- ggplot(glef,
               aes(x=log(abundance),y=z_land))+
    #geom_vline(xintercept = 0,linewidth=0.2,color="gray60")+
    geom_hline(yintercept = 0,linewidth=0.2,color="gray60")+
    # geom_vline(data = data.frame(plant2=sprintf("*%s*",
    #                                             inf_border_l[which(inf_border_l[,2]==ra),3]),
    #                              val=as.numeric(inf_border_l[which(inf_border_l[,2]==ra),4])),
    #            aes(xintercept = as.numeric(val)),
    #            linetype="dashed",color="gray40")+
    
    geom_star(#data=function(x){x[x$landscape=="No re-organization",]},
              aes(size=landscape,fill=target_p,starshape=landscape))+
    # geom_star(data=function(x){x[x$landscape!="No re-organization",]},
    #           aes(size=landscape,fill=tag_ord,starshape=landscape))+
    labs(y="Z-standardized \u0394*topography*\n( Fungal community )",
         x="log (mean relative read count)",
         fill="Genus",
         size="Energy landscape re-organization",
         starshape="Energy landscape re-organization")+
    theme_bw()+
    scale_fill_manual(values = col_phyl)+
    scale_size_manual(values = c("No re-organization"=1,
                                 "without \u0394*evenness* change"=3,
                                 "with \u0394*evenness* reduction"=3,
                                 "with \u0394*evenness* increase"=3))+
    scale_starshape_manual(values = c("No re-organization"=15,
                                      "without \u0394*evenness* change"=28,
                                      "with \u0394*evenness* reduction"=23,
                                      "with \u0394*evenness* increase"=11))+
    theme(aspect.ratio = 1,
          legend.position = "none",
          strip.text = element_markdown(size = 12),
          axis.text = element_text(size = 10),
          axis.title = element_markdown(size = 12),
          legend.text = element_markdown(size=9))+
    guides(fill=guide_legend(override.aes = list(starshape=15,size=2),
                             ncol=2))+
    facet_grid(plant2~target_k,scales = "free")
  
  g1
  
  ggsave(sprintf("%s/Z_land_vs_RA_nolab_All_%s_%s.pdf",save.dir,fb,ra),g1,h=14,w=5,
         device = cairo_pdf)
  
  ggsave(sprintf("%s/Z_land_vs_RA_All_%s_%s.pdf",save.dir,fb,ra),g1+
           geom_text_repel(data=function(x){x[x$Top,]},
                           min.segment.length = 0.1,
                           aes(label=target_label),size=3,
                           parse = TRUE,
                           seed = 715240),h=14,w=5,
         device = cairo_pdf)
  
  ggsave(plot=g_legend(g1+theme(legend.position="right")),
         sprintf("%s/legend_Z_land_vs_RA_All_%s_%s.pdf",save.dir,fb,ra),width=6,height=9,
         device = cairo_pdf)
  
  g1 <- ggplot(glef,
               aes(x=occurence,y=z_land))+
    #geom_vline(xintercept = 0,linewidth=0.2,color="gray60")+
    geom_hline(yintercept = 0,linewidth=0.2,color="gray60")+
    # geom_vline(data = data.frame(plant2=sprintf("*%s*",
    #                                             inf_border_l[which(inf_border_l[,2]==ra),3]),
    #                              val=as.numeric(inf_border_l[which(inf_border_l[,2]==ra),4])),
    #            aes(xintercept = as.numeric(val)),
    #            linetype="dashed",color="gray40")+
    
    geom_star(#data=function(x){x[x$landscape=="No re-organization",]},
      aes(size=landscape,fill=target_p,starshape=landscape))+
    # geom_star(data=function(x){x[x$landscape!="No re-organization",]},
    #           aes(size=landscape,fill=tag_ord,starshape=landscape))+
    labs(y="Z-standardized \u0394*topography*\n( Fungal community )",
         x="Occupancy (%)",
         fill="Genus",
         size="Energy landscape re-organization",
         starshape="Energy landscape re-organization")+
    theme_bw()+
    scale_fill_manual(values = col_phyl)+
    scale_size_manual(values = c("No re-organization"=1,
                                 "without \u0394*evenness* change"=3,
                                 "with \u0394*evenness* reduction"=3,
                                 "with \u0394*evenness* increase"=3))+
    scale_starshape_manual(values = c("No re-organization"=15,
                                      "without \u0394*evenness* change"=28,
                                      "with \u0394*evenness* reduction"=23,
                                      "with \u0394*evenness* increase"=11))+
    theme(aspect.ratio = 1,
          legend.position = "none",
          strip.text = element_markdown(size = 12),
          axis.text = element_text(size = 10),
          axis.title = element_markdown(size = 12),
          legend.text = element_markdown(size=9))+
    guides(fill=guide_legend(override.aes = list(starshape=15,size=2),
                             ncol=2))+
    facet_wrap(~plant2,scales = "free")
  
  g1
  
  ggsave(sprintf("%s/Z_land_vs_Occ_nolab_All_%s_%s.pdf",save.dir,fb,ra),g1,
         h=6,w=8,
         device = cairo_pdf)
  
  ggsave(sprintf("%s/Z_land_vs_Occ_All_%s_%s.pdf",save.dir,fb,ra),g1+
           geom_text_repel(data=function(x){x[x$Top,]},
                           min.segment.length = 0.1,
                           aes(label=target_label),size=3,
                           parse = TRUE,
                           seed = 715240),width=8,height=6,
         device = cairo_pdf)
  
  ggsave(plot=g_legend(g1+theme(legend.position="right")),
         sprintf("%s/legend_Z_land_vs_Occ_All_%s_%s.pdf",save.dir,fb,ra),width=6,height=9,
         device = cairo_pdf)
  
  
  fb <- "Prokaryote"
  glef <- na.omit(lef[which(lef$fb==fb&lef$plant%in%plants),])
  cr_ab <- c()
  cr_ab2 <- c()
  cr_oc <- c()
  p_ab <- c()
  p_ab2 <- c()
  p_oc <- c()
  
  for(pl in 1:length(plants)){#pl <- 1
    crab <-  cor.test(glef[which(glef$plant==plants[pl]&glef$target_k=="Fungi"),
                           "abundance"],
                      glef[which(glef$plant==plants[pl]&glef$target_k=="Fungi"),
                           "z_land"],method="kendall")
    crab2 <-  cor.test(glef[which(glef$plant==plants[pl]&glef$target_k=="Prokaryote"),
                            "abundance"],
                       glef[which(glef$plant==plants[pl]&glef$target_k=="Prokaryote"),
                            "z_land"],method="kendall")
    
    
    cr_ab[pl] <- crab$estimate
    p_ab[pl] <- crab$p.value
    
    cr_ab2[pl] <- crab2$estimate
    p_ab2[pl] <- crab2$p.value
    
    croc <-  cor.test(glef[which(glef$plant==plants[pl]),"occurence"],
                      glef[which(glef$plant==plants[pl]),"z_land"],method="kendall")
    cr_oc[pl] <- croc$estimate
    p_oc[pl] <- croc$p.value
  }
  
  write.csv(cbind(plant=plants,
                  cr_ab_f=cr_ab,
                  p_ab_f=p.adjust(p_ab,method="BH"),
                  cr_ab_p=cr_ab2,
                  p_ab_p=p.adjust(p_ab2,method="BH"),
                  cr_oc=cr_oc,
                  p_oc=p.adjust(p_oc,method="BH")),
            sprintf("%s/cor_test_%s_%s.csv",save.dir,ra,fb))
  
  g1 <- ggplot(glef,
               aes(x=log(abundance),y=z_land))+
    #geom_vline(xintercept = 0,linewidth=0.2,color="gray60")+
    geom_hline(yintercept = 0,linewidth=0.2,color="gray60")+
    # geom_vline(data = data.frame(plant2=sprintf("*%s*",
    #                                             inf_border_l[which(inf_border_l[,2]==ra),3]),
    #                              val=as.numeric(inf_border_l[which(inf_border_l[,2]==ra),4])),
    #            aes(xintercept = as.numeric(val)),
    #            linetype="dashed",color="gray40")+
    
    geom_star(#data=function(x){x[x$landscape=="No re-organization",]},
      aes(size=landscape,fill=target_p,starshape=landscape))+
    # geom_star(data=function(x){x[x$landscape!="No re-organization",]},
    #           aes(size=landscape,fill=tag_ord,starshape=landscape))+
    labs(y="Z-standardized \u0394*topography*\n( Prokaryotic community )",
         x="log (mean relative read count)",
         fill="Genus",
         size="Energy landscape re-organization",
         starshape="Energy landscape re-organization")+
    theme_bw()+
    scale_fill_manual(values = col_phyl)+
    scale_size_manual(values = c("No re-organization"=1,
                                 "without \u0394*evenness* change"=3,
                                 "with \u0394*evenness* reduction"=3,
                                 "with \u0394*evenness* increase"=3))+
    scale_starshape_manual(values = c("No re-organization"=15,
                                      "without \u0394*evenness* change"=28,
                                      "with \u0394*evenness* reduction"=23,
                                      "with \u0394*evenness* increase"=11))+
    theme(aspect.ratio = 1,
          legend.position = "none",
          strip.text = element_markdown(size = 12),
          axis.text = element_text(size = 10),
          axis.title = element_markdown(size = 12),
          legend.text = element_markdown(size=9))+
    guides(fill=guide_legend(override.aes = list(starshape=15,size=2),
                             ncol=2))+
    facet_grid(plant2~target_k,scales = "free")
  
  g1
  
  ggsave(sprintf("%s/Z_land_vs_RA_nolab_All_%s_%s.pdf",save.dir,fb,ra),g1,h=14,w=5,
         device = cairo_pdf)
  
  ggsave(sprintf("%s/Z_land_vs_RA_All_%s_%s.pdf",save.dir,fb,ra),g1+
           geom_text_repel(data=function(x){x[x$Top,]},
                           min.segment.length = 0.1,
                           aes(label=target_label),size=3,
                           parse = TRUE,
                           seed = 715240),h=14,w=5,
         device = cairo_pdf)
  
  ggsave(plot=g_legend(g1+theme(legend.position="right")),
         sprintf("%s/legend_Z_land_vs_RA_All_%s_%s.pdf",save.dir,fb,ra),width=6,height=9,
         device = cairo_pdf)
  
  g1 <- ggplot(glef,
               aes(x=occurence,y=z_land))+
    #geom_vline(xintercept = 0,linewidth=0.2,color="gray60")+
    geom_hline(yintercept = 0,linewidth=0.2,color="gray60")+
    # geom_vline(data = data.frame(plant2=sprintf("*%s*",
    #                                             inf_border_l[which(inf_border_l[,2]==ra),3]),
    #                              val=as.numeric(inf_border_l[which(inf_border_l[,2]==ra),4])),
    #            aes(xintercept = as.numeric(val)),
    #            linetype="dashed",color="gray40")+
    
    geom_star(#data=function(x){x[x$landscape=="No re-organization",]},
      aes(size=landscape,fill=target_p,starshape=landscape))+
    # geom_star(data=function(x){x[x$landscape!="No re-organization",]},
    #           aes(size=landscape,fill=tag_ord,starshape=landscape))+
    labs(y="Z-standardized \u0394*topography*\n( Prokaryotic community )",
         x="Occupancy (%)",
         fill="Genus",
         size="Energy landscape re-organization",
         starshape="Energy landscape re-organization")+
    theme_bw()+
    scale_fill_manual(values = col_phyl)+
    scale_size_manual(values = c("No re-organization"=1,
                                 "without \u0394*evenness* change"=3,
                                 "with \u0394*evenness* reduction"=3,
                                 "with \u0394*evenness* increase"=3))+
    scale_starshape_manual(values = c("No re-organization"=15,
                                      "without \u0394*evenness* change"=28,
                                      "with \u0394*evenness* reduction"=23,
                                      "with \u0394*evenness* increase"=11))+
    theme(aspect.ratio = 1,
          legend.position = "none",
          strip.text = element_markdown(size = 12),
          axis.text = element_text(size = 10),
          axis.title = element_markdown(size = 12),
          legend.text = element_markdown(size=9))+
    guides(fill=guide_legend(override.aes = list(starshape=15,size=2),
                             ncol=2))+
    facet_wrap(~plant2,scales = "free")
  
  g1
  
  ggsave(sprintf("%s/Z_land_vs_Occ_nolab_All_%s_%s.pdf",save.dir,fb,ra),g1,
         h=6,w=8,
         device = cairo_pdf)
  
  ggsave(sprintf("%s/Z_land_vs_Occ_All_%s_%s.pdf",save.dir,fb,ra),g1+
           geom_text_repel(data=function(x){x[x$Top,]},
                           min.segment.length = 0.1,
                           aes(label=target_label),size=3,
                           parse = TRUE,
                           seed = 715240),width=8,height=6,
         device = cairo_pdf)
  
  ggsave(plot=g_legend(g1+theme(legend.position="right")),
         sprintf("%s/legend_Z_land_vs_Occ_All_%s_%s.pdf",save.dir,fb,ra),width=6,height=9,
         device = cairo_pdf)
  
  
}

