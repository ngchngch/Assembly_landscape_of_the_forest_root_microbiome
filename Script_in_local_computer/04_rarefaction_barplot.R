  
library(ggplot2)
library(tidyr)

dir_01 <- "Output/01_caverage_rarefaction"
dir_02_01 <- "Output_supercomputer/02_01_ELA_prep_abundance_threshold"
dir_03_07 <- "Output/03_07_graphics_Zconv_landchanges_each_biplot_250312"
#########################################################################
save.dir <- "Output/04_rarecurve_barplot"
dir.create(save.dir)

col_genus <- readRDS(sprintf("%s/taxcol_genus.rds",dir_03_07))

fix_col_gns <- gsub("*","",gsub("Fng_|Prok_","",rownames(col_genus)),fixed = TRUE)
names(fix_col_gns) <- rownames(col_genus)

info <- readRDS(sprintf("%s/comp_sample_info_plant2.rds",dir_02_01))

tx <- list(Fungi=readRDS("Base_data/Fungi/taxa_list_mod.rds"),
           Prokaryote=readRDS("Base_data/Bacteria/240404_Megre/merged/taxonomy_list_dada2.rds"))

df <- list(Fungi=readRDS(sprintf("%s/covrarefy_sqtb_fungi.rds",dir_01))$table,
           Prokaryote=readRDS(sprintf("%s/covrarefy_sqtb_bacteria.rds",dir_01))$table)

info_sl <- info[intersect(rownames(df$Fungi),rownames(df$Prokaryote)),]
s_nam <- info_sl[!is.na(info_sl[,"plant"]),"ID"]
length(s_nam)
#########
#--rarecurve

##randamly select 200 samples
fb <- "Fungi"

pdf(sprintf("%s/rarecurve_%s.pdf",save.dir,fb),h=5,w=5)
rarecurve(df[[fb]][sample(s_nam,200),],
          col=color,xlab="Read count",ylab="Number of OTUs",
          label=FALSE)
dev.off()

fb <- "Prokaryote"

pdf(sprintf("%s/rarecurve_%s.pdf",save.dir,fb),h=5,w=5)
rarecurve(df[[fb]][sample(s_nam,200),],
          col=color,xlab="Read count",ylab="Number of OTUs",
          label=FALSE)
dev.off()

#####################
#--Family
taxa <- "Family"
fb <- "Fungi"

tdf <- Taxa.mat(df[[fb]][s_nam,],tx[[fb]],taxa)

#--set color
self <- c(setdiff(colnames(tdf),"Unidentified")[order(colSums(tdf[,setdiff(colnames(tdf),"Unidentified")]),decreasing=TRUE)][1:20],"Others","Unidentified")

colf <- c("#00005B","#FFFF00","#AB005D","#FF8DD7","#008203","#242424","#00D087","#E425A6","#E8FFFF","#0000E1","#049CC6","#3E5025", 
          "#FAAD3B","#AD5700","#31A800","#4561FE" ,"#F06A00","#18E6FF","#69323D","#00657A","gray40","gray80"
)

names(colf) <- self

#--barplot

dat1 <- as.data.frame(cbind(tdf[,setdiff(self,c("Others","Unidentified"))],
                           Others=rowSums(tdf[,
                                                        which(!colnames(tdf) %in% c(setdiff(self,"Others")))]),
                           Unidentified=tdf[,"Unidentified"]))

grel1 <- gather(cbind(ID=rownames(dat1),
                     as.data.frame(dat1)),
               taxa,RA,-c(1))

taxa <- "Family"
fb <- "Prokaryote"

tdf <- Taxa.mat(df[[fb]][s_nam,],tx[[fb]],taxa)

#--set color
selp <- c(setdiff(colnames(tdf),"Unidentified")[order(colSums(tdf[,setdiff(colnames(tdf),"Unidentified")]),decreasing=TRUE)][1:20],"Others","Unidentified")

colp <- c("#F546FF" ,"#00FB00","#00BBFF","#A5788D","#FDBBB1","#5A7966","#303E50","#D3C700","#975B6E","#C2FBAD","#004F5D","#FFBDFF",
          "#FC2200","#FF788C","#5C2C00","#A40600","#99754E","#8083BB","#BD9512","#38340A","gray50","gray70" )


names(colp) <- selp

#--barplot

dat2 <- as.data.frame(cbind(tdf[,setdiff(selp,c("Others","Unidentified"))],
                            Others=rowSums(tdf[,
                                               which(!colnames(tdf) %in% c(setdiff(selp,"Others")))]),
                            Unidentified=tdf[,"Unidentified"]))

grel2 <- gather(cbind(ID=rownames(dat2),
                      as.data.frame(dat2)),
                taxa,RA,-c(1))

##################
#--Genus
#####################
gcolor <- setdiff(color,c(colf,colp,colf_g,colp_g,col_genus[,2]))

taxa <- "Genus"
fb <- "Fungi"

tdf <- Taxa.mat(df[[fb]][s_nam,],tx[[fb]],taxa)

#--set color
self_g <- c(setdiff(colnames(tdf),"Unidentified")[order(colSums(tdf[,setdiff(colnames(tdf),"Unidentified")]),decreasing=TRUE)][1:20],"Others","Unidentified")


colf_g <- ifelse(self_g %in% setdiff(fix_col_gns,"Others"),
                 col_genus[,2][match(self_g,fix_col_gns)],
                 c("#2A1F44","#CF004B","#08B3A2","#640058","#008F62","#FE9200","#C1FFFE","#AAE557","#01D0F5","#004700",
                   "#0044A6","#AA2CF6","#D6E4D2","#696500","#694849","#CB00A3","#97D15E","#450057","#62F2D5","#FFFF93",
                   "gray40","gray80"))

names(colf_g) <- ifelse(self_g %in% c("Others","Unidentified"),
                        self_g,sprintf("*%s*",self_g))

#--barplot

dat3 <- as.data.frame(cbind(tdf[,setdiff(self_g,c("Others","Unidentified"))],
                            Others=rowSums(tdf[,
                                               which(!colnames(tdf) %in% c(setdiff(self_g,"Others")))]),
                            Unidentified=tdf[,"Unidentified"]))

grel3 <- gather(cbind(ID=rownames(dat3),
                      as.data.frame(dat3)),
                taxa,RA,-c(1))


fb <- "Prokaryote"

tdf <- Taxa.mat(df[[fb]][s_nam,],tx[[fb]],taxa)

#--set color
selp_g <- c(setdiff(colnames(tdf),"Unidentified")[order(colSums(tdf[,setdiff(colnames(tdf),"Unidentified")]),decreasing=TRUE)][1:20],"Others","Unidentified")


colp_g <- ifelse(selp_g %in% setdiff(fix_col_gns,"Others"),
                 col_genus[,2][match(selp_g,fix_col_gns)],
                 c("#002100","#E719FF","#876436","#3D200E","#F7E8B2","#003F2E","#62A5B5","#00E300","#2F88FF","#E6D0BE",
                   "#340033","#9DEBFF","#1B988D","#FFC55E","#4201F8","#BFF825","#F69D68","#8C47C5","#E1556B","#C57300",
                   "gray40","gray80"))

names(colp_g) <- ifelse(selp_g %in% c("Others","Unidentified"),
                        selp_g,sprintf("*%s*",selp_g))

#--barplot

dat4 <- as.data.frame(cbind(tdf[,setdiff(selp_g,c("Others","Unidentified"))],
                            Others=rowSums(tdf[,
                                               which(!colnames(tdf) %in% c(setdiff(selp_g,"Others")))]),
                            Unidentified=tdf[,"Unidentified"]))

grel4 <- gather(cbind(ID=rownames(dat4),
                      as.data.frame(dat4)),
                taxa,RA,-c(1))

#--community dissimilarity for ordering samples in barplot

hc <- cbind(dat1/rowSums(dat1),
      2*dat2/rowSums(dat2),
      2*dat3/rowSums(dat3),
      2*dat4/rowSums(dat4)) %>%
  vegdist(method="bray") %>%
  hclust

ord <- s_nam[hc$order]

#--draw 

grel1$ID2 <- factor(grel1$ID,levels=ord)

grel1$taxa2 <- factor(grel1$taxa,levels=self)

g_f_f <- ggplot(grel1,
            aes(x=ID2,y=RA))+
  geom_bar(aes(fill=taxa2,color=taxa2),stat="identity",position = "fill")+
  labs(y="",fill="Family ( Fungi )",color="Family ( Fungi )")+
  theme(aspect.ratio =0.5,
        legend.text = element_markdown(size=8),
        legend.key.size = unit(0.2,"cm"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=8),
        axis.title.y = element_text(size=15),
        axis.title.x = element_blank())+
  scale_fill_manual(values=colf)+
  scale_color_manual(values=colf)

g_f_f

grel2$ID2 <- factor(grel2$ID,levels=ord)

grel2$taxa2 <- factor(grel2$taxa,levels=selp)

g_p_f <- ggplot(grel2,
                aes(x=ID2,y=RA))+
  geom_bar(aes(fill=taxa2,color=taxa2),stat="identity",position = "fill")+
  labs(y="",fill="Family ( Prokaryotes )",color="Family ( Prokaryotes )")+
  theme(aspect.ratio =0.5,
        axis.text.x = element_blank(),
        legend.key.size = unit(0.2,"cm"),
        legend.text = element_markdown(size=8),
        axis.text.y = element_text(size=8),
        axis.title.y = element_text(size=15),
        axis.title.x = element_blank())+
  scale_fill_manual(values=colp)+
  scale_color_manual(values=colp)

g_p_f

grel3$ID2 <- factor(grel3$ID,levels=ord)

grel3$taxa2 <- factor(ifelse(grel3$taxa %in% c("Others","Unidentified"),
                             grel3$taxa,sprintf("*%s*",grel3$taxa)),
                             levels=names(colf_g))

g_f_g <- ggplot(grel3,
                aes(x=ID2,y=RA))+
  geom_bar(aes(fill=taxa2,color=taxa2),stat="identity",position = "fill")+
  labs(y="Relative read count",fill="Genus ( Fungi )",color="Genus ( Fungi )")+
  theme(aspect.ratio =0.5,
        legend.text = element_markdown(size=8),
        axis.text.x = element_blank(),
        legend.key.size = unit(0.2,"cm"),
        axis.text.y = element_text(size=8),
        axis.title.y = element_text(size=15),
        axis.title.x = element_blank())+
  scale_fill_manual(values=colf_g)+
  scale_color_manual(values=colf_g)

g_f_g


grel4$ID2 <- factor(grel4$ID,levels=ord)

grel4$taxa2 <- factor(ifelse(grel4$taxa %in% c("Others","Unidentified"),
                             grel4$taxa,sprintf("*%s*",grel4$taxa)),
                      levels=names(colp_g))

g_p_g <- ggplot(grel4,
                aes(x=ID2,y=RA))+
  geom_bar(aes(fill=taxa2,color=taxa2),stat="identity",position = "fill")+
  labs(x="Root-tips",y="",fill="Genus ( Prokaryotes )",color="Genus ( Prokaryotes )")+
  theme(aspect.ratio =0.5,
        legend.text = element_markdown(size=8),
        legend.key.size = unit(0.2,"cm"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=8),
        axis.title.y = element_text(size=15))+
  scale_fill_manual(values=colp_g)+
  scale_color_manual(values=colp_g)

g_p_g

####

mg <- g_f_f+g_f_g+g_p_f+g_p_g+plot_layout(ncol=1,
                                    axis_titles = "collect")

ggsave(sprintf("%s/FigS_barplot_merge.pdf",save.dir),mg,h=8,w=8)

###
#mean abundanse
colMeans(dat1/rowSums(dat1))

colMeans(dat3/rowSums(dat3))

colMeans(dat2/rowSums(dat2))

colMeans(dat4/rowSums(dat4))

####
#--sample number each host

gt <- as.data.frame(table(info_sl[s_nam,"plant"]))
g <- ggplot(gt[order(gt$Freq,decreasing = TRUE),],
            aes(x=reorder(sprintf("*%s*",Var1),-Freq),
                y=Freq))+
  geom_bar(stat="identity",aes(fill=Var1))+
  labs(x="Host plant",y="No. of samples")+
  theme_bw()+
  theme(axis.text.x = element_markdown(size=12,angle=45,hjust=1),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=15),
        legend.position = "none",
        aspect.ratio = 1.2)

g
ggsave(sprintf("%s/FigS3_barplot_sample_number.pdf",save.dir),g,h=8,w=8)
