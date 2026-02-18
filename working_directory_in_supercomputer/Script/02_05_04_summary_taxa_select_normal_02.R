set.seed(1234)
library(rELA)
library(vegan)
library(foreach)
library(doParallel)
library(tidyr)
library(ggplot2)
library(ggstar)

source("packages/01_1_function.R")
#########################################################################
save.dir <- "02_05_04_summary_taxa_select_normal_02"
dir.create(save.dir)







#####################
#parameters
taxa <- "Family"

##R2 of each taxa
#



ltau <- list.files(dir_02_04_04,"mantel_tau_Fungi",full.names = TRUE,recursive = TRUE)

df <- cbind(fb="Fungi",as.data.frame(do.call(rbind,lapply(ltau,readRDS))))

df$best <- FALSE

df[order(df$mean_tau,-df$nSp,-df$rel_th,decreasing = TRUE)[1],"best"] <- TRUE

ltau2 <- list.files(dir_02_04_04,"mantel_tau_Prokaryote",full.names = TRUE,recursive = TRUE)

df2 <- cbind(fb="Prokaryote",as.data.frame(do.call(rbind,lapply(ltau2,readRDS))))

df2$best <- FALSE

df2[order(df2$mean_tau,-df2$nSp,-df2$rel_th,decreasing = TRUE)[1],"best"] <- TRUE

df3 <- rbind(df,df2) 

saveRDS(df3,sprintf("%s/df_heatmap_taxa_select_correlation.rds",save.dir))

g <- ggplot(df,aes(x=nSp,y=rel_th))+
  geom_tile(aes(fill=mean_tau))+
  geom_star(data=function(x){x[x$best,]},
            starshape=1,fill="darkorange",size=3)+
  labs(x="No. of family",y="Binarization thresholds")+
  theme_bw()+
  theme(aspect.ratio = 1,
        strip.text = element_text(size=15),
        axis.title = element_text(size=15),
        axis.text = element_text(size=13))+
  scale_fill_viridis_c()

ggsave(sprintf("%s/taxa_select_correlation_fungi.pdf",save.dir),g,h=6,w=6)

g <- ggplot(df2,aes(x=nSp,y=rel_th))+
  geom_tile(aes(fill=mean_tau))+
  geom_star(data=function(x){x[x$best,]},
            starshape=1,fill="darkorange",size=3)+
  labs(x="No. of family",y="Binarization thresholds")+
  theme_bw()+
  theme(aspect.ratio = 1,
        strip.text = element_text(size=15),
        axis.title = element_text(size=15),
        axis.text = element_text(size=13))+
  scale_fill_viridis_c()

ggsave(sprintf("%s/taxa_select_correlation_bac.pdf",save.dir),g,h=6,w=6)

g <- ggplot(df3,aes(x=nSp,y=rel_th))+
  geom_tile(aes(fill=mean_tau))+
  geom_star(data=function(x){x[x$best,]},
            starshape=1,fill="darkorange",size=3)+
  labs(x="Number of families",y="Binarization thresholds (relative abundance)",
       fill="Mean Kendall's rank correlation coefficients\n(Abundance-based community dissimilarity vs membership dissimilarity in selected families)")+
  theme_bw()+
  theme(aspect.ratio = 1,
        legend.position="bottom",
        legend.title=element_text(size=10),
        legend.text=element_text(size=9),
        strip.text = element_text(size=15),
        axis.title = element_text(size=15),
        axis.text = element_text(size=13))+
  scale_fill_viridis_c()+
  facet_wrap(~fb)

ggsave(sprintf("%s/taxa_select_correlation_merge.pdf",save.dir),g,h=6,w=12)

###################
fb <- "Fungi"

ocmat <- list(Fungi=readRDS(sprintf("%s/ocmat_remove_M2SD_Fungi_Family_relth%s.rds",ELA_prep_dir,df[df$best,"rel_th"])),
              Prokaryote=readRDS(sprintf("%s/ocmat_remove_M2SD_Prokaryota_Family_relth%s.rds",ELA_prep_dir,df2[df2$best,"rel_th"])))

r_res <- readRDS(sprintf("%s/Permanova_sig_taxa_R2_relth%s.rds",dir_02_03_04,df[df$best,"rel_th"]))
r_res1 <- r_res[r_res$fb==fb,]
fdf <- r_res1[order(r_res1[,"R2"],decreasing=TRUE),]

om <- ocmat[[fb]][,fdf[1:df[df$best,"nSp"],"taxa"]]

saveRDS(om,file = paste0(save.dir,"/ELA_input_ocmat_Fungi.rds"))
write.csv(cbind(Family=colnames(om),Occurence=colSums(om),fdf[1:df[df$best,"nSp"],]),file = paste0(save.dir,"/Selected_Family_Fungi.csv"))

fb <- "Prokaryote"

r_res <- readRDS(sprintf("%s/Permanova_sig_taxa_R2_relth%s.rds",dir_02_03_04,df2[df2$best,"rel_th"]))
r_res1 <- r_res[r_res$fb==fb,]
fdf <- r_res1[order(r_res1[,"R2"],decreasing=TRUE),]

om <- ocmat[[fb]][,fdf[1:df2[df2$best,"nSp"],"taxa"]]

saveRDS(om,file = paste0(save.dir,"/ELA_input_ocmat_Prokaryote.rds"))
write.csv(cbind(Family=colnames(om),Occurence=colSums(om),fdf[1:df2[df2$best,"nSp"],]),file = paste0(save.dir,"/Selected_Family_Prok.csv"))