
##########################################################################
set.seed(1234)

library(ggplot2)
library(ggtext)

#install.packages("/Users/ngch/Desktop/Now_using/data/sugaFB_ELA/packages/rELA.v0.51.tar.gz")





#########################################################################
save.dir <- "06_02_gLV_graphics_samplesize_validation"
dir.create(save.dir)

ls0 <- list.files(dir_06_01,"summary_ELA_result",recursive = TRUE,full.names = TRUE)

seed <- sapply(strsplit(ls0,"/"),"[",2)

Nsp <- sapply(strsplit(ls0,"_"),function(x){
  as.numeric(gsub("Nsp|.rds","",x[length(x)]))
})

connect <- sapply(strsplit(ls0,"_"),function(x){
  as.numeric(gsub("conn","",x[length(x)-1]))
})

sm <- NULL
for(l in 1:length(ls0)){
  sm <- rbind(sm,cbind(seed=seed[l],Nsp=Nsp[l],connect=connect[l],
              readRDS(ls0[l])))
}

saveRDS(sm,sprintf("%s/gLV_sample_size_validation_summary.rds",save.dir))

g <- ggplot(sm,aes(x=as.factor(nsample),y=nbasin))+
  geom_boxplot(outlier.color = NA)+
  geom_jitter(width = 0.2, height = 0,alpha=0.5)+
  facet_grid(seed~Nsp,scales = "free")+
  theme_bw()+
  labs(x="Number of samples used in energy landscape analysis",
       y="Number of basins identified in ELA")+
  theme(aspect.ratio=1,
        axis.text.x = element_text(size=10),
        axis.title = element_text(size=12),
        axis.text.y=element_text(size=10),
        strip.text = element_text(size=12))

ggsave(sprintf("%s/gLV_sample_size_nbasin.pdf",save.dir),
       plot=g,width=16,height=8)
  
g <- ggplot(unique(sm[,c("Nsp","nsample","JSdist_bfreq")]),
            aes(x=nsample,y=JSdist_bfreq,color=seed))+
  geom_line()+
  geom_point()+
  facet_wrap(~Nsp,scales = "free",nrow = 2)+
  theme_bw()+
  labs(x="Number of samples used in energy landscape analysis",
       y="Jansen-Shannon distance of basins propotion")+
  theme(aspect.ratio=1,
        axis.text.x = element_text(size=10),
        axis.title = element_text(size=12),
        axis.text.y=element_text(size=10),
        strip.text = element_text(size=12))

ggsave(sprintf("%s/gLV_sample_size_basin_freq_dist.pdf",save.dir),
       plot=g,width=16,height=8)
