library(ggplot2)
library(ggrepel)
library(vegan)
library(ggtext)
library(rELA)
library(doParallel)
#########################################################################
save.dir <- "02_06_05_graphics_summarize_ELA_bootstrap"
dir.create(save.dir)





#dir_02_06 <- "/Users/ngch/Desktop/Now_using/data/sugaFB_ELA/revise_251212/output_spacom/02_06_ELA"
#dir_02_06_03 <- "/Users/ngch/Desktop/Now_using/data/sugaFB_ELA/revise_251212/output_spacom/02_06_03_summarize_ELA_bootstrap"
########################################################################
#read original functions
source("packages/01_1_function.R")
sscolor <- readRDS(sprintf("%s/states_colvector.rds",dir_03_10_02))

fb <- "Fungi"
nsp <- 34
a <- readRDS(sprintf("%s/basin_summary_bootstrap_%s.rds",dir_02_06_03,fb))
true_basin <- readRDS(sprintf("%s/%s/detected_SS_%s.rds",dir_02_06,fb,fb))
gr_para <- readRDS(sprintf("%s/%s/graph_param_SStab_full_%s_Family.rds",dir_02_06,fb,fb))

cmat <- NULL
for(r in 1:nrow(a)){
  cmat <- rbind(cmat,
                cbind(Plant=a$pnam[r],
                      Freq=a$Freq[r],
                      B_name=ifelse(as.character(a$Var1[r]) %in% gr_para[,1],
                                    gr_para[as.character(a$Var1[r]),"rename_SS"],NA),
                      B_color=ifelse(as.character(a$Var1[r]) %in% gr_para[,1],
                                     gr_para[as.character(a$Var1[r]),"color"],NA),
                      B_Full=a$Var1[r]%in%true_basin[true_basin[,"plant"] == a$pnam[r],"SS"],
                      as.data.frame(
                                    matrix(as.numeric(id2bin(as.character(a$Var1[r]),nsp)),
                                           nrow=1))))
  
  
}

#check
tab_check <- c()
true_numbers <- c()
for(p in 1:length(unique(cmat[,"Plant"]))){
  pl <- unique(cmat[,"Plant"])[p]
  true_numbers[p] <- sum(true_basin[,"plant"]==pl)
  tab_check[p] <- sum(cmat[cmat[,"Plant"]==pl,"B_Full"])
}
all(tab_check == true_numbers)

d <- vegdist(cmat[,-c(1:5)],method="jaccard")

pcoa_d <- cmdscale(d, k=2)

df_pcoa <- cbind(cmat[,1:5],as.data.frame(pcoa_d))
colnames(df_pcoa)[6:7] <- c("PCoA1","PCoA2")


g <- ggplot(df_pcoa, aes(x=PCoA1, y=PCoA2)) +
  geom_point(aes(alpha=as.numeric(Freq)),
             color="black",size=2) +
  geom_point(data=function(x){x[x$B_Full,]},
             aes(color=B_name),stroke=1,
             fill="transparent",
             shape=21,size=2) +
  geom_text_repel(data=function(x){x[x$B_Full,]},aes(label=B_name), size=5,
                  alpha=1,
                  box.padding = 0.35,
                  point.padding = 0.5,
                  segment.color = 'grey50') +
  theme_bw() +
  labs(x="PCoA 1",
       y="PCoA 2",
       alpha="Frequency",
       color="Basin detected\nin full ELA") +
  theme(aspect.ratio=1,
        strip.text = element_markdown(size = 14),
        axis.text=element_text(size=12),
        axis.title = element_text(size=14))+
  facet_wrap(~sprintf("*%s*",Plant))+
  scale_color_manual(values=sscolor)

g
ggsave(sprintf("%s/basin_distribution_among_bootstrap_%s.pdf",save.dir,fb),
       plot=g, width=10, height=8)


fb <- "Prokaryota"
nsp <- 45
a <- readRDS(sprintf("%s/basin_summary_bootstrap_%s.rds",dir_02_06_03,fb))
true_basin <- readRDS(sprintf("%s/%s/detected_SS_%s.rds",dir_02_06,fb,fb))
gr_para <- readRDS(sprintf("%s/%s/graph_param_SStab_full_%s_Family.rds",dir_02_06,fb,fb))

cmat <- NULL
for(r in 1:nrow(a)){
  cmat <- rbind(cmat,
                cbind(Plant=a$pnam[r],
                      Freq=a$Freq[r],
                      B_name=ifelse(as.character(a$Var1[r]) %in% gr_para[,1],
                                    gr_para[as.character(a$Var1[r]),"rename_SS"],NA),
                      B_color=ifelse(as.character(a$Var1[r]) %in% gr_para[,1],
                                     gr_para[as.character(a$Var1[r]),"color"],NA),
                      B_Full=a$Var1[r]%in%true_basin[true_basin[,"plant"] == a$pnam[r],"SS"],
                      as.data.frame(
                        matrix(as.numeric(id2bin(as.character(a$Var1[r]),nsp)),
                               nrow=1))))
  
  
}

#check
tab_check <- c()
true_numbers <- c()
for(p in 1:length(unique(cmat[,"Plant"]))){
  pl <- unique(cmat[,"Plant"])[p]
  true_numbers[p] <- sum(true_basin[,"plant"]==pl)
  tab_check[p] <- sum(cmat[cmat[,"Plant"]==pl,"B_Full"])
}
all(tab_check == true_numbers)

d <- vegdist(cmat[,-c(1:5)],method="jaccard")

pcoa_d <- cmdscale(d, k=2)

df_pcoa <- cbind(cmat[,1:5],as.data.frame(pcoa_d))
colnames(df_pcoa)[6:7] <- c("PCoA1","PCoA2")


g <- ggplot(df_pcoa, aes(x=PCoA1, y=PCoA2)) +
  geom_point(aes(alpha=as.numeric(Freq)),
             color="black",size=2) +
  geom_point(data=function(x){x[x$B_Full,]},
             aes(color=B_name),stroke=1,
             fill="transparent",
             shape=21,size=2) +
  geom_text_repel(data=function(x){x[x$B_Full,]},aes(label=B_name), size=5,
                  alpha=1,
                  box.padding = 0.35,
                  point.padding = 0.5,
                  segment.color = 'grey50') +
  theme_bw() +
  labs(x="PCoA 1",
       y="PCoA 2",
       alpha="Frequency",
       color="Basin detected\nin full ELA") +
  theme(aspect.ratio=1,
        strip.text = element_markdown(size = 14),
        axis.text=element_text(size=12),
        axis.title = element_text(size=14))+
  facet_wrap(~sprintf("*%s*",Plant))+
  scale_color_manual(values=sscolor)


g
ggsave(sprintf("%s/basin_distribution_among_bootstrap_%s.pdf",save.dir,fb),
       plot=g, width=10, height=8)
