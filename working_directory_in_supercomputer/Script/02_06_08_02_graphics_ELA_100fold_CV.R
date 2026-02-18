
##########################################################################
set.seed(1234)

library(ggplot2)
library(ggtext)

#install.packages("/Users/ngch/Desktop/Now_using/data/sugaFB_ELA/packages/rELA.v0.51.tar.gz")





#########################################################################
save.dir <- "02_06_08_02_graphics_ELA_100fold_CV"
dir.create(save.dir)

fb <- "Fungi"
ls0 <- list.files(dir_02_06_08,"CV_results",recursive = TRUE,full.names = TRUE)

ls <- do.call(rbind,lapply(ls0[grep(fb,ls0)],readRDS))

ls$basin_match <- ls$Basin_CV==ls$Basin_full

res <- NULL
for(i in 1:length(unique(ls$plant))){
  pi <- unique(ls$plant)[i]
  lsi <- sum(ls[ls$plant==pi,"basin_match"])/sum(ls$plant==pi)
  res <- rbind(res,data.frame(plant=sprintf("*%s*",pi),
                       match=c(TRUE,FALSE),
                       count=c(sum(ls[ls$plant==pi,"basin_match"]),
                                    sum(!ls[ls$plant==pi,"basin_match"])),
                       ratio=paste(sum(ls[ls$plant==pi,"basin_match"]),sum(ls$plant==pi),sep="/")
                       ))
}

g <- ggplot(res,aes(x=plant,y=count))+
  geom_bar(aes(fill=match),stat="identity",position="fill")+
  geom_text(data=res[res$match==TRUE,],y=0.95,aes(x=plant,label=ratio))+
  labs(x="Host plant background",y="Proportion of samples",
       fill="same basin\nin full & CV ELA")+
  theme_bw()+
  theme(aspect.ratio=0.6,
        axis.text.x = element_markdown(size=10),
        axis.title = element_text(size=12),
        axis.text.y=element_text(size=10))
g

ggsave(sprintf("%s/ELA_100fold_CV_basinmatch_barplot_%s.pdf",save.dir,fb),
       plot=g,width=6,height=4)

g <- ggplot(ls,aes(x=sprintf("*%s*",plant),y=dist_jaccard))+
  geom_violin()+
  labs(x="Host plant background",y="Jaccard distance between full vs CV ELA")+
  theme_bw()+
  theme(aspect.ratio=0.6,
        axis.text.x = element_markdown(size=10),
        axis.title = element_text(size=12),
        axis.text.y=element_text(size=10))
g

ggsave(sprintf("%s/ELA_100fold_CV_jaccard_violin_%s.pdf",save.dir,fb),
       plot=g,width=6,height=4)

fb <- "Prokaryota"
ls0 <- list.files(dir_02_06_08,"CV_results",recursive = TRUE,full.names = TRUE)

ls <- do.call(rbind,lapply(ls0[grep(fb,ls0)],readRDS))

ls <- rbind(readRDS("/Users/ngch/Desktop/Now_using/data/sugaFB_ELA/revise_251212/output_spacom/CV_results_Fungi_Family_fold1.rds"),
            readRDS("/Users/ngch/Desktop/Now_using/data/sugaFB_ELA/revise_251212/output_spacom/CV_results_Fungi_Family_fold2.rds"))
ls$basin_match <- ls$Basin_CV==ls$Basin_full

res <- NULL
for(i in 1:length(unique(ls$plant))){
  pi <- unique(ls$plant)[i]
  lsi <- sum(ls[ls$plant==pi,"basin_match"])/sum(ls$plant==pi)
  res <- rbind(res,data.frame(plant=sprintf("*%s*",pi),
                              match=c(TRUE,FALSE),
                              count=c(sum(ls[ls$plant==pi,"basin_match"]),
                                      sum(!ls[ls$plant==pi,"basin_match"])),
                              ratio=paste(sum(ls[ls$plant==pi,"basin_match"]),sum(ls$plant==pi),sep="/")
  ))
}


g <- ggplot(res,aes(x=plant,y=count))+
  geom_bar(aes(fill=match),stat="identity",position="fill")+
  geom_text(data=res[res$match==TRUE,],y=0.95,aes(x=plant,label=ratio))+
  labs(x="Host plant background",y="Proportion of samples",
       fill="same basin\nin full & CV ELA")+
  theme_bw()+
  theme(aspect.ratio=0.6,
        axis.text.x = element_markdown(size=8),
        axis.title = element_text(size=12),
        axis.text.y=element_text(size=10))
g

ggsave(sprintf("%s/ELA_100fold_CV_basinmatch_barplot_%s.pdf",save.dir,fb),
       plot=g,width=6,height=4,dpi=300)


g <- ggplot(ls,aes(x=sprintf("*%s*",plant),y=dist_jaccard))+
  geom_violin()+
  labs(x="Host plant background",y="Jaccard distance between full vs CV ELA")+
  theme_bw()+
  theme(aspect.ratio=0.6,
        axis.text.x = element_markdown(size=10),
        axis.title = element_text(size=12),
        axis.text.y=element_text(size=10))
g

ggsave(sprintf("%s/ELA_100fold_CV_jaccard_violin_%s.pdf",save.dir,fb),
       plot=g,width=6,height=4)
