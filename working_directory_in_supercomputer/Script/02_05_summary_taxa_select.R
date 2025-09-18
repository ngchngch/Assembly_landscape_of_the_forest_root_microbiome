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
save.dir <- "02_05_summary_taxa_select"
dir.create(save.dir)







#####################
#parameters
taxa <- "Family"

##R2 of each taxa
r_res <- readRDS(sprintf("%s/Permanova_sig_taxa_R2.rds",dir_02_03))

ocmat <- list(Fungi=readRDS(sprintf("%s/ocmat_remove_M2SD_Fungi_Family.rds",ELA_prep_dir)),
              Prokaryote=readRDS(sprintf("%s/ocmat_remove_M2SD_Prokaryota_Family.rds",ELA_prep_dir)))


ltau <- list.files(dir_02_04,"mantel_tau_Fungi",full.names = TRUE,recursive = TRUE)

df <- cbind(fb="Fungi",as.data.frame(do.call(rbind,lapply(ltau,readRDS))))

df$best <- FALSE

df[order(df$mean_tau,-df$nSp,decreasing = TRUE)[1],"best"] <- TRUE

ltau2 <- list.files(dir_02_04,"mantel_tau_Prokaryote",full.names = TRUE,recursive = TRUE)

df2 <- cbind(fb="Prokaryote",as.data.frame(do.call(rbind,lapply(ltau2,readRDS))))

df2$best <- FALSE

df2[order(df2$mean_tau,-df2$nSp,decreasing = TRUE)[1],"best"] <- TRUE

df3 <- rbind(df,df2) 


g <- ggplot(df3,aes(x=nSp,y=mean_tau))+
  geom_line()+
  geom_point()+
  geom_star(data=function(x){x[x$best,]},
            starshape=1,fill="darkorange",size=4)+
  labs(x="No. of family",y="Mean kendall's correlation coefficient")+
  theme_bw()+
  theme(aspect.ratio = 1,
        strip.text = element_text(size=15),
        axis.title = element_text(size=15),
        axis.text = element_text(size=13))+
  facet_wrap(~fb,nrow=1,scales = "free_y")

ggsave(sprintf("%s/taxa_select_correlation.pdf",save.dir),g,h=6,w=12)

###################
fb <- "Fungi"

r_res1 <- r_res[r_res$fb==fb,]
fdf <- r_res1[order(r_res1[,"R2"],decreasing=TRUE),]

om <- ocmat[[fb]][,fdf[1:df[df$best,"nSp"],"taxa"]]

saveRDS(om,file = paste0(save.dir,"/ELA_input_ocmat_Fungi.rds"))

fb <- "Prokaryote"

r_res1 <- r_res[r_res$fb==fb,]
fdf <- r_res1[order(r_res1[,"R2"],decreasing=TRUE),]

om <- ocmat[[fb]][,fdf[1:df[df$best,"nSp"],"taxa"]]

saveRDS(om,file = paste0(save.dir,"/ELA_input_ocmat_Prokaryote.rds"))
