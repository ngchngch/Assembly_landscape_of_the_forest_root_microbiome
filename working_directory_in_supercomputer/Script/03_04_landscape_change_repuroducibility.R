#setwd("/Users/ngch/Desktop/Now_using/data/sugaFB_ELA/analysis_series_SpCom/02_ELA_240914")  

set.seed(1234)

####################
library(ggplot2)



save.dir <- "03_04_landscape_change_repuroducibility"
dir.create(save.dir)
###################

ls <- list.files(pattern="ELA_SSprob_diff",
                 dir_03_03,
                 full.names = TRUE, recursive = TRUE)

fl <- lapply(ls[grep("Fungi",ls)], readRDS)

fdf <- do.call(rbind,fl)

var_mat <- NULL
for(i in 1:length(unique(fdf$plant))){
  for(j in 1:length(unique(fdf$ra))){#i <- 1;j <- 1;k <- 1
    for(k in 1:length(unique(fdf$Taxa))){
      reps <- fdf[fdf$ra == unique(fdf$ra)[j] & fdf$plant == unique(fdf$plant)[i] & fdf$Taxa == unique(fdf$Taxa)[k],]
      
      var_mat <- rbind(var_mat,
                       data.frame(ra= unique(fdf$ra)[j],
                                  plant= unique(fdf$plant)[i],
                                  Taxa= unique(fdf$Taxa)[k],
                                  var_dland=var(reps$d_land),var_deven=var(reps$d_even),
                                  CV_dland=sd(reps$d_land)/mean(reps$d_land),
                                  var_deven=var(reps$d_even),
                                  CV_deven=sd(reps$d_even)/mean(reps$d_even))
                       )
    }
     }
}


var_mat[var_mat$ra == "perc25","ra"] <- "25%"
var_mat[var_mat$ra == "median","ra"] <- "50%"
var_mat[var_mat$ra == "perc75","ra"] <- "75%"


saveRDS(var_mat,sprintf("%s/var_mat_Fungi.rds",save.dir))

g <- ggplot(var_mat,aes(x=var_dland))+
  geom_histogram()+
  facet_grid(plant~ra)+
  theme_bw()+
  theme(aspect.ratio = 1,
        axis.text.x = element_text(angle=45),
        axis.text = element_text(size=12),
        axis.title = element_text(size=15))

ggsave(sprintf("%s/histogram_var_d_land_rep_F.pdf",save.dir),g,width=10,height=10)

g <- ggplot(var_mat,aes(x=var_deven))+
  geom_histogram()+
  facet_grid(plant~ra)+
  theme_bw()+
  theme(aspect.ratio = 1,
        axis.text.x = element_text(angle=45),
        axis.text = element_text(size=12),
        axis.title = element_text(size=15))

ggsave(sprintf("%s/histogram_var_d_even_rep_F.pdf",save.dir),g,width=10,height=10)


g <- ggplot(var_mat,aes(x=CV_dland))+
  geom_histogram()+
  facet_grid(plant~ra)+
  theme_bw()+
  theme(aspect.ratio = 1,
        axis.text.x = element_text(angle=45),
        axis.text = element_text(size=12),
        axis.title = element_text(size=15))

ggsave(sprintf("%s/histogram_CV_d_land_rep_F.pdf",save.dir),g,width=10,height=10)

g <- ggplot(var_mat,aes(x=CV_deven))+
  geom_histogram()+
  facet_grid(plant~ra)+
  theme_bw()+
  theme(aspect.ratio = 1,
        axis.text.x = element_text(angle=45),
        axis.text = element_text(size=12),
        axis.title = element_text(size=15))

ggsave(sprintf("%s/histogram_CV_d_even_rep_F.pdf",save.dir),g,width=10,height=10)

#########################################

bl <- lapply(ls[grep("Prok",ls)], readRDS)

bdf <- do.call(rbind,bl)


var_mat <- NULL
for(i in 1:length(unique(bdf$plant))){
  for(j in 1:length(unique(bdf$ra))){#i <- 1;j <- 1;k <- 1
    for(k in 1:length(unique(bdf$Taxa))){
      reps <- bdf[bdf$ra == unique(bdf$ra)[j] & bdf$plant == unique(bdf$plant)[i] & bdf$Taxa == unique(bdf$Taxa)[k],]
      
      var_mat <- rbind(var_mat,
                       data.frame(ra= unique(bdf$ra)[j],
                                  plant= unique(bdf$plant)[i],
                                  Taxa= unique(bdf$Taxa)[k],
                                  var_dland=var(reps$d_land),var_deven=var(reps$d_even),
                                  CV_dland=sd(reps$d_land)/mean(reps$d_land),
                                  var_deven=var(reps$d_even),
                                  CV_deven=sd(reps$d_even)/mean(reps$d_even)))
    }
  }
}


var_mat[var_mat$ra == "perc25","ra"] <- "25%"
var_mat[var_mat$ra == "median","ra"] <- "50%"
var_mat[var_mat$ra == "perc75","ra"] <- "75%"

saveRDS(var_mat,sprintf("%s/var_mat_Prokaryote.rds",save.dir))

g <- ggplot(var_mat,aes(x=var_dland))+
  geom_histogram()+
  facet_grid(plant~ra)+
  theme_bw()+
  theme(aspect.ratio = 1,
        axis.text = element_text(size=12),
        axis.text.x = element_text(angle=45),
        axis.title = element_text(size=15))

ggsave(sprintf("%s/histogram_var_d_land_rep_B.pdf",save.dir),g,width=10,height=10)

g <- ggplot(var_mat,aes(x=var_deven))+
  geom_histogram()+
  facet_grid(plant~ra)+
  theme_bw()+
  theme(aspect.ratio = 1,
        axis.text.x = element_text(angle=45),
        axis.text = element_text(size=12),
        axis.title = element_text(size=15))

ggsave(sprintf("%s/histogram_var_d_even_rep_B.pdf",save.dir),g,width=10,height=10)


g <- ggplot(var_mat,aes(x=CV_dland))+
  geom_histogram()+
  facet_grid(plant~ra)+
  theme_bw()+
  theme(aspect.ratio = 1,
        axis.text.x = element_text(angle=45),
        axis.text = element_text(size=12),
        axis.title = element_text(size=15))

ggsave(sprintf("%s/histogram_CV_d_land_rep_B.pdf",save.dir),g,width=10,height=10)

g <- ggplot(var_mat,aes(x=CV_deven))+
  geom_histogram()+
  facet_grid(plant~ra)+
  theme_bw()+
  theme(aspect.ratio = 1,
        axis.text.x = element_text(angle=45),
        axis.text = element_text(size=12),
        axis.title = element_text(size=15))

ggsave(sprintf("%s/histogram_CV_d_even_rep_B.pdf",save.dir),g,width=10,height=10)
