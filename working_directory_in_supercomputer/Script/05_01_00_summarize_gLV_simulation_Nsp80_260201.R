set.seed(1234)
#########################################################################
save.dir <- "05_01_00_summarize_gLV_simulation_Nsp80_260201"
dir.create(save.dir)

########################################################################



#dir_05_01 <- "05_01_gLV_simulation_Nsp40_260107"

res <- NULL
for(connect in c(0.2,0.3,0.4)){
  #connect <- 0.2
  lf <- list.files(dir_05_01, pattern = sprintf("community_matrix_conn%s",connect),
                   recursive=TRUE,full.names = TRUE)
  
  dlf <- sapply(lf,function(x){
    x2 <- readRDS(x)
    return(colMeans(x2$bin_mat))
    })
  freq <- colSums(dlf>0.1&dlf<0.9)
  top3 <- freq[order(freq,decreasing=TRUE)][1:5]
  
  
  res <- rbind(res,data.frame(connect,as.numeric(gsub("seed","",sapply(strsplit(names(top3),"/"),"[",2))),
                         names(top3),
                         top3
                         ))
}

saveRDS(res,file.path(save.dir,"top3_seed_across.rds"))
