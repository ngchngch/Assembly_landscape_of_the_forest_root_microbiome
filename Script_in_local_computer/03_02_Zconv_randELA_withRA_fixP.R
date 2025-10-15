library(ggplot2)
library(ggstar)

setwd("/Volumes/8TBHDD_NGCH/sugadaira_bacteria_2023/02_ELA_241211")
set.seed(1234)

dir_03_01="Output_supercomputer/03_01_ELA_withRA_4step_250227"
dir_03_02="Output_supercomputer/03_02_summarize_randELA_withRA_fixP"







#########################################################################
save.dir <- "03_02_02_Zconv_dland_fixP"
dir.create(save.dir)

########################################################################

org_files <- list.files(path=dir_03_01,
                  pattern = "ELA_SSprob_diff_", full.names = TRUE, recursive = TRUE)
names(org_files) <- sapply(strsplit(org_files,"/"),function(x){x[3]})
rd_summary <- list.dirs(path=dir_03_02,
                        recursive = F)
  
gns <- unique(list.files(rd_summary))
　　　　　　　　　　　　
zres <- NULL
zres2 <- NULL
nofiles <- c()
nofiles2 <- c()

f_files <- list.files(pattern = "^d",rd_summary[grep("Fungi",rd_summary)], full.names = TRUE,recursive = TRUE)
fnam <- sapply(strsplit(f_files,"/"),function(x){x[3]})

p_files <- list.files(pattern = "^d",rd_summary[grep("Prokaryota",rd_summary)], full.names = TRUE,recursive = TRUE)
pnam <- sapply(strsplit(p_files,"/"),function(x){x[3]})

dir.create(file.path(save.dir,"hist"))
dir.create(file.path(save.dir,"hist","Prokaryota"))
dir.create(file.path(save.dir,"hist","Fungi"))

for (i in c(1:length(gns))){#i <- 153
  cat(sprintf("%s %s / %s \n",gns[i],i,length(gns)))
  fb <- "Fungi"
  
  dir.create(file.path(save.dir,"hist",fb,gns[i]))
  org_files2 <- org_files[grep(fb,org_files)]
  
  if(length(org_files2[grep(gns[i],org_files2,fixed = TRUE)])>0){
    
    org <- readRDS(org_files2[gns[i]])
    
    f2 <- f_files[fnam == gns[i]]
    
    rland <- do.call(cbind,lapply(f2[grep("dland_1_3000_6001_10500.rds|dland_1_10500",f2)],readRDS))[,1:10000]
    
    # hist(x[-1],breaks=100)
    # x[1]
    # 0.103579+0.2062075*0.1700627
    zdland <- apply(cbind(org[,"d_land"],rland),1,function(x){#x <- cbind(org[,"d_land"],rland)[5,]
      (x[1]-mean(x[-1]))/sd(x[-1])
    })
    
    pdland <- apply(cbind(org[,"d_land"],rland),1,function(x){#x <- cbind(org[,"d_even"],rland)[1,]
      sum(x[-1] >= x[1])/length(x[-1])
    })
    
    
    reven <- do.call(cbind,lapply(f2[grep("deven_1_3000_6001_10500.rds|deven_1_10500",f2)],readRDS))[,1:10000]
    
    zdeven <- apply(cbind(org[,"d_even"],reven),1,function(x){#x <- cbind(org[,"d_even"],rland)[1,]
      (x[1]-mean(x[-1]))/sd(x[-1])
    })
    
    pdeven <- apply(cbind(z=zdeven,org=org[,"d_even"],reven),1,function(x){
      if(!is.na(x[1])){
        if(x[1]>0){
          sum(x[-c(1,2)] >= x[2])/length(x[-c(1,2)])
        }else{
          sum(x[-c(1,2)] <= x[2])/length(x[-c(1,2)])
        }}else{
          NA
        }
    })
    
    zres <- rbind(zres,data.frame(org[,c("Taxa","plant","ra")],z_land=zdland,p_land=pdland,z_even=zdeven,p_even=pdeven))
    
  }else{
    nofiles <- c(nofiles,gns[i])
  }
  
  for(ja in 1:nrow(org)){
    if(!is.na(rland[ja,1])){
      png(sprintf("%s/hist/%s/%s/hist_test_%s_%s.png",save.dir,fb,gns[i],org[ja,"plant"],org[ja,"ab"]))  
      hist(rland[ja,],breaks=100,
           xlab="randam_dland",col="blue")
      abline(v=org[ja,"d_land"],col="red")
      dev.off()
    }    
  }
  
  
  fb <- "Prokaryota"
  dir.create(file.path(save.dir,"hist",fb,gns[i]))
  org_files2 <- org_files[grep(fb,org_files)]
  
  if(length(org_files2[grep(gns[i],org_files2,fixed = TRUE)])>0){
    
    org <- readRDS(org_files2[gns[i]])
    f2 <- p_files[pnam == gns[i]]
    
    rland <- do.call(cbind,lapply(f2[grep("dland",f2)],readRDS))[,1:10000]
    
    zdland <- apply(cbind(org[,"d_land"],rland),1,function(x){
      (x[1]-mean(x[-1]))/sd(x[-1])
    })
    
    pdland <- apply(cbind(org[,"d_land"],rland),1,function(x){
      sum(x[-1] >= x[1])/length(x[-1])
    })
    
    
    reven <- do.call(cbind,lapply(f2[grep("deven",f2)],readRDS))[,1:10000]
    
    zdeven <- apply(cbind(org[,"d_even"],reven),1,function(x){
      (x[1]-mean(x[-1]))/sd(x[-1])
    })
    
    pdeven <- apply(cbind(z=zdeven,org=org[,"d_even"],reven),1,function(x){
      if(!is.na(x[1])){
        if(x[1]>0){
          sum(x[-c(1,2)] >= x[2])/length(x[-c(1,2)])
        }else{
          sum(x[-c(1,2)] <= x[2])/length(x[-c(1,2)])
        }}else{
          NA
        }
    })
    
    zres2 <- rbind(zres2,data.frame(org[,c("Taxa","plant","ra")],z_land=zdland,p_land=pdland,z_even=zdeven,p_even=pdeven))
    
  }else{
    nofiles2 <- c(nofiles2,gns[i])
  }
  
  for(ja in 1:nrow(org)){
    if(!is.na(rland[ja,1])){
    png(sprintf("%s/hist/%s/%s/hist_test_%s_%s.png",save.dir,fb,gns[i],org[ja,"plant"],org[ja,"ab"]))  
    hist(rland[ja,],breaks=100,
         xlab="randam_dland",col="blue")
    abline(v=org[ja,"d_land"],col="red")
    dev.off()
    }
  }
  
  }

zres$p_land[zres$p_land==0] <- 1/10000
zres$p_even[zres$p_even==0] <- 1/10000
zres2$p_land[zres2$p_land==0] <- 1/10000
zres2$p_even[zres2$p_even==0] <- 1/10000

saveRDS(zres,file.path(save.dir,"zvalue_Fungi.rds"))
saveRDS(zres2,file.path(save.dir,"zvalue_Prokaryote.rds"))
