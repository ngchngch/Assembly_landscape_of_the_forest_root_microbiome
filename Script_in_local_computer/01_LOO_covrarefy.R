
library(renv)

library(vegan)

library(parallel)

set.seed(123)

detectCores()
nthred <- 128

#########################################################################
save.dir <- "Output/01_caverage_rarefaction"
dir.create(save.dir)
#==========rarefaction=================#

covrarefy <- function(df,readth=0,covth = 0,ncore){
  
  OTU_table <- df[rowSums(df)>readth,]
  s_comp <- list(NULL)
  for(i in 1:nrow(OTU_table)){
    s_comp[[i]]<-OTU_table[i,]
  }
  rareslopelist<-mclapply(s_comp,function(x){
    rareslope(x,1:(sum(x)-1))
  },mc.cores=ncore)
  
  getmincov<-unlist(mclapply(rareslopelist,function(x){
    x[length(x)]
  },mc.cores=ncore))
  
  #histogram(getmincov)
  maxslope <- max(getmincov)
  
  if(covth<(1-maxslope)){
    cov_th <- maxslope
  }else{
    cov_th <- 1-covth
  }

  #指定したカバレッジに到達した（＝傾きが指定値を下回る）瞬間のリード数をサンプルごとに採ってくる
  cvrfun<-function(x){min(which(x<=max(getmincov[getmincov<=cov_th])))+1} #関数を設定。上記で1を引いた分を足し戻す
  cvrrare<-unlist(mclapply(rareslopelist,cvrfun,mc.cores = ncore))　#lapply+unlistでベクトル形式にして一括で値を取得
  
  cvrrare2 <- cvrrare[getmincov<=cov_th]
  set.seed(123) #再現性をとるためにランダム変数を固定（数字は何でもいい）
  OTU_covrared<-rrarefy(OTU_table[getmincov<=cov_th,],cvrrare2) #得られたリード数に沿って各サンプルからリサンプリング
  
  
  OTU_covrared2 <- OTU_covrared[,colSums(OTU_covrared)>0]
  return(list(table=OTU_covrared2,cov_th=(1-cov_th)))
}

#########################################################################
#read table (seqtab before rarefaction)
df_f <- readRDS("Base_data/OTU_basedata_set/No_rarefy_sqtb_rootF.rds")
df_b <- readRDS("Base_data/OTU_basedata_set/seqOTUtab.rds")

tx_f <- readRDS("Base_data/Fungi/taxa_list_mod.rds")
tx_b <- as.data.frame(readRDS("Base_data/OTU_basedata_set/taxonomy_list_dada2.rds"))

info <- readRDS("Base_data/comp_sample_info.rds")

#remove non Fungi/BacteriaOrArcheae OTU
#Fungi
df_f2 <- df_f[which(rownames(df_f) %in% rownames(info)[which(info$plant != "Unidentified")]),
              which(colnames(df_f) %in% rownames(tx_f)[which(tx_f$Kingdom == "Fungi")])]

#Bacteria
df_b2 <- df_b[which(rownames(df_b) %in% rownames(info)[which(info$plant != "Unidentified")]),
              which(colnames(df_b) %in% rownames(tx_b)[which(tx_b$Kingdom %in% c("Bacteria","Archeae") &
                                                               tx_b$Order != "Chloroplast" &
                                                               tx_b$Family != "Mitochondria")])]


#read count threshold

readth <- 1000
df_f3 <- df_f2[rowSums(df_f2)>readth,]

df_b3 <- df_b2[rowSums(df_b2)>readth,]
dim(df_f3);dim(df_b3)

cons_samp <- intersect(rownames(df_f3),rownames(df_b3))
length(cons_samp)

df_f4 <- df_f3[cons_samp,]
df_b4 <- df_b3[cons_samp,]

##coverage based rarefaction
rdf_f <- covrarefy(df_f4,readth = 1000,ncore = nthred)
rdf_b <- covrarefy(df_b4,readth = 1000,ncore = nthred)

saveRDS(rdf_f,sprintf("%s/covrarefy_sqtb_fungi.rds",save.dir))
saveRDS(rdf_b,sprintf("%s/covrarefy_sqtb_bacteria.rds",save.dir))
#rdf_f <- readRDS("analysis_series/240424_caverage_rarefaction/covrarefy_sqtb_fungi.rds")
#rdf_b <- readRDS("analysis_series/240424_caverage_rarefaction/covrarefy_sqtb_bacteria.rds")

rtab_f <- rdf_f$table/rowSums(rdf_f$table)
rtab_b <- rdf_b$table/rowSums(rdf_b$table)

##check taxa (Family?? Genus??) abundance & occurence
#OTU
ocab_f_o <- cbind(RA=colMeans(rtab_f),Oc=colSums(rtab_f>0))
ocab_b_o <- cbind(RA=colMeans(rtab_b),Oc=colSums(rtab_b>0))

#Genus
rtf_g <- Taxa.mat(rtab_f,tx_f,"Genus")
rtb_g <- Taxa.mat(rtab_b,tx_b,"Genus")

ocab_f_g <- cbind(RA=colMeans(rtf_g),Oc=colSums(rtf_g>0))
ocab_b_g <- cbind(RA=colMeans(rtb_g),Oc=colSums(rtb_g>0))
#Family

rtf_f <- Taxa.mat(rtab_f,tx_f,"Family")
rtb_f <- Taxa.mat(rtab_b,tx_b,"Family")

ocab_f_f <- cbind(RA=colMeans(rtf_f),Oc=colSums(rtf_f>0))
ocab_b_f <- cbind(RA=colMeans(rtb_f),Oc=colSums(rtb_f>0))

#check threshhold
occth <- 100
sum(ocab_f_f[,2]>occth);sum(ocab_b_f[,2]>occth)
sum(ocab_f_g[,2]>occth);sum(ocab_b_g[,2]>occth)
sum(ocab_f_o[,2]>occth);sum(ocab_b_o[,2]>occth)

#select "Genus"
tg_tx_f <- setdiff(rownames(ocab_f_g[ocab_f_g[,2]>occth,]),"Unidentified")
tg_tx_b <- setdiff(rownames(ocab_b_g[ocab_b_g[,2]>occth,]),"Unidentified")

#remove target taxa
#Fungi
fb <- "Fungi"
dir.create(sprintf("%s/%s",save.dir,fb))
for(i in 1:length(tg_tx_f)){#i <-1
  show.progress(i,1:length(tg_tx_f))
  tg <- tg_tx_f[i]
  save.dir2 <- sprintf("%s/%s/%s",save.dir,fb,tg)
  dir.create(save.dir2)
  #use df before read count threshold
  df_pre <- df_f3[cons_samp,which(colnames(df_f3) %in% rownames(tx_f[tx_f$Genus != tg,]))]
  rtb <- covrarefy(df_pre,readth = 1000,ncore = nthred)
  
  if(length(setdiff(rownames(df_pre),rownames(rtb$table)))>0 ){
  png(sprintf("%s/target Taxa RA in delete samples.png",save.dir2),
      res=300,h=1200,w=1500)
  hist(Taxa.mat(df_f3/rowSums(df_f3),tx_f,"Genus")[setdiff(rownames(df_pre),rownames(rtb$table)),tg],
       main=sprintf("relative abundance of %s\nin deleted sample",tg),
       xlab=sprintf("relative abundance of %s",tg))
  dev.off()
  }
  
  png(sprintf("%s/target Taxa RA in samples.png",save.dir2),
      res=300,h=1200,w=1500)
  hist(Taxa.mat(df_f3/rowSums(df_f3),tx_f,"Genus")[rownames(rtb$table),tg],
       main=sprintf("relative abundance of %s in sample\nOccurence after rarefaction = %s",tg,
                    sum(Taxa.mat(rtab_f/rowSums(rtab_f),tx_f,"Genus")[rownames(rtb$table),tg]>0)),
       xlab=sprintf("relative abundance of %s",tg))
  dev.off()
  
     matlist <- list(Fungi=rtb$table/rowSums(rtb$table),Bacteria=rtab_b[rownames(rtb$table),])
    saveRDS(matlist,sprintf("%s/matrix_list_%s.rds",save.dir2,tg))
  }

#Bacteria
fb <- "Bacteria"
dir.create(sprintf("%s/%s",save.dir,fb))
for(i in 1:length(tg_tx_b)){#i <-1
  show.progress(i,1:length(tg_tx_b))
  tg <- tg_tx_b[i]
  save.dir2 <- sprintf("%s/%s/%s",save.dir,fb,tg)
  dir.create(save.dir2)
  #use df before read count threshold
  df_pre <- df_b3[cons_samp,which(colnames(df_b3) %in% rownames(tx_b[tx_b$Genus != tg,]))]
  rtb <- covrarefy(df_pre,readth = 1000,ncore = nthred)
  
  if(length(setdiff(rownames(df_pre),rownames(rtb$table)))>0 ){
    png(sprintf("%s/target Taxa RA in deleted samples.png",save.dir2),
        res=300,h=1200,w=1500)
    hist(Taxa.mat(df_b3/rowSums(df_b3),tx_b,"Genus")[setdiff(rownames(df_pre),rownames(rtb$table)),tg],
         main=sprintf("relative abundance of %s\nin deleted sample",tg),
         xlab=sprintf("relative abundance of %s",tg))
    dev.off()
    
  }
  
  png(sprintf("%s/target Taxa RA in samples.png",save.dir2),
      res=300,h=1200,w=1500)
  hist(Taxa.mat(df_b3/rowSums(df_b3),tx_b,"Genus")[rownames(rtb$table),tg],
       main=sprintf("relative abundance of %s in sample\nOccurence after rarefaction = %s",tg,
                    sum(Taxa.mat(rtab_b/rowSums(rtab_b),tx_b,"Genus")[rownames(rtb$table),tg]>0)),
       xlab=sprintf("relative abundance of %s",tg))
  dev.off()
  

    matlist <- list(Fungi=rtab_f[rownames(rtb$table),],Bacteria=rtb$table/rowSums(rtb$table))
    saveRDS(matlist,sprintf("%s/matrix_list_%s.rds",save.dir2,tg))

}



