##read packages
lib <- 'tidyr';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'vegan';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'stringr';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'parallel';library(package = lib, character.only=TRUE);packageVersion(lib)
lib <- 'dplyr';library(package = lib, character.only=TRUE);packageVersion(lib)


####
setwd('/media/evolution/Transcend/data/sugadaira_bacteria_2023/240404_Megre')

source("/media/evolution/Transcend/data/Sugadaira_sequence/Final_merge/analysis_sequence/01_1_function.R")
## -- Load data tables


##set directory

save.dir <- "01_rarefaction"

dir.create(save.dir)

#==========read data====#

st <- readRDS("merged/OTU_basedata_set/seqOTUtab.rds")
tx <- readRDS("merged/OTU_basedata_set/taxonomy_list_dada2.rds")

#remove mtcondria & chloroplast & non bac/archeae otu
mt_cl <- rownames(tx[which(tx[,"Order"] != "Chloroplast" &
                             tx[,"Family"] != "Mitochondria" &
                             tx[,"Kingdom"] %in% c("Bacteria", "Archaea")), ])

st2 <- st[,which(colnames(st) %in% mt_cl)]

#separate root/soil sample
st_r <- st2[grep("NF",rownames(st2)),]
st_s <- st2[grep("NS",rownames(st2)),]

#################
#rarefactioncurve
par(mar = c(6.0, 12.0, 4.1, 2))
par(mgp = c(2,1,1))
png(sprintf("%s/rarecurve_Root.png",save.dir),
    height = 3200, width = 3000, res = 600)
rarecurve(st_r[sample(rownames(st_r),size=200),],label=F,col=color,
          xlab="Read counts",y="Number of OTUs",
          cex.lab  = 1.2, 
          cex.axis = 1.1)

dev.off()

par(mar = c(6.0, 12.0, 4.1, 2))
par(mgp = c(2,1,1))
png(sprintf("%s/rarecurve_soil.png",save.dir),
    height = 3200, width = 3000, res = 600)
rarecurve(st_s,label=F,col=color,
          xlab="Read counts",y="Number of OTUs",
          cex.lab  = 1.2, 
          cex.axis = 1.1)

dev.off()

#==========rarefaction=================#


covrarefy <- function(df,readth=0,ncore){
  
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
  cov_th <- max(getmincov)
  #histogram(getmincov[getmincov<=cov_th])
  
  #pass_sample <- rownames(OTU_table)[getmincov<=cov_th]
  
  #length(pass_sample)/nrow(OTU_table)
  
  
  #指定したカバレッジに到達した（＝傾きが指定値を下回る）瞬間のリード数をサンプルごとに採ってくる
  cvrfun<-function(x){min(which(x<=max(getmincov[getmincov<=cov_th])))+1} #関数を設定。上記で1を引いた分を足し戻す
  cvrrare<-unlist(mclapply(rareslopelist,cvrfun,mc.cores = ncore))　#lapply+unlistでベクトル形式にして一括で値を取得
  
  cvrrare2 <- cvrrare[getmincov<=cov_th]
  set.seed(123) #再現性をとるためにランダム変数を固定（数字は何でもいい）
  OTU_covrared<-rrarefy(OTU_table[getmincov<=cov_th,],cvrrare2) #得られたリード数に沿って各サンプルからリサンプリング
  
  
  OTU_covrared2 <- OTU_covrared[,colSums(OTU_covrared)>0]
  return(list(table=OTU_covrared2,cov_th=(1-cov_th)))
}

################
#===root=========#
rdf_root <- covrarefy(st_r[rowSums(st_r) > 2000 , ],ncore=128)

saveRDS(rdf_root,sprintf("%s/covrarefy_sqtb_rootB.rds",save.dir))
#============soil=====================#
rdf_soil <- covrarefy(st_s[rowSums(st_s) > 5000, ],ncore=128)

saveRDS(rdf_soil,sprintf("%s/covrarefy_sqtb_soilB.rds",save.dir))


