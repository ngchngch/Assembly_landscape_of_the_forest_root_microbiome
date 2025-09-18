###########################
#shimapallete by Nishino

color <- c('#00005B','#FFFF00','#AB005D','#FF8DD7','#008203',
           '#242424','#00D087','#E425A6','#E8FFFF','#0000E1',
           '#049CC6','#3E5025','#FAAD3B','#AD5700',
           '#31A800','#4561FE',
           '#F06A00','#18E6FF','#69323D','#00657A','#F546FF',
           '#00FB00','#00BBFF','#A5788D','#FDBBB1','#5A7966','#303E50','#B8A55F',
           '#975B6E','#C2FBAD','#AA8BCC','#818200','#9EC3B7','#85512A','#2A1F44','#0D41F5',
           '#73976F','#FFE1FF','#880000','#CF004B','#2F7047','#08B3A2','#640058','#C5633F',
           '#BED1FD','#008F62','#F15F8D','#BDACB8','#6776B3','#712D85','#BA8F79','#00C000',
           '#B5ADEA','#A6DF9B','#754E97','#004700','#0044A6','#231C02','#89AE9B','#D6E4D2',
           '#FFDC00','#004F5D','#9975FF','#85856B','#696500','#362B36','#FFBDFF','#FC2200','#694849','#FF788C','#5C2C00','#00242B','#A40600','#99754E','#CB00A3','#706A55','#8083BB','#BD9512','#38340A','#C4BA67','#0D071B','#363636','#7E7291','#D7B6D0','#002100','#4E2F27','#E719FF','#1E332F','#876436','#D3C700','#3D200E','#F7E8B2','#003F2E','#62A5B5','#00E300','#2F88FF','#E6D0BE','#340033','#8A9DA9','#9DEBFF','#1B988D','#494763','#DFDFDF','#AB979F','#93D8D3','#50973D','#6E618A','#773A00','#A5826A','#005E4F','#1C3093','#446D74','#919405','#FAE8E0','#C8A800','#4D2D52','#4F6CA7','#9F985A','#3A8AA8','#211C2A','#79766A','#0C9000','#C4844A','#5C4022','#4DFF9C','#00A7FB','#6B6B6B','#DB916B','#66AF54','#B4B4B4','#474D92','#324535','#CBA48F','#C1FFFE','#58783B','#6C5F6A','#4D3F39','#3F611F','#3E3D78','#81B97F','#894861','#838A98','#C57300','#183257','#9BB0BD','#00A774','#CF89A0','#007300','#6E8253','#FFFF93','#8F7A7B','#00E993','#807205','#005900','#3D5DAD','#E6CCE8','#C5797E','#576448','#679092','#16C5B9','#4C4F3E','#B4B298','#B9ECCC','#400005','#BC70C5','#98907D','#8767F9','#BDD0C6','#6D808E','#918847','#75B9CC','#D6E4FF','#9C8AA8','#007D6F','#8791F6','#01D0F5','#D09EBB','#B4594A','#833E36','#A79F83','#73001F','#A3A0C9','#D0BDB2','#5E7486','#9BC3FF','#BD6295','#C393EB','#6C436F','#685200','#B2C0D8','#511431','#2DDBC4','#090733','#A26972','#97D15E','#450057','#62F2D5','#A86D00','#8AA6FB','#896B60','#CDCB9C','#675538','#91BD3C','#A74300','#6900A3','#95CB97','#435D61','#DEB590','#AB367A','#000087','#027BEE','#5746FF','#A15FBE','#BA022B','#0000B4','#007F9F','#58525D','#95579F','#D50400','#7F13BC','#825A52','#810050','#FACE9C','#DE3300','#A214A7','#A44D43','#980235','#FFC55E','#4201F8','#BFF825','#FFFFD6','#F69D68','#8C47C5','#E6DF7E','#E1556B','#FFA6C8','#FE9200','#FF72FF','#AAE557','#AA2CF6','#FA9599','#E2006C','#F26EC7','#FF7934','#D200F8','#FE4840','#F92FAA','#F5044A','#F0705E')
##############From fujita package "AnalysisHelper"
#' Make color coresspondece table

#' @description You can make color coresspondece table for ggplot2.

#' 

#' @param data matrix or dataframe containing sample x abundance

#' @param color a vector contained color information

#' @param others Not dominant color

#' @param specificName You can specify the color to specific species (like 'Unidentified')

#' @param specificColor For the specific species color

#' 

#' 

#' @examples

#' data(test)

#' col <- palettes()

#' colpal <- makeColPalette(data=test[,-1], 

#'                             color=col, 

#'                             others='grey30', 

#'                             specificName = 'unidentified', 

#'                             specificColor = 'grey90',

#'                             sort=sum)

#' 

#' @export







makeColPalette <- function (data, color = NULL, othersCol = "grey30", specificName = NULL, 

    specificColor = "grey90", sortFun = sum, na.rm = TRUE) 

{

    if (na.rm) {

        total.abundance <- apply(data, 2, sortFun, na.rm = TRUE)

    }

    else {

        total.abundance <- apply(data, 2, sortFun)

    }

    df <- data.frame(taxa = colnames(data), total.abundance = total.abundance, 

        stringsAsFactors = FALSE)

    df <- df[order(df$total.abundance, decreasing = TRUE), ]

    df$color <- NA

    if (!is.null(specificName)) {

        if (is.null(specificColor)) {

            specificColor <- "grey90"

        }

        if (any(df$taxa %in% specificName)) {

            specific <- which(df$taxa %in% specificName)

            df <- rbind(df[-specific, ], df[specific, ])

            specific <- -(length(specific) - 1):0 + nrow(df)

            df[c(1:length(color)), "color"] <- color

            df[specific, "color"] <- specificColor

            if (any(is.na(df$color))) {

                df$color[is.na(df$color)] <- othersCol

            }

            col <- df$color

            names(col) <- df$taxa

        }

        else {

            warning("Ignored specificName species because missing the species")

            df[c(1:length(color)), "color"] <- color

            if (any(is.na(df$color))) {

                df$color[is.na(df$color)] <- othersCol

            }

            col <- df$color

            names(col) <- df$taxa

        }

    }

    else {

        df[c(1:length(color)), "color"] <- color

        if (any(is.na(df$color))) {

            df$color[is.na(df$color)] <- othersCol

        }

        col <- df$color

        names(col) <- df$taxa

    }

    return(col[!is.na(names(col))])

}



################################################################################################
####
#### R script for Fujita (2019)
####
#### Functions for 01_1... R script
#### 2019.11.23Fujita
#### R 3.6.0
#### Set working directory of 'MAFF' folder -- setwd('~/Desktop/MAFF')
#### 
################################################################################################

Fillmerge <- function(gth,sep,fill,value,facet1,facet2){
  a <- 0
  sepnames <- unique(gth[,sep])
  fillnames <- unique(gth[,fill])
  if(is.na(facet1)==F){
    f1nam <- unique(gth[,facet1])
    if(is.na(facet2)==F){
      f2nam <- unique(gth[,facet2])
      mtab <- data.frame(NA,nrow=length(sepnames)*length(fillnames)*length(f1names)*length(f2names),ncol=5)
      colnames(mtab) <- c(sep,fill,value)
      for(i in 1:length(sepnames)){
        for(j in 1:length(fillnames)){
          a <- a+1
          mtab[a,c(1:2)] <-  c(sepnames[i],fillnames[j])
          mtab[a,value] <- sum(gth[which(gth[,sep]==sepnames[i] & gth[,fill]==fillnames[j]),value])
        }
      }
    }else{
      mtab <- data.frame(NA,nrow=length(sepnames)*length(fillnames)*length(f1names),ncol=4)
      colnames(mtab) <- c(sep,fill,value)
      for(i in 1:length(sepnames)){
        for(j in 1:length(fillnames)){
          a <- a+1
          mtab[a,c(1:2)] <-  c(sepnames[i],fillnames[j])
          mtab[a,value] <- sum(gth[which(gth[,sep]==sepnames[i] & gth[,fill]==fillnames[j]),value])
        }
      }
    }
  }else{
    mtab <- data.frame(NA,nrow=length(sepnames)*length(fillnames),ncol=3)
    colnames(mtab) <- c(sep,fill,value)
    for(i in 1:length(sepnames)){
      for(j in 1:length(fillnames)){
        a <- a+1
        mtab[a,c(1:2)] <-  c(sepnames[i],fillnames[j])
        mtab[a,value] <- sum(gth[which(gth[,sep]==sepnames[i] & gth[,fill]==fillnames[j]),value])
      }
    }
  }
  mtab
}


# -- Decide parameters
Load.information <- function(Region=NULL) {
	
	if(Region=='ITS1F'){
		pattern.f <- "__ITS1F.forward.fastq.gz"
		reference1 <- "function_refrence/sh_general_release_dynamic_10.10.2017_2.fasta" 
		fwd <- 'NNNNNNCTHGGTCATTTAGAGGAASTAA'
		
		return(list(pattern=pattern.f, reference=reference1, primer=fwd))
	}
	
	if(Region=='515f'){
		pattern.f <- "515f.forward.fastq.gz"
		reference1 <- "function_refrence/silva_nr_v132_train_set_with_std.fasta" 
		fwd <- 'NNNNNNGTGYCAGCMGCCGCGGTAA'
		
		return(list(pattern=pattern.f, reference=reference1, primer=fwd))
	}
}


Check.primer = function(fwd= fwd,
						path=path,
						pattern=pattern){

  allOrients <- function(primer) {
    # Create all orientations of the input sequence
    require(Biostrings)
    dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
    orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
    return(sapply(orients, toString))  # Convert back to character vector
  }
  FWD.orients <- allOrients(fwd)

  fnFs <- sort(list.files(path, pattern = pattern, full.names = TRUE))
  fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
  filterAndTrim(fwd=fnFs, filt=fnFs.filtN, maxN = 0, multithread = TRUE)

  primerHits <- function(primer, fn) {
   # Counts number of reads in which the primer is found
   nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
    return(sum(nhits > 0))
  }

  tmp <- sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]])
  tmp
}


show.progress <- function(i, x){
  # i??x?ɂ?for?̃??[?v?ϐ??ƃ??X?g???^????
  prg <- round(which(x==i)/length(x)*100)
  done   <- paste(rep('#', round(prg/2.5)),    collapse='')
  remain <- paste(rep('-', 40-round(prg/2.5)), collapse='')
  prg.bar <- paste('|', done, remain, '|  ', as.character(prg), '%', sep='')
  message('\r', prg.bar, appendLF=FALSE)
}

############################################################################
####
#### R script for Fujita (2019)
####
#### Function for 02_data_shaping scripts
#### 2019.11.04 Fujita
#### R 3.6.0
#### Set working directory of 'MTS' folder -- setwd('~/Desktop/Microbiome_TimeSeries/MTS')
####
############################################################################

# -- Make each taxnomy matrix
Taxa.mat <- function(x,y, taxaLabel){
  
  
  # x is sample-ASV matrix, y is taxonomy table, taxaLabel is a taxanomy level
  # You need to change rownames(y) to ASV/OTU name
  
  Taxa.ID <- colnames(x)
  matTT <- t(x)
  mergeMat<- cbind(y[as.matrix(rownames(matTT)),], matTT)
  
  TaxaList <- unique(mergeMat[, taxaLabel])
  mat.Genus <- matrix(NA, nrow=length(TaxaList), ncol=nrow(x))
  rownames(mat.Genus) <- TaxaList
  colnames(mat.Genus) <- rownames(x)
  
  for (i in 1:length(TaxaList)) {
    if(length(which(mergeMat[,taxaLabel] == TaxaList[i]))!=1){
      mat.Genus[i,] <- colSums(matTT[which(mergeMat[,taxaLabel] == TaxaList[i]),])
    }else{
      mat.Genus[i,] <- matTT[which(mergeMat[,taxaLabel] == TaxaList[i]),]
    }
  }
  
  t(mat.Genus[order(rowSums(mat.Genus), decreasing=TRUE),])
}

# -- Make each taxnomy matrix
Taxa.mat2 <- function(x, y, taxaLabel){
  
  # x is sample-ASV matrix, y is taxonomy table, taxaLabel is a taxanomy level
  # You need to change rownames(y) to ASV/OTU name
  
  Taxa.ID <- colnames(x)
  matTT <- t(x)
  mergeMat<- cbind(y[as.matrix(rownames(matTT)),], matTT)
  
  TaxaList <- unique(mergeMat[, taxaLabel])
  mat.Genus <- matrix(NA, nrow=length(TaxaList), ncol=nrow(x))
  rownames(mat.Genus) <- TaxaList
  colnames(mat.Genus) <- rownames(x)
  
  for (i in 1:length(TaxaList)) {
    
    if(length(which(mergeMat[,taxaLabel] == TaxaList[i]))!=1){
      mat.Genus[i,] <- colSums(matTT[which(mergeMat[,taxaLabel] == TaxaList[i]),])/colSums(matTT)
    }else{
      mat.Genus[i,] <- matTT[which(mergeMat[,taxaLabel] == TaxaList[i]),]/colSums(matTT)
    }
    
  }
  
  t(mat.Genus[order(rowSums(mat.Genus), decreasing=TRUE),]*100)
}


# -- Extract Legend 
g_legend<-function(a.gplot){ 
  
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)}

# -- Show palette of RColorBrewer
rcolor <- function(x){display.brewer.all()}

# -- Extract each treatment
Subset <- function(data,info,row){ 
  
  k <- colnames(info) # get colnames to subset each treatment
  for(i in 1:length(k)) { data <- data[which(data[,k[i]] == info[row,k[i]]),] }
  data
  
}

###################
# -- Colorpalette
lib <- 'RColorBrewer';library(package = lib, character.only=TRUE);packageVersion(lib) 
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col.vector = unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))


palettes <- function(x){ col.vector[1:length(unique(x))]}

# -- Extract each treatment
Subset <- function(data,info,row){ 
  
  k <- colnames(info) # get colnames to subset each treatment
  for(i in 1:length(k)) { data <- data[which(data[,k[i]] == info[row,k[i]]),] }
  data
  
}

##########https://qiita.com/taka-ohi/items/dcc7115cc1b5ccda65a7
sig <- function(a,bin=F,star=T) {
  if(star==T){
    
    if(bin==T){
      if (a >= 0.025) {
        return("")
      } else {
        
        if ((a < 0.025)&&(a >= 0.005)) {
          return("*")
        } else {
          if ((a < 0.005)&&(a >= 0.0005)) {
            return("**")
          } else return("***")
        }
      }
    }else{
      if (a >= 0.05) {
        return("")
      } else {
        
        if ((a < 0.05)&&(a >= 0.01)) {
          return("*")
        } else {
          if ((a < 0.01)&&(a >= 0.001)) {
            return("**")
          } else return("***")
        }
      }
    } 
  }else{
    
    if(bin==T){
      if (a >= 0.025) {
        return("")
      } else {
        
        if ((a < 0.025)&&(a >= 0.005)) {
          return("p<0.05")
        } else {
          if ((a < 0.005)&&(a >= 0.0005)) {
            return("p<0.01")
          } else return("p<0.001")
        }
      }
    }else{
      if (a >= 0.05) {
        return("")
      } else {
        
        if ((a < 0.05)&&(a >= 0.01)) {
          return("p<0.05")
        } else {
          if ((a < 0.01)&&(a >= 0.001)) {
            return("p<0.01")
          } else return("p<0.001")
        }
      }
    }
  }}
  
  
  ######################################################################
##noguchi 20221114 to make color pallete for taxonomic barplot 

Colpal_clt <- function(seqtab,
                       class="Phylum",
                       class2,
                       taxalist,
                       palette,
                       unident_tag = "Unidentified",
                       th=FALSE){
  unident_cls <- unique(taxalist[grep(unident_tag,taxalist[,class]),class2])
  
  unident_cls2 <- setdiff(unique(taxalist[grep(unident_tag,taxalist[,class2]),class2]),unident_cls)
  
  unident_list <- taxalist[which(taxalist[,class2] %in% unident_cls2),c(class,class2)]
  
  txtab <- Taxa.mat(seqtab,taxalist,class2)
  #make color pallete
  #color <- makeColPalette(txtab[,setdiff(colnames(txtab)[order(setdiff(colSums(txtab),"unidentified"),decreasing=T)],"unidentified")],color=palettes(x=30))
  taxorder <-c()
  for(i in 1:(length(unique(setdiff(taxalist[,class],unident_cls))))){#i <- 1
    tlist <- setdiff(unique(taxalist[which(taxalist[,class] == unique(setdiff(taxalist[,class],unident_cls))[i]),class2]),unident_cls2)
    utlist <- intersect(unique(taxalist[which(taxalist[,class] == unique(setdiff(taxalist[,class],unident_cls))[i]),class2]),unident_cls2)
    
    taxorder <-c(taxorder,
                 colnames(txtab)[which(colnames(txtab) %in% tlist)][order(colSums(txtab)[which(colnames(txtab) %in% tlist)],decreasing = TRUE)],
                 colnames(txtab)[which(colnames(txtab) %in% utlist)][order(colSums(txtab)[which(colnames(txtab) %in% utlist)],decreasing = TRUE)])
  }
  taxorder2 <- c(taxorder,colnames(txtab)[which(colnames(txtab) %in% unident_cls)][order(colSums(txtab)[which(colnames(txtab) %in% unident_cls)],decreasing = TRUE)])
  
  
  color <- makeColPalette(txtab[,setdiff(colnames(txtab),unident_cls)],color=palette)
  
  if(length(unident_cls)<5){
    unident_col <- paste0("gray",seq((80-10*(length(unident_cls)-1)),80,by=10))
  }else{print("too many unidentfied taxa in class")}
  
  names(unident_col) <- unident_cls
  color2 <- c(color,unident_col)[taxorder2]
  return(color2)
}

#####################
#covarage base rarefaction 
#modified 

covrarefy <- function(df,readth,covarage_th=0.99,ncore){
  covth <- 1-covarage_th
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
  if(max(getmincov)<covth){
    cov_th <- max(getmincov)
  }else{
    cov_th <- max(covth)
  }
  
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
  return(list(tab=OTU_covrared2,
              min_cov=cov_th))
}

######################################
#remove contami OTU
Decontam_df <- function(df,info,negacon="NegaCon",th=0.1){
  require("phyloseq")
  taxa <- matrix("dami",nrow=ncol(df),ncol=4)
  rownames(taxa) <- colnames(df)
  ps <- phyloseq(otu_table(df, taxa_are_rows = FALSE),
                 sample_data(info[rownames(df),]),
                 tax_table(as.matrix(taxa[colnames(df),])))
  sample_data(ps)$is.neg <- sample_data(ps)$sample_type == "NegaCon"
  
  contamdf.prev05 <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=th)
  
  no_contami <- df[,!contamdf.prev05$contaminant]
  return(no_contami)
}

#####################################

th_str <- function(tab,th1,th2,tf){
  data.frame(name=rownames(tab)[which(ifelse(tf,rowSums(tab)>th1,rowSums(tab)>th2))],
             reads=rowSums(tab)[which(ifelse(tf,rowSums(tab)>th1,rowSums(tab)>th2))])
}

###############
rowSelect <- function(table,colName,sp_Names,invert=F){
  if(invert==F){
    
    table[which(table[,colName] %in% sp_Names),]
  }else{
    table[which(!table[,colName] %in% sp_Names),]
  }
}




#####
ggBar <- function(df,taxa,info,class,
                  lim=0,spNames=NA,Unident_Names="Unidentified",
                  color,order=NA){
  tax_tab <- Taxa.mat(df[,which(colnames(df) %in% rownames(taxa_f4))],taxa_f4,class)
  
  if(lim != 0){
    tax_tab2 <- tax_tab[,setdiff(colnames(tax_tab),Unident_Names)]
    tax_tab3 <- tax_tab2[,order(colSums(tax_tab2),decreasing=T)[1:lim]]
    tax_tab4 <- cbind(tax_tab3,
                      Unidentified=tax_tab[,Unident_Names],
                      others=rowSums(tax_tab2[,setdiff(colnames(tax_tab2),colnames(tax_tab3))]))
    
    if(is.na(spNames[1])){
      taxa <- c(setdiff(colnames(tax_tab4),c(Unident_Names,"others")),Unident_Names,"others")
    }else{
      taxa <- c(setdiff(colnames(tax_tab4),c(spNames,Unident_Names,"others")),spNames,Unident_Names,"others")
    }
    
  }else{
    tax_tab4 <- tax_tab
    if(is.na(spNames[1])){
      taxa <- c(setdiff(colnames(tax_tab4),c(Unident_Names)),Unident_Names)
    }else{
      taxa <- c(setdiff(colnames(tax_tab4),c(spNames,Unident_Names)),spNames,Unident_Names)
    }
    
  }
  #apply(df_f3[setdiff(rownames(tax_tab),rownames(info)),],1,function(x){taxa_f[names(x[x>0]),]})
  
  #class <- "ASV"#"Order"#"Family"#"Genus"
  
  #tax_tab4 <- log10(tax_tab4_)*tax_tab4_/rowSums(tax_tab4_)
  gtab <- gather(cbind(info[rownames(tax_tab4),],as.data.frame(tax_tab4)),class,abundance,-c(1:(ncol(info))))
  
  
  colfng <- cbind(taxa=taxa,color=c(color[1:(length(taxa)-2)],"gray80","gray40"))
  rownames(colfng) <- colfng[,1]
  
  gtab$z <- factor(gtab$class,levels=taxa)
  
  d <- vegdist(tax_tab4)
  
  a <- hclust(d)
  
  if(!is.na(order[1])){
    gtab$ID2 <- factor(gtab$ID,levels=order)
  }else{
    nam <- rownames(tax_tab4)[a$order]
    gtab$ID2 <- factor(gtab$ID,levels=nam)
    return(nam)
  }
  
  
  gg.bar.asv <- ggplot()+
    geom_bar(data=gtab, 
             aes(x=ID2, y=abundance, fill=z,color=z),position="fill",
             stat='identity',width=0.95) +
    #scale_y_log10()+
    labs(title = paste("Top", lim, class,  "Bar plot"),
         fill = class, color = class,y = "relative read count",x="sample")+ 
    theme_light()+
    theme(axis.text.x = element_blank(),axis.text.y = element_text(size=10),
          legend.text = element_markdown(),
          legend.key.size=unit(3, "mm")#,strip.background = element_blank()
          #,strip.text = element_blank()
    )+
    #facet_wrap(~sampling,scale="free")+
    coord_cartesian(ylim = c(0.04,0.98))+
    scale_fill_manual(values=c(colfng[levels(gtab$z),2]),
                      guide=guide_legend(ncol=1))+
    scale_color_manual(values=c(colfng[levels(gtab$z),2]),
                       guide=guide_legend(ncol=1))
  
  return(gg.bar.asv)
  
}

sampleOrder <- function(df,taxa,class){
  n_class <- length(class)
  tab <- matrix(1,nrow=nrow(df),ncol=1)
  for(i in 1:n_class){
    tab <- cbind(tab,Taxa.mat(df,taxa,class[i]))
  }
  
  d <- vegdist(tab)
  
  a <- hclust(d)
  
  nam <- rownames(df)[a$order]
  return(nam)
}


#function

blockSample <- function(mat,el,name){
  block_rand <- c()
  names <- c()
  randnames <- c()
    for(i in 1:length(unique(el))){#i <- 1
      names <- c(names,name[which(el %in% unique(el)[i])])
      randnames <- c(randnames,sample(name[which(el %in% unique(el)[i])]))
    }
  if(is.vector(mat)){
    block_rand <- mat[names,]
    names(block_rand) <- randnames
    block_rand <- block_rand[names(mat),]
  }else{
    block_rand <- mat[names,]
    rownames(block_rand) <- randnames
    block_rand <- block_rand[rownames(mat),]
  }
 
  return(block_rand)
}

##########################################
#detect outliers

detectOuter <- function(vec,p_value=0.05,method="BH",with_P=TRUE){
  #vec <- mdf$Chisq
  require("outliers")
  df <- vec
  t <- grubbs.test(df)
  outer <- c()
  p <- t$p.value
    
  
  while(t$p.value<p_value){
    r <- sum(vec==df[which.max(abs(df))])
    while(r>0 && t$p.value<p_value){
     #r <-1
      outer <- c(outer,which(vec==df[which.max(abs(df))])[r])
      df <- vec[-outer]
      t <- grubbs.test(df)
      p <- c(p,t$p.value)
      r <- r-1
     
    }}
  r_end <- sum(vec==df[which.max(abs(df))])
  outer <- c(outer,which(vec==df[which.max(abs(df))])[r_end])
  
  if(with_P){
    return(list(which_BH=c(outer[1:(which(p.adjust(p,method=method)>=p_value)[1]-1)]),
                which_all=c(outer),
                p_BH=c(p.adjust(p,method=method)),
                raw_p=c(p)))
  }else{
    return(outer[1:(which(p.adjust(p,method=method)>=p_value)[1]-1)])
  }

}

###################
check_pal <- function(col,n_color){#col <- color;n_color=20
  plot(1:n_color,1:n_color,col=col[1:n_color], pch = 19)
}

##############3

#merge row abundant OTU as Others

mergeOther <- function(df,lim,Unident_names="Unidentified",Others_names="Others"){
  df2 <- df[,which(colnames(df) != Unident_names)]
  major_tab <- df2[,order(colSums(df2),decreasing = TRUE)[1:lim]]
  others <- rowSums(df2[,setdiff(colnames(df2),colnames(major_tab))])
  
  if(sum(colnames(df) == Unident_names)>0){
    df3 <- cbind(major_tab,others,
                 Unidentified=df[,which(colnames(df) == Unident_names)])
    colnames(df3) <- c(colnames(major_tab),Others_names,Unident_names)
    
  }else{
    df3 <- cbind(major_tab,others)
    colnames(df3) <- c(colnames(major_tab),Others_names)
    
  }
  return(df3)
}


################


multiComp_SD <- function(data,groups,value){
  
  # #input
  # data <- df_cor3
  # groups <- df_cor3$Guild2
  # value <- "jsdm_low"
  # ###
  alp <- c("a","b","c","d","e","f","g","h","i","j")
  
  gdf <- data.frame(group=unique(groups),
                    med=NA,alpha="")
  
  for(i in 1:nrow(gdf)){
    gdf[which(gdf[,1] == gdf[i,1]),2] <- median(data[which(groups == gdf[i,1]),value])
  }
  
  set.seed(0)
  p_exact2 = NSM3::pSDCFlig(data[,value],
                            groups)
  
  
  lab_mat <- matrix(unlist(strsplit(p_exact2$labels," - ")),ncol=2,byrow = T)
  colnames(lab_mat) <- c("Group1","Group2")
  res <- data.frame(labels=p_exact2$labels,lab_mat,
                    P=p_exact2$p.val)
  def_gdf <- gdf
  
  if(!all(res$P>0.05)){
    
    i <- 1
    r <-1
    while(r>0){
      
      def_gdf[which(def_gdf$group == gdf[which.max(gdf[,"med"]),"group"]),
              "alpha"] <- alp[i]
      
      top_g <- as.character(def_gdf[which(def_gdf$group == gdf[which.max(gdf[,"med"]),"group"]),
                                    "group"])
      
      sel_res <- res[grep(top_g,res$labels),]
      
      nosig1 <- setdiff(unique(c(sel_res[sel_res[,"P"]>0.05,"Group1"],
                                 sel_res[sel_res[,"P"]>0.05,"Group2"])),top_g)
      
      def_gdf[which(def_gdf[,"group"] %in% nosig1),"alpha"] <- paste(def_gdf[which(def_gdf[,"group"] %in% nosig1),"alpha"],alp[i],sep="")
      
      gdf <- gdf[which(!gdf[,"group"] %in% c(nosig1,top_g)),]
      r <- sum(!gdf[,"group"] %in% c(nosig1,top_g))
      i <- i+1
      
    }
    
  }
  rownames(def_gdf) <- def_gdf$group
  return(list(result_test=res,alphabet=def_gdf))
}
