
library(seqinr)
library(dada2)

# if(!requireNamespace("BiocManager")){
#   install.packages("BiocManager")
# }
# BiocManager::install("phyloseq")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("decontam")

library(decontam)


setwd('Base_data/Bioinfomatics')


source("../packages/01_1_function.R")
## -- Load data tables

outputdir <- "merged"

dir.create(outputdir)

####

tablist1 <- list.files("seqtab/lib1",full.names = T)

mtab1 <- mergeSequenceTables(tables=tabs1,repeats = "sum")


####

tablist2 <- list.files("seqtab/lib2",full.names = T)

mtab2 <- mergeSequenceTables(tables=tabs2,repeats = "sum")


####decontam

info <- read.csv("sample_info.csv",row.names = 1)
info <- data.frame(sample=union(rownames(mtab1),rownames(mtab2)),
                   nonw=NA
                   )

info$sample_type <- "sample"
info[grep("BLANK",info$sample),"sample_type"] <- "NegaCon"
rownames(info) <- info[,1]

test <-Decontam_df(mtab1,info,th=0.1)

test2 <-Decontam_df(mtab2,info,th=0.1)


dups <- intersect(rownames(test),rownames(test2))


tf <-sapply(dups,function(x){
  sum(test[x,])>sum(test2[x,])
})

lib1 <- test[which(!rownames(test) %in% names(tf[which(!tf)])),]

lib2 <- test2[which(!rownames(test2) %in% names(tf[which(tf)])),]

mtab <- mergeSequenceTables(tables=list(lib1,lib2),repeats = "error")




# Remove chimeras
seqtab <- mtab

seq.mat <- cbind(colnames(seqtab),sprintf('ASVB_%s', formatC(1:ncol(seqtab), width = nchar(ncol(seqtab)), flag = "0"))) 
write.fasta(as.list(seq.mat[,1]), seq.mat[,2], sprintf("%s/nonchim_seq.fasta", outputdir) )

seqtab2 <- seqtab
colnames(seqtab2) <- seq.mat[,2]
saveRDS(seqtab2, sprintf('%s/seqtab_rmChimera.rds', outputdir))
write.csv(cbind(sample=rownames(seqtab2), seqtab2), sprintf('%s/seqtab_rmChimera.csv', outputdir), row.names=FALSE)


###############################################################33
# check denoising result

track <- cbind(rowSums(mtab),rowSums(seqtab2))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("denoisedF","nonchim")
rownames(track) <- rownames(seqtab2)
#write.csv(track, sprintf('%s/denoising_result.csv', outputdir))

##################################################################################################
minident <- "0.97"
system2(command = vsearchpath, 
        args = c("--cluster_fast $input", sprintf("%s/nonchim_seq.fasta", outputdir),
                 "--id", minident,
                 "--mothur_shared_out", sprintf("%s/ASV_OTU_corestab_%s.txt", outputdir, minident),
                 "--centroids", sprintf("%s/OTUseq_%s.fasta", outputdir, minident),
                 "--msaout",  sprintf("%s/seqAlign_%s.txt", outputdir, minident) ) )

otu <- read.table(sprintf("%s/ASV_OTU_corestab_%s.txt", outputdir, minident),
                  header=TRUE, row.names=2)[,-c(1:2)]

## ========================================== ##

otutab <- matrix(0, ncol=ncol(otu), nrow=nrow(seqtab2),
                 dimnames=list(rownames(seqtab2), colnames(otu)))

for(i in 1:ncol(otu)){
  
  if( sum(otu[,i])>1 ){
    memberSeq <- rownames(otu)[which(otu[,i]>0)]
    otutab[,i] <- rowSums(seqtab2[,which(colnames(seqtab2) %in% memberSeq) ])
  }else{
    centroidSeq <- colnames(otu)[i]
    otutab[,i] <- seqtab2[, which(colnames(seqtab2) == centroidSeq) ]
  }
  
}

saveRDS(otutab, sprintf("%s/seqOTUtab.rds", outputdir))
write.csv(cbind(sample=rownames(otutab), otutab), sprintf('%s/seqOTUtab.csv', outputdir), row.names=FALSE)

saveRDS(otutab, "../OTU_basedata_set/seqOTUtab.rds")
write.csv(cbind(sample=rownames(otutab), otutab),
          '../OTU_basedata_set/seqOTUtab.csv', row.names=FALSE)

seqfasta <- read.fasta(sprintf('%s/OTUseq_0.97.fasta',outputdir))

seq.mat <- cbind(seqfasta,
                 sprintf('B_%s', formatC(1:length(seqfasta), width = nchar(length(seqfasta)), flag = "0"))) 
write.fasta(as.list(seq.mat[,1]), seq.mat[,2], sprintf("%s/OTUseq_0.97_rename.fasta", outputdir) )

write.fasta(as.list(seq.mat[,1]), seq.mat[,2], "../OTU_basedata_set/OTUseq_0.97_rename.fasta")

##################################################################################################

finish <- Sys.time()
print(finish-start); cat("\n")
##################################################################################################
sink("log_and_script/Version.txt", append=TRUE)
sessionInfo()
sink()

#############################
library(seqinr)

referencepath='reference/silva_nr99_v138_train_set_append_stdDNA.fa'

seqfasta <- read.fasta(sprintf('%s/OTUseq_0.97_rename.fasta',outputdir))

seqlist <- sapply(seqfasta, paste, collapse="")
taxa <- assignTaxonomy(seqlist, referencepath, multithread=TRUE)

#if(primer=="515f")taxa <- addSpecies(taxa,  sprintf('%s/silva_nr99_v138_train_set_append_stdDNA.fa', referencepath ) )


taxa.print <- taxa  # Removing sequence rownames for display only

## ===================================================== ##
# -- Compile table

# -- Compile table
if( max( sapply(strsplit(taxa.print[,"Kingdom"], "k__"), length), na.rm=TRUE ) > 1){
  taxa.print[,"Kingdom"] <- sapply(strsplit(taxa.print[,"Kingdom"], "k__"),  "[" ,2)
}

taxa.print[,"Phylum"] <- gsub("p__", "", taxa.print[,"Phylum"])
taxa.print[,"Class"] <- gsub("c__", "", taxa.print[,"Class"])
taxa.print[,"Order"] <- gsub("o__", "", taxa.print[,"Order"])
taxa.print[,"Family"] <- gsub("f__", "", taxa.print[,"Family"])
taxa.print[,"Genus"] <- gsub("g__", "", taxa.print[,"Genus"])
taxa.print[,"Species"] <- gsub("s__", "", taxa.print[,"Species"])

taxa.print[is.na(taxa.print)] <- "Unidentified"
taxa.print <- gsub("unidentified", "Unidentified", taxa.print)


rownames(taxa.print) <-names(seqlist)

write.table(cbind(ID=names(seqlist) ,taxa.print), sprintf("%s/taxonomy_list_dada2.txt",outputdir), row.names=FALSE, sep="\t", quote=FALSE)
saveRDS(cbind(ID=names(seqlist) ,taxa.print), sprintf("%s/taxonomy_list_dada2.rds",outputdir))

write.table(cbind(ID=names(seqlist) ,taxa.print), "../OTU_basedata_set/taxonomy_list_dada2.txt", row.names=FALSE, sep="\t", quote=FALSE)
saveRDS(cbind(ID=names(seqlist) ,taxa.print), "../OTU_basedata_set/taxonomy_list_dada2.rds")

sink("log_and_script/Version.txt", append=TRUE)
sessionInfo()
sink()
