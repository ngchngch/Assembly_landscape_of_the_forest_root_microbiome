####################################################################
## 											
## ----- Taxonomy Annotation with 
##       RDP naive Bayesian classifier method  in R ------------- ##
##
## 					2022. 02. 28. by Fujita
####################################################################

## == Setting parameter
inputdir="04_Denoising"
referencepath="/home/toju/Desktop/Scripts/noguchi/referenceDB/silva_nr99_v138_train_set_append_stdDNA.fa"

## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ##
piplineID="05_TaxonomyAnnotation"


####################################################################
dir.create(piplineID, showWarnings = FALSE)


## ===================================================== ##
## -- Load packages
library(dada2)
library(seqinr)

## -- Load data tables
seqtab <- as.matrix(readRDS( sprintf("%s/seqtab_rmChimera.rds",inputdir)))
cols <- colnames(seqtab)
seqfasta <- read.fasta(sprintf("%s/nonchim_seq.fasta",inputdir))

## ===================================================== ##
# -- Assign taxonomy
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


rownames(taxa.print) <- taxa.print[,1]

write.table(cbind(ID=cols ,taxa.print), sprintf("%s/taxonomy_list.txt",piplineID), row.names=FALSE, sep="\t", quote=FALSE)
saveRDS(cbind(ID=cols ,taxa.print), sprintf("%s/taxonomy_list.rds",piplineID))

sink("log_and_script/Version.txt", append=TRUE)
sessionInfo()
sink()
