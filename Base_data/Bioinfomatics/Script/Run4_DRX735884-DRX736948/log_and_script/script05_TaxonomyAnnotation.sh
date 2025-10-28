#!/bin/bash

####################################################################
## 											
## ------ Taxonomy Annotation with 
##	  RDP naive Bayesian classifier method      ------------- ##
##
## 					2022. 02. 28. by Fujita
####################################################################

inputdir=04_Denoising
outputdir=05_TaxonomyAnnotation
referencepath=/home/toju/Desktop/Scripts/noguchi/referenceDB/silva_nr99_v138_train_set_append_stdDNA.fa

clustering=TRUE

## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ##
piplineID="05_TaxonomyAnnotation"

####################################################################


## ============= Remove and Make directories ==================== ##

## --Remove older directory
if ls ${piplineID}_QCreport* >/dev/null 2>&1 ; then
	rm -r ${piplineID}
fi

if ls log_and_script/${piplineID}_log* >/dev/null 2>&1 ; then
	rm log_and_script/log${piplineID}
fi

## ------------------------------------------------------------- ##
## -- Making directory to save results
mkdir -p ${piplineID}

## ------------------------------------------------------------- ##
## -- Version check
echo "## ++++++ ${piplineID} +++++ ##" >> log_and_script/Version.txt

####################################################################

cat <<RRR > log_and_script/script${piplineID}.R
####################################################################
## 											
## ----- Taxonomy Annotation with 
##       RDP naive Bayesian classifier method  in R ------------- ##
##
## 					2022. 02. 28. by Fujita
####################################################################

## == Setting parameter
inputdir="${inputdir}"
referencepath="${referencepath}"

## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ ##
piplineID="$piplineID"

RRR

cat <<'RRR' >> log_and_script/script${piplineID}.R

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
RRR

Rscript log_and_script/script${piplineID}.R 2>&1 | tee log_and_script/log${piplineID}.txt


