
#!/bin/bash
## ============================================= ##
## -- Making scripts
sh /home/toju/Desktop/Scripts/noguchi/Filtering2Annotation_support_single/02_cutadapt_R.sh NNNNNNGTGYCAGCMGCCGCGGTAA NNNNNNGGACTACNVGGGTWTCTAAT 32 /home/toju/miniconda3/bin/cutadapt

sh /home/toju/Desktop/Scripts/noguchi/Filtering2Annotation_support_single/03_readQC_dada2.sh 200 10 0.1 5 2 240 15 32
sh /home/toju/Desktop/Scripts/noguchi/Filtering2Annotation_support_single/04_DenoisingRemoveCHIMERA_dada2.sh 32 0.99 /home/toju/miniconda3/bin/vsearch
sh /home/toju/Desktop/Scripts/noguchi/Filtering2Annotation_support_single/05_TaxonomyAnnotation_dada2.sh /home/toju/Desktop/Scripts/noguchi/referenceDB/silva_nr99_v138_train_set_append_stdDNA.fa 32 0.05 19 5,90%
sh /home/toju/Desktop/Scripts/noguchi/Filtering2Annotation_support_single/DRA_helper_script01_rename_fastq.sh
sh /home/toju/Desktop/Scripts/noguchi/Filtering2Annotation_support_single/06_QualityCheck.sh 32

## ============================================= ##

