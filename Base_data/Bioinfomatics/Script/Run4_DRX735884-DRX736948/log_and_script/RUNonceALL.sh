
			## ++++++++++++++++++++++++++++++++++++++ ##	
			## -- Summary directory
			filename=FALSE
			i=515f
			l=806rB
			filename=FALSE
			if [ $filename = FALSE ]; then
				filename=""
				output=${i}-${l}_summary
			else
				filename=${filename}_
				output=${i}-${l}_${filename}summary
			fi
			mkdir -p $output
sh log_and_script/script02_Cutadaptor.sh  2>&1 | tee log_and_script/log02_Cutadaptor.txt
sh log_and_script/script03_FilterTrimming.sh 2>&1 | tee log_and_script/log03_FilterTrimming.txt
sh log_and_script/script04_Denoising.sh  2>&1 | tee log_and_script/log04_Denoising.txt
cp 04_Denoising/*fasta 04_Denoising/seqtab_rmChimera.rds 04_Denoising/seqOTUtab.rds $output
sh log_and_script/script05_TaxonomyAnnotation.sh 2>&1 | tee log_and_script/log05_TaxonomyAnnotation.txt
cp 05_TaxonomyAnnotation/taxonomy_list.txt $output/${filename}taxonomyList.txt
cp 05_TaxonomyAnnotation/taxonomy_list.rds $output/${filename}taxonomyList.rds
