#!/bin/bash
set -e

#####################################################
# 2022. 02. 28.
# 
# Read me ->
# This script is for setting the parameter.
# When you run this script, you set current directory at 'Sample_lane1" and
# put 'the index list file' and 'primer list file' into current directory.
#
####################################################
### -- Setting parameter -- ###

## ============================================== ##
# -- Input file names
expname=No410
index1=index1.txt
index2=index2.txt

prameter=parameterList.tsv
## ============================================== ##
# -- Set global options
thread=32

## ============================================== ##
# -- Select run
pairedEND=FALSE
skipDemultiplex=FALSE # Whether skip the demultiplex process. If skip, you input [TRUE]
skipRun=FALSE # Whether skip from fltering process to annotation. If skip, you input [TRUE]
skipCutadaptor=FALSE
skipFilt=FALSE # Whether skip fltering process. If skip, you input [TRUE]
skipClus=FALSE # Whether skip clustering process. If skip, you input [TRUE]
skipAssign=FALSE # Whether skip annotation process. If skip, you input [TRUE]
filename=FALSE # You can decide the file name of last ouput files. 

#####################################################
# -- Path of script and reference
if [ $pairedEND = TRUE ]; then
scriptpath='/home/toju/Desktop/Scripts/noguchi/Filtering2Annotation_support_pair'
else
scriptpath='/home/toju/Desktop/Scripts/noguchi/Filtering2Annotation_support_single'
fi

referencepath='/home/toju/Desktop/Scripts/noguchi/referenceDB'
cutadaptpath=`which cutadapt` # check which cutadapt
vsearchpath=`which vsearch` # check which cutadapt
demultidir=demultiplex
######################################################
## |||||||||||||||||||||||||||||||||||||||||||||||| ##
######################################################
## -- Make save directory and remove old analyis

## -- Log diirectory
if [ -d log_and_script ]; then
  rm -r log_and_script
fi
mkdir -p log_and_script

## -- move to new directory

## -- move to new directory
for i in `cut -f 1 ${prameter}`
do

	if ls ${i}-*_analysis >/dev/null 2>&1
	then
		timestamp=`date '+%F_%R'`

		mkdir -p old_${timestamp}

		mv ${i}-*_analysis old_${timestamp}

	fi
	
done

###########################################################
## -- Make Primer files
Rscript ${scriptpath}/makePrimerFile.R ${prameter}

###########################################################
## -- Run the demuliplex script
if [ $skipDemultiplex = FALSE ]; then
		
		## -- Make script
		sh $scriptpath/01_Demultiplex.sh $expname $index1 $index2 \
									  	 Fprimer.txt Rprimer.txt \
									  	 $demultidir ${thread}		
									  	 
fi

###########################################################
## === Splitting each primer result

while true
do
	read -r f1 <&3 || break
	read -r f2 <&4 || break

	if [ `echo "$f1" | grep ">"` ] ; then
  		
		i=` echo "$f1" | awk -F ">" '{ print $2 }' `
    	l=` echo "$f2" | awk -F ">" '{ print $2 }' `
	
  		## ############################################### ##
	   	if ls $demultidir/*${i}* >/dev/null 2>&1
	   	then    	
  			## =========================================== ##
  			## -- Take specific primer
	      	echo "The reads contained ${i}-${l} primer were detected."
	
			## =========================================== ##
			## -- Move to analysis directory
			mkdir -p ${i}-${l}_analysis/01_Demultiplexed_fastaq
			cp -r log_and_script ${i}-${l}_analysis/
		 
			cp $demultidir/*$i* ${i}-${l}_analysis/01_Demultiplexed_fastaq 
				    
			echo "The reads were moved to new directory"

			cd ${i}-${l}_analysis    
			## =========================================== ##
			## -- Make script files
			Rscript ${scriptpath}/makeScriptFile.R ../${prameter} ${cutadaptpath} ${vsearchpath} ${i} \
													${thread} ${referencepath} ${scriptpath}
			sh log_and_script/makeScript.sh
			## =========================================== ##
			## -- Run the scripts
			
			echo "" > log_and_script/RUNonceALL.sh
						
			cat << EOF >> log_and_script/RUNonceALL.sh
			## ++++++++++++++++++++++++++++++++++++++ ##	
			## -- Summary directory
			filename=$filename
			i=${i}
			l=${l}
EOF
			cat << 'EOF' >> log_and_script/RUNonceALL.sh			
			filename=FALSE
			if [ $filename = FALSE ]; then
				filename=""
				output=${i}-${l}_summary
			else
				filename=${filename}_
				output=${i}-${l}_${filename}summary
			fi
			mkdir -p $output
EOF

			echo "sh log_and_script/script02_Cutadaptor.sh  2>&1 | tee log_and_script/log02_Cutadaptor.txt" >> log_and_script/RUNonceALL.sh
			echo "sh log_and_script/script03_FilterTrimming.sh 2>&1 | tee log_and_script/log03_FilterTrimming.txt" >> log_and_script/RUNonceALL.sh
			
			if [ $skipClus = FALSE ]; then
				echo "sh log_and_script/script04_Denoising.sh  2>&1 | tee log_and_script/log04_Denoising.txt" >> log_and_script/RUNonceALL.sh
				echo 'cp 04_Denoising/*fasta 04_Denoising/seqtab_rmChimera.rds 04_Denoising/seqOTUtab.rds $output' >> log_and_script/RUNonceALL.sh
			fi
			
			if [ $skipAssign = FALSE ]; then
				echo "sh log_and_script/script05_TaxonomyAnnotation.sh 2>&1 | tee log_and_script/log05_TaxonomyAnnotation.txt" >> log_and_script/RUNonceALL.sh
				echo 'cp 05_TaxonomyAnnotation/taxonomy_list.txt $output/${filename}taxonomyList.txt' >> log_and_script/RUNonceALL.sh
				echo 'cp 05_TaxonomyAnnotation/taxonomy_list.rds $output/${filename}taxonomyList.rds' >> log_and_script/RUNonceALL.sh
			fi		
			
			## ========================================== ##

			if [ $skipRun = FALSE ]; then				
			    sh log_and_script/RUNonceALL.sh
			fi
			## =========================================== ##

			cd ../
		else

      			echo "The reads contained ${i}-${l} primer were NOT detected."

		fi
		## ############################################### ##
	fi

done 3<Fprimer.txt 4<Rprimer.txt

