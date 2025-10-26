#!bin/bash
que=APC
ncpu=6
nmem=24
data=241113

mainDIR="XXXXXXXX"
LOGPATH="LOGs2"
MODULE="R"

fb="all"
## -----------
if [ ${que} = "SMALL" ]; then
    walltime="#PBS -l walltime=12:00:00"
else
    walltime=""
fi
## -----------
cd $mainDIR
mkdir -p scl_Rscripts
mkdir -p QSUBs
mkdir -p LOGs2
mkdir -p Rscripts

################################################
###02_01

TASK="02_01_ELA_prep_abundance_threshold_250906"
QSUB="QSUBs/$TASK" 
data=250226
LOGPATH="LOGs_${TASK}_${data}"
mkdir -p LOGs2/$LOGPATH


cat << EOF > ${QSUB}".qsub"
#!/bin/bash
#PBS -q ${que}
#PBS -N $TASK_$TASK
#PBS -l select=1:ncpus=${ncpu}:mem=${nmem}gb
#PBS -e $mainDIR/$LOGPATH/$TASK.log
#PBS -j eo
$walltime

## define variables        
mainDIR=${mainDIR}
cd ${mainDIR}

ncpu=$ncpu; nmem=$nmem
QSUB=$QSUB
data=$data

## load app
source /etc/profile.d/modules.sh
module load $MODULE

echo "setwd(\"${mainDIR}\")" > "scl_Rscripts/${TASK}_$data.R"
cat "XXX/Script/${TASK}.R" >> "scl_Rscripts/${TASK}_$data.R"

sed -i "38s/.*/n.core=${ncpu}/" "scl_Rscripts/${TASK}_$data.R"

Rscript --vanilla "scl_Rscripts/${TASK}_$data.R"
EOF

qsub ${QSUB}".qsub"


################################################
###02_02

TASK="02_02_ELA_SSaccumulation_curve_test"
QSUB="QSUBs/$TASK" 


for f in 20 30 40 45 50 55 60
do

cat << EOF > ${QSUB}$f".qsub"
#!/bin/bash
#PBS -q ${que}
#PBS -N $TASK$f
#PBS -l select=1:ncpus=${ncpu}:mem=${nmem}gb
#PBS -e $mainDIR/$LOGPATH/$TASK$f.log
#PBS -j eo
$walltime

## define variables        
mainDIR=${mainDIR}
cd ${mainDIR}

ncpu=$ncpu; nmem=$nmem
QSUB=$QSUB
data=$data
f=$f

## load app
source /etc/profile.d/modules.sh
module load $MODULE

echo "setwd(\"${mainDIR}\")" > "scl_Rscripts/${TASK}$data$f.R"
cat "XXX/Script/${TASK}.R" >> "scl_Rscripts/${TASK}$data$f.R"

sed -i "22s/.*/n.core=${ncpu}/" "scl_Rscripts/${TASK}$data$f.R"
sed -i "23s/.*/ELA_prep_dir='02_01_ELA_prep_abundance_threshold'/" "scl_Rscripts/${TASK}$data$f.R"
sed -i "24s/.*/nSp=$f/" "scl_Rscripts/${TASK}$data$f.R"

Rscript --vanilla "scl_Rscripts/${TASK}$data$f.R"
EOF

qsub ${QSUB}$f".qsub"

done


################################################
###02_02###67

TASK="02_02_correspond_Comm_taxaOcc_Fungi"
QSUB="QSUBs/$TASK" 


for f in `seq 3 67`
do

cat << EOF > ${QSUB}$f".qsub"
#!/bin/bash
#PBS -q ${que}
#PBS -N $TASK$f
#PBS -l select=1:ncpus=${ncpu}:mem=${nmem}gb
#PBS -e $mainDIR/$LOGPATH/$TASK$f.log
#PBS -j eo
$walltime

## define variables        
mainDIR=${mainDIR}
cd ${mainDIR}

ncpu=$ncpu; nmem=$nmem
QSUB=$QSUB
data=$data
f=$f

## load app
source /etc/profile.d/modules.sh
module load $MODULE

echo "setwd(\"${mainDIR}\")" > "scl_Rscripts/${TASK}$data$f.R"
cat "XXX/Script/${TASK}.R" >> "scl_Rscripts/${TASK}$data$f.R"

sed -i "12s/.*/n.core=${ncpu}/" "scl_Rscripts/${TASK}$data$f.R"
sed -i "13s/.*/ELA_prep_dir='02_01_ELA_prep_abundance_threshold'/" "scl_Rscripts/${TASK}$data$f.R"
sed -i "14s/.*/tx=$f/" "scl_Rscripts/${TASK}$data$f.R"

Rscript --vanilla "scl_Rscripts/${TASK}$data$f.R"
EOF

qsub ${QSUB}$f".qsub"

done


################################################
###02_02###175

TASK="02_02_correspond_Comm_taxaOcc_Prokarote"
QSUB="QSUBs/$TASK" 


for f in `seq 3 175`
do

cat << EOF > ${QSUB}$f".qsub"
#!/bin/bash
#PBS -q ${que}
#PBS -N $TASK$f
#PBS -l select=1:ncpus=${ncpu}:mem=${nmem}gb
#PBS -e $mainDIR/$LOGPATH/$TASK$f.log
#PBS -j eo
$walltime

## define variables        
mainDIR=${mainDIR}
cd ${mainDIR}

ncpu=$ncpu; nmem=$nmem
QSUB=$QSUB
data=$data
f=$f

## load app
source /etc/profile.d/modules.sh
module load $MODULE

echo "setwd(\"${mainDIR}\")" > "scl_Rscripts/${TASK}$data$f.R"
cat "XXX/Script/${TASK}.R" >> "scl_Rscripts/${TASK}$data$f.R"

sed -i "12s/.*/n.core=${ncpu}/" "scl_Rscripts/${TASK}$data$f.R"
sed -i "13s/.*/ELA_prep_dir='02_01_ELA_prep_abundance_threshold'/" "scl_Rscripts/${TASK}$data$f.R"
sed -i "14s/.*/tx=$f/" "scl_Rscripts/${TASK}$data$f.R"

Rscript --vanilla "scl_Rscripts/${TASK}$data$f.R"
EOF

qsub ${QSUB}$f".qsub"

done



################################################
###02_03

TASK="02_03_summarize_basedR2"
QSUB="QSUBs/$TASK" 


cat << EOF > ${QSUB}".qsub"
#!/bin/bash
#PBS -q ${que}
#PBS -N $TASK
#PBS -l select=1:ncpus=${ncpu}:mem=${nmem}gb
#PBS -e $mainDIR/$LOGPATH/$TASK.log
#PBS -j eo
$walltime

## define variables        
mainDIR=${mainDIR}
cd ${mainDIR}

ncpu=$ncpu; nmem=$nmem
QSUB=$QSUB
data=$data

## load app
source /etc/profile.d/modules.sh
module load $MODULE

echo "setwd(\"${mainDIR}\")" > "scl_Rscripts/${TASK}$data.R"
cat "XXX/Script/${TASK}.R" >> "scl_Rscripts/${TASK}$data.R"

sed -i "18s/.*/n.core=${ncpu}/" "scl_Rscripts/${TASK}_$data.R"
sed -i "19s/.*/ELA_prep_dir='02_01_ELA_prep_abundance_threshold'/" "scl_Rscripts/${TASK}$data.R"
sed -i "20s/.*/dir_02_02_F='02_02_correspond_Comm_taxaOcc_Fungi'/" "scl_Rscripts/${TASK}$data.R"
sed -i "21s/.*/dir_02_02_P='02_02_correspond_Comm_taxaOcc_Prokarote'/" "scl_Rscripts/${TASK}$data.R"

Rscript --vanilla "scl_Rscripts/${TASK}$data.R"

EOF

qsub ${QSUB}".qsub"


################################################
###02_04

TASK="02_04_taxa_select_basedR2"
QSUB="QSUBs/$TASK" 


for f in `seq 23 50`
do

cat << EOF > ${QSUB}$f".qsub"
#!/bin/bash
#PBS -q ${que}
#PBS -N $TASK$f
#PBS -l select=1:ncpus=${ncpu}:mem=${nmem}gb
#PBS -e $mainDIR/$LOGPATH/$TASK$f.log
#PBS -j eo
$walltime

## define variables        
mainDIR=${mainDIR}
cd ${mainDIR}

ncpu=$ncpu; nmem=$nmem
QSUB=$QSUB
data=$data
f=$f

## load app
source /etc/profile.d/modules.sh
module load $MODULE

echo "setwd(\"${mainDIR}\")" > "scl_Rscripts/${TASK}$data$f.R"
cat "XXX/Script/${TASK}.R" >> "scl_Rscripts/${TASK}$data$f.R"

sed -i "18s/.*/dir_02_03='02_03_summarize_basedR2'/" "scl_Rscripts/${TASK}$data$f.R"
sed -i "19s/.*/ELA_prep_dir='02_01_ELA_prep_abundance_threshold'/" "scl_Rscripts/${TASK}$data$f.R"
sed -i "20s/.*/nSp=$f/" "scl_Rscripts/${TASK}$data$f.R"
sed -i "21s/.*/n.core=1/" "scl_Rscripts/${TASK}$data$f.R"

Rscript --vanilla "scl_Rscripts/${TASK}$data$f.R"
EOF

qsub ${QSUB}$f".qsub"

done

################################################
###02_04

TASK="02_04_taxa_select_basedR2"
QSUB="QSUBs/$TASK" 


cat << EOF > ${QSUB}".qsub"
#!/bin/bash
#PBS -q ${que}
#PBS -N $TASK
#PBS -l select=1:ncpus=${ncpu}:mem=${nmem}gb
#PBS -e $mainDIR/$LOGPATH/$TASK.log
#PBS -j eo
$walltime

## define variables        
mainDIR=${mainDIR}
cd ${mainDIR}

ncpu=$ncpu; nmem=$nmem
QSUB=$QSUB
data=$data

## load app
source /etc/profile.d/modules.sh
module load $MODULE

echo "setwd(\"${mainDIR}\")" > "scl_Rscripts/${TASK}$data.R"
cat "XXX/Script/${TASK}.R" >> "scl_Rscripts/${TASK}$data.R"

sed -i "18s/.*/n.core=${ncpu}/" "scl_Rscripts/${TASK}_$data.R"
sed -i "19s/.*/ELA_prep_dir='02_01_ELA_prep_abundance_threshold'/" "scl_Rscripts/${TASK}$data.R"
sed -i "20s/.*/dir_02_02_F='02_02_correspond_Comm_taxaOcc_Fungi'/" "scl_Rscripts/${TASK}$data.R"
sed -i "21s/.*/dir_02_02_P='02_02_correspond_Comm_taxaOcc_Prokarote'/" "scl_Rscripts/${TASK}$data.R"

Rscript --vanilla "scl_Rscripts/${TASK}$data.R"

EOF

qsub ${QSUB}".qsub"


################################################
###02_05

TASK="02_05_summary_taxa_select"
QSUB="QSUBs/$TASK" 


cat << EOF > ${QSUB}".qsub"
#!/bin/bash
#PBS -q ${que}
#PBS -N $TASK
#PBS -l select=1:ncpus=${ncpu}:mem=${nmem}gb
#PBS -e $mainDIR/$LOGPATH/$TASK.log
#PBS -j eo
$walltime

## define variables        
mainDIR=${mainDIR}
cd ${mainDIR}

ncpu=$ncpu; nmem=$nmem
QSUB=$QSUB
data=$data

## load app
source /etc/profile.d/modules.sh
module load $MODULE

echo "setwd(\"${mainDIR}\")" > "scl_Rscripts/${TASK}$data.R"
cat "XXX/Script/${TASK}.R" >> "scl_Rscripts/${TASK}$data.R"

sed -i "15s/.*/n.core=${ncpu}/" "scl_Rscripts/${TASK}_$data.R"
sed -i "16s/.*/ELA_prep_dir='02_01_ELA_prep_abundance_threshold'/" "scl_Rscripts/${TASK}$data.R"
sed -i "17s/.*/dir_02_03='02_03_summarize_basedR2'/" "scl_Rscripts/${TASK}$data.R"
sed -i "18s/.*/dir_02_04='02_04_taxa_select_basedR2'/" "scl_Rscripts/${TASK}$data.R"

Rscript --vanilla "scl_Rscripts/${TASK}$data.R"

EOF

qsub ${QSUB}".qsub"


################################################
###02_06

TASK="02_06_ELA"
QSUB="QSUBs/$TASK" 


cat << EOF > ${QSUB}".qsub"
#!/bin/bash
#PBS -q ${que}
#PBS -N $TASK
#PBS -l select=1:ncpus=${ncpu}:mem=${nmem}gb
#PBS -e $mainDIR/$LOGPATH/$TASK.log
#PBS -j eo
$walltime

## define variables        
mainDIR=${mainDIR}
cd ${mainDIR}

ncpu=$ncpu; nmem=$nmem
QSUB=$QSUB
data=$data

## load app
source /etc/profile.d/modules.sh
module load $MODULE

echo "setwd(\"${mainDIR}\")" > "scl_Rscripts/${TASK}$data.R"
cat "XXX/Script/${TASK}.R" >> "scl_Rscripts/${TASK}$data.R"

sed -i "39s/.*/n.core=${ncpu}/" "scl_Rscripts/${TASK}$data.R"
sed -i "40s/.*/ELA_prep_dir='02_01_ELA_prep_abundance_threshold'/" "scl_Rscripts/${TASK}$data.R"
sed -i "41s/.*/dir_02_05='02_05_summary_taxa_select'/" "scl_Rscripts/${TASK}$data.R"

Rscript --vanilla "scl_Rscripts/${TASK}$data.R"

EOF

qsub ${QSUB}".qsub"


################################################
###02_06

TASK="02_07_assemblygraph_onlyBasin"
QSUB="QSUBs/$TASK" 
data=250311
LOGPATH="LOGs_${TASK}_${data}"
mkdir -p LOGs2/$LOGPATH
mkdir -p Rscripts/${TASK}

cat << EOF > ${QSUB}".qsub"
#!/bin/bash
#PBS -q ${que}
#PBS -N $TASK
#PBS -l select=1:ncpus=${ncpu}:mem=${nmem}gb
#PBS -e $mainDIR/LOGs2/$LOGPATH/$TASK.log
#PBS -j eo
$walltime

## define variables        
mainDIR=${mainDIR}
cd ${mainDIR}

ncpu=$ncpu; nmem=$nmem
QSUB=$QSUB
data=$data

## load app
source /etc/profile.d/modules.sh
module load $MODULE

echo "setwd(\"${mainDIR}\")" > "Rscripts/${TASK}/${TASK}$data.R"
cat "XXX/Script/${TASK}.R" >> "Rscripts/${TASK}/${TASK}$data.R"

sed -i "183s/.*/n.core=${ncpu}/" "Rscripts/${TASK}/${TASK}$data.R"
sed -i "184s/.*/ELA_prep_dir='02_01_ELA_prep_abundance_threshold'/" "Rscripts/${TASK}/${TASK}$data.R"
sed -i "185s/.*/dir_02_05='02_05_summary_taxa_select'/" "Rscripts/${TASK}/${TASK}$data.R"
sed -i "186s/.*/dir_02_06='02_06_ELA'/" "Rscripts/${TASK}/${TASK}$data.R"

Rscript --vanilla "Rscripts/${TASK}/${TASK}$data.R"

EOF

qsub ${QSUB}".qsub"


################################################
###02_06

TASK="02_08_basin_barplot"
QSUB="QSUBs/$TASK" 


cat << EOF > ${QSUB}".qsub"
#!/bin/bash
#PBS -q ${que}
#PBS -N $TASK
#PBS -l select=1:ncpus=${ncpu}:mem=${nmem}gb
#PBS -e $mainDIR/$LOGPATH/$TASK.log
#PBS -j eo
$walltime

## define variables        
mainDIR=${mainDIR}
cd ${mainDIR}

ncpu=$ncpu; nmem=$nmem
QSUB=$QSUB
data=$data

## load app
source /etc/profile.d/modules.sh
module load $MODULE

echo "setwd(\"${mainDIR}\")" > "scl_Rscripts/${TASK}$data.R"
cat "XXX/Script/${TASK}.R" >> "scl_Rscripts/${TASK}$data.R"

sed -i "39s/.*/n.core=${ncpu}/" "scl_Rscripts/${TASK}$data.R"
sed -i "40s/.*/ELA_prep_dir='02_01_ELA_prep_abundance_threshold'/" "scl_Rscripts/${TASK}$data.R"
sed -i "41s/.*/dir_02_05='02_05_summary_taxa_select'/" "scl_Rscripts/${TASK}$data.R"
sed -i "42s/.*/dir_02_06='02_06_ELA'/" "scl_Rscripts/${TASK}$data.R"

Rscript --vanilla "scl_Rscripts/${TASK}$data.R"

EOF

qsub ${QSUB}".qsub"


################################################
###02_09

TASK="02_09_landscape_distance_eachPlant"
QSUB="QSUBs/$TASK" 


cat << EOF > ${QSUB}".qsub"
#!/bin/bash
#PBS -q ${que}
#PBS -N $TASK
#PBS -l select=1:ncpus=${ncpu}:mem=${nmem}gb
#PBS -e $mainDIR/$LOGPATH/$TASK.log
#PBS -j eo
$walltime

## define variables        
mainDIR=${mainDIR}
cd ${mainDIR}

ncpu=$ncpu; nmem=$nmem
QSUB=$QSUB
data=$data

## load app
source /etc/profile.d/modules.sh
module load $MODULE

echo "setwd(\"${mainDIR}\")" > "scl_Rscripts/${TASK}$data.R"
cat "XXX/Script/${TASK}.R" >> "scl_Rscripts/${TASK}$data.R"

sed -i "155s/.*/n.core=${ncpu}/" "scl_Rscripts/${TASK}$data.R"
sed -i "156s/.*/ELA_prep_dir='02_01_ELA_prep_abundance_threshold'/" "scl_Rscripts/${TASK}$data.R"
sed -i "158s/.*/dir_02_06='02_06_ELA'/" "scl_Rscripts/${TASK}$data.R"

Rscript --vanilla "scl_Rscripts/${TASK}$data.R"

EOF

qsub ${QSUB}".qsub"


###############################################
###03_01

TASK="03_01_ELA_withRA_4step"
QSUB="QSUBs/$TASK" 
data=250307
LOGPATH="LOGs_${TASK}_${data}"
mkdir -p LOGs2/$LOGPATH



for f in 3 8 12 14 15 16 21 22 24 81 83 84 111 121 134
do

cat << EOF > ${QSUB}$f".qsub"
#!/bin/bash
#PBS -q ${que}
#PBS -N $TASK$f
#PBS -l select=1:ncpus=${ncpu}:mem=${nmem}gb
#PBS -e $mainDIR/LOGs2/$LOGPATH/$TASK$f.log
#PBS -j eo
$walltime

## define variables        
mainDIR=${mainDIR}
cd ${mainDIR}

ncpu=$ncpu; nmem=$nmem
QSUB=$QSUB
data=$data
f=$f

## load app
source /etc/profile.d/modules.sh
module load $MODULE

echo "setwd(\"${mainDIR}\")" > "scl_Rscripts/${TASK}$data$f.R"
cat "XXX/Script/${TASK}.R" >> "scl_Rscripts/${TASK}$data$f.R"
sed -i "40s/.*/ELA_prep_dir='02_01_ELA_prep_abundance_threshold'/" "scl_Rscripts/${TASK}$data$f.R"
sed -i "41s/.*/n.core=${ncpu}/" "scl_Rscripts/${TASK}$data$f.R"
sed -i "42s/.*/dir_02_05='02_05_summary_taxa_select'/" "scl_Rscripts/${TASK}$data$f.R"
sed -i "43s/.*/sp=$f/" "scl_Rscripts/${TASK}$data$f.R"


Rscript --vanilla "scl_Rscripts/${TASK}$data$f.R"
EOF

qsub ${QSUB}$f".qsub"

done

#######################################
###############################################
###03_02_02 summarize

TASK="03_02_summarize_randELA_withRA_fixP_250306"
QSUB="QSUBs/$TASK" 
data=250306
LOGPATH="LOGs_${TASK}_${data}"
mkdir -p LOGs2/$LOGPATH



for f in `seq 1 176`
do

cat << EOF > ${QSUB}$f".qsub"
#!/bin/bash
#PBS -q ${que}
#PBS -N $TASK$f
#PBS -l select=1:ncpus=${ncpu}:mem=${nmem}gb
#PBS -e $mainDIR/LOGs2/$LOGPATH/$TASK$f.log
#PBS -j eo
$walltime

## define variables        
mainDIR=${mainDIR}
cd ${mainDIR}

ncpu=$ncpu; nmem=$nmem
QSUB=$QSUB
data=$data
f=$f

## load app
source /etc/profile.d/modules.sh
module load $MODULE

echo "setwd(\"${mainDIR}\")" > "scl_Rscripts/${TASK}$data$f.R"
cat "XXX/Script/${TASK}.R" >> "scl_Rscripts/${TASK}$data$f.R"
sed -i "40s/.*/ELA_prep_dir='02_01_ELA_prep_abundance_threshold'/" "scl_Rscripts/${TASK}$data$f.R"
sed -i "41s/.*/n.core=${ncpu}/" "scl_Rscripts/${TASK}$data$f.R"
sed -i "42s/.*/sp=$f/" "scl_Rscripts/${TASK}$data$f.R"


Rscript --vanilla "scl_Rscripts/${TASK}$data$f.R"
EOF

qsub ${QSUB}$f".qsub"

done


###############################################
###03_03

TASK="03_03_ELA_withRA_4step_rep"
QSUB="QSUBs/$TASK"  
data=250306
LOGPATH="LOGs_${TASK}_${data}"
mkdir -p LOGs2/$LOGPATH



for f in `seq 1 176`
do

cat << EOF > ${QSUB}$f".qsub"
#!/bin/bash
#PBS -q ${que}
#PBS -N $TASK$f
#PBS -l select=1:ncpus=${ncpu}:mem=${nmem}gb
#PBS -e $mainDIR/LOGs2/$LOGPATH/$TASK$f.log
#PBS -j eo
$walltime

## define variables        
mainDIR=${mainDIR}
cd ${mainDIR}

ncpu=$ncpu; nmem=$nmem
QSUB=$QSUB
data=$data
f=$f

## load app
source /etc/profile.d/modules.sh
module load $MODULE

echo "setwd(\"${mainDIR}\")" > "scl_Rscripts/${TASK}$data$f.R"
cat "XXX/Script/${TASK}.R" >> "scl_Rscripts/${TASK}$data$f.R"
sed -i "40s/.*/ELA_prep_dir='02_01_ELA_prep_abundance_threshold'/" "scl_Rscripts/${TASK}$data$f.R"
sed -i "41s/.*/n.core=${ncpu}/" "scl_Rscripts/${TASK}$data$f.R"
sed -i "42s/.*/dir_02_05='02_05_summary_taxa_select'/" "scl_Rscripts/${TASK}$data$f.R"
sed -i "43s/.*/spi=$f/" "scl_Rscripts/${TASK}$data$f.R"


Rscript --vanilla "scl_Rscripts/${TASK}$data$f.R"
EOF

qsub ${QSUB}$f".qsub"

done


###############################################
###03_04

TASK="03_04_landscape_change_repuroducibility"
QSUB="QSUBs/$TASK" 
data=250306
LOGPATH="LOGs_${TASK}_${data}"
mkdir -p LOGs2/$LOGPATH




cat << EOF > ${QSUB}$f".qsub"
#!/bin/bash
#PBS -q ${que}
#PBS -N $TASK$f
#PBS -l select=1:ncpus=${ncpu}:mem=${nmem}gb
#PBS -e $mainDIR/LOGs2/$LOGPATH/$TASK$f.log
#PBS -j eo
$walltime

## define variables        
mainDIR=${mainDIR}
cd ${mainDIR}

ncpu=$ncpu; nmem=$nmem
QSUB=$QSUB
data=$data

## load app
source /etc/profile.d/modules.sh
module load $MODULE

echo "setwd(\"${mainDIR}\")" > "scl_Rscripts/${TASK}$data$f.R"
cat "XXX/Script/${TASK}.R" >> "scl_Rscripts/${TASK}$data$f.R"
sed -i "8s/.*/dir_03_03='03_03_ELA_withRA_4step_rep_241212'/" "scl_Rscripts/${TASK}$data$f.R"


Rscript --vanilla "scl_Rscripts/${TASK}$data$f.R"
EOF

qsub ${QSUB}$f".qsub"




###############################################
###03_03

TASK="03_05_states_flow_diagram"
QSUB="QSUBs/$TASK" 
data=250306
LOGPATH="LOGs_${TASK}_${data}"
mkdir -p LOGs2/$LOGPATH



for f in `seq 1 176`
do

cat << EOF > ${QSUB}$f".qsub"
#!/bin/bash
#PBS -q ${que}
#PBS -N $TASK$f
#PBS -l select=1:ncpus=${ncpu}:mem=${nmem}gb
#PBS -e $mainDIR/LOGs2/$LOGPATH/$TASK$f.log
#PBS -j eo
$walltime

## define variables        
mainDIR=${mainDIR}
cd ${mainDIR}

ncpu=$ncpu; nmem=$nmem
QSUB=$QSUB
data=$data
f=$f

## load app
source /etc/profile.d/modules.sh
module load $MODULE

echo "setwd(\"${mainDIR}\")" > "scl_Rscripts/${TASK}$data$f.R"
cat "XXX/Script/${TASK}.R" >> "scl_Rscripts/${TASK}$data$f.R"
sed -i "219s/.*/ELA_prep_dir='02_01_ELA_prep_abundance_threshold'/" "scl_Rscripts/${TASK}$data$f.R"
sed -i "220s/.*/n.core=${ncpu}/" "scl_Rscripts/${TASK}$data$f.R"
sed -i "221s/.*/dir_02_05='02_05_summary_taxa_select'/" "scl_Rscripts/${TASK}$data$f.R"
sed -i "222s/.*/dir_03_01='03_01_ELA_withRA_4step_250227'/" "scl_Rscripts/${TASK}$data$f.R"
sed -i "223s/.*/spi=$f/" "scl_Rscripts/${TASK}$data$f.R"
sed -i "224s/.*/dir_02_06='02_06_ELA'/" "scl_Rscripts/${TASK}$data$f.R"


Rscript --vanilla "scl_Rscripts/${TASK}$data$f.R"
EOF

qsub ${QSUB}$f".qsub"

done







###############################################
###03_03

TASK="03_05_states_flow_diagram_select_250312"
QSUB="QSUBs/$TASK" 
data=250311
LOGPATH="LOGs_${TASK}_${data}"
mkdir -p LOGs2/$LOGPATH


#42

for f in `seq 1 74`
do

cat << EOF > ${QSUB}$f".qsub"
#!/bin/bash
#PBS -q ${que}
#PBS -N $TASK$f
#PBS -l select=1:ncpus=${ncpu}:mem=${nmem}gb
#PBS -e $mainDIR/LOGs2/$LOGPATH/$TASK$f.log
#PBS -j eo
$walltime

## define variables        
mainDIR=${mainDIR}
cd ${mainDIR}

ncpu=$ncpu; nmem=$nmem
QSUB=$QSUB
data=$data
f=$f

## load app
source /etc/profile.d/modules.sh
module load $MODULE

echo "setwd(\"${mainDIR}\")" > "scl_Rscripts/${TASK}$data$f.R"
cat "XXX/Script/${TASK}.R" >> "scl_Rscripts/${TASK}$data$f.R"
sed -i "227s/.*/ELA_prep_dir='02_01_ELA_prep_abundance_threshold'/" "scl_Rscripts/${TASK}$data$f.R"
sed -i "228s/.*/n.core=${ncpu}/" "scl_Rscripts/${TASK}$data$f.R"
sed -i "229s/.*/dir_02_05='02_05_summary_taxa_select'/" "scl_Rscripts/${TASK}$data$f.R"
sed -i "230s/.*/dir_03_01='03_01_ELA_withRA_4step_250227'/" "scl_Rscripts/${TASK}$data$f.R"
sed -i "231s/.*/spi=$f/" "scl_Rscripts/${TASK}$data$f.R"
sed -i "232s/.*/dir_02_06='02_06_ELA'/" "scl_Rscripts/${TASK}$data$f.R"
sed -i "233s/.*/dir_02_07='02_07_assemblygraph_250311'/" "scl_Rscripts/${TASK}$data$f.R"


Rscript --vanilla "scl_Rscripts/${TASK}$data$f.R"
EOF

qsub ${QSUB}$f".qsub"

done
