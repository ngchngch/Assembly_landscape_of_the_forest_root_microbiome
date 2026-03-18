#!bin/bash
que=APC
ncpu=8
ncpu_each=7
nmem=24
data=251213


mainDIR="/aptmp/mikihito/02_ELA_v060_241113_2"
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

TASK="02_01_03_ELA_prep_abundance_threshold_genuslevel"
QSUB="QSUBs/$TASK" 
data=251231
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
export OMP_NUM_THREADS=1

echo "setwd(\"${mainDIR}\")" > "scl_Rscripts/${TASK}$data.R"
cat "/aptmp/mikihito/02_ELA_v060_241113_2/Script/${TASK}.R" >> "scl_Rscripts/${TASK}$data.R"

Rscript --vanilla "scl_Rscripts/${TASK}$data.R"

EOF

qsub ${QSUB}".qsub"

################################################
###02_02###60

TASK="02_02_03_correspond_Comm_taxaOcc_Fungi_genuslevel"
QSUB="QSUBs/$TASK" 
data=251227
LOGPATH="LOGs_${TASK}_${data}"
mkdir -p LOGs2/$LOGPATH
mkdir -p Rscripts/${TASK}

for f in `seq 1 86`
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
export OMP_NUM_THREADS=1

echo "setwd(\"${mainDIR}\")" > "scl_Rscripts/${TASK}$data$f.R"
cat "/aptmp/mikihito/02_ELA_v060_241113_2/Script/${TASK}.R" >> "scl_Rscripts/${TASK}$data$f.R"

sed -i "12s/.*/n.core=${ncpu}/" "scl_Rscripts/${TASK}$data$f.R"
sed -i "13s/.*/ELA_prep_dir='02_01_03_ELA_prep_abundance_threshold_genuslevel'/" "scl_Rscripts/${TASK}$data$f.R"
sed -i "14s/.*/tx=$f/" "scl_Rscripts/${TASK}$data$f.R"

Rscript --vanilla "scl_Rscripts/${TASK}$data$f.R"
EOF

qsub ${QSUB}$f".qsub"

done


################################################
###02_02###175

TASK="02_02_03_correspond_Comm_taxaOcc_Prokarote_genuslevel"
QSUB="QSUBs/$TASK" 
data=260102
LOGPATH="LOGs_${TASK}_${data}"
mkdir -p LOGs2/$LOGPATH
mkdir -p Rscripts/${TASK}

for f in `seq 1 209`
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
export OMP_NUM_THREADS=1

echo "setwd(\"${mainDIR}\")" > "scl_Rscripts/${TASK}$data$f.R"
cat "/aptmp/mikihito/02_ELA_v060_241113_2/Script/${TASK}.R" >> "scl_Rscripts/${TASK}$data$f.R"

sed -i "12s/.*/n.core=${ncpu}/" "scl_Rscripts/${TASK}$data$f.R"
sed -i "13s/.*/ELA_prep_dir='02_01_03_ELA_prep_abundance_threshold_genuslevel'/" "scl_Rscripts/${TASK}$data$f.R"
sed -i "14s/.*/tx=$f/" "scl_Rscripts/${TASK}$data$f.R"

Rscript --vanilla "scl_Rscripts/${TASK}$data$f.R"
EOF

qsub ${QSUB}$f".qsub"

done



################################################
###02_03

TASK="02_03_03_summarize_basedR2_genuslevel"
QSUB="QSUBs/$TASK" 
data=250104
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
export OMP_NUM_THREADS=1

echo "setwd(\"${mainDIR}\")" > "scl_Rscripts/${TASK}$data.R"
cat "/aptmp/mikihito/02_ELA_v060_241113_2/Script/${TASK}.R" >> "scl_Rscripts/${TASK}$data.R"

sed -i "18s/.*/n.core=${ncpu}/" "scl_Rscripts/${TASK}_$data.R"
sed -i "19s/.*/ELA_prep_dir='02_01_03_ELA_prep_abundance_threshold_genuslevel'/" "scl_Rscripts/${TASK}$data.R"
sed -i "20s/.*/dir_02_02_03_F='02_02_03_correspond_Comm_taxaOcc_Fungi_genuslevel'/" "scl_Rscripts/${TASK}$data.R"
sed -i "21s/.*/dir_02_02_03_P='02_02_03_correspond_Comm_taxaOcc_Prokarote_genuslevel'/" "scl_Rscripts/${TASK}$data.R"

Rscript --vanilla "scl_Rscripts/${TASK}$data.R"

EOF

qsub ${QSUB}".qsub"


################################################
###02_04

TASK="02_04_03_taxa_select_basedR2_genuslevel"
QSUB="QSUBs/$TASK" 
data=260104
LOGPATH="LOGs_${TASK}_${data}"
mkdir -p LOGs2/$LOGPATH
mkdir -p Rscripts/${TASK}

for f in `seq 20 50`
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
export OMP_NUM_THREADS=1

echo "setwd(\"${mainDIR}\")" > "scl_Rscripts/${TASK}$data$f.R"
cat "/aptmp/mikihito/02_ELA_v060_241113_2/Script/${TASK}.R" >> "scl_Rscripts/${TASK}$data$f.R"

sed -i "18s/.*/dir_02_03_03='02_03_03_summarize_basedR2_genuslevel'/" "scl_Rscripts/${TASK}$data$f.R"
sed -i "19s/.*/ELA_prep_dir='02_01_03_ELA_prep_abundance_threshold_genuslevel'/" "scl_Rscripts/${TASK}$data$f.R"
sed -i "20s/.*/nSp=$f/" "scl_Rscripts/${TASK}$data$f.R"
sed -i "21s/.*/n.core=1/" "scl_Rscripts/${TASK}$data$f.R"

Rscript --vanilla "scl_Rscripts/${TASK}$data$f.R"
EOF

qsub ${QSUB}$f".qsub"

done


################################################
###02_05

TASK="02_05_03_summary_taxa_select_genuslevel_250108"
QSUB="QSUBs/$TASK" 
data=260108
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
export OMP_NUM_THREADS=1

echo "setwd(\"${mainDIR}\")" > "scl_Rscripts/${TASK}$data.R"
cat "/aptmp/mikihito/02_ELA_v060_241113_2/Script/${TASK}.R" >> "scl_Rscripts/${TASK}$data.R"

sed -i "15s/.*/n.core=${ncpu}/" "scl_Rscripts/${TASK}_$data.R"
sed -i "16s/.*/ELA_prep_dir='02_01_03_ELA_prep_abundance_threshold_genuslevel'/" "scl_Rscripts/${TASK}$data.R"
sed -i "17s/.*/dir_02_03_03='02_03_03_summarize_basedR2_genuslevel'/" "scl_Rscripts/${TASK}$data.R"
sed -i "18s/.*/dir_02_04_03='02_04_03_taxa_select_basedR2_genuslevel'/" "scl_Rscripts/${TASK}$data.R"

Rscript --vanilla "scl_Rscripts/${TASK}$data.R"

EOF

qsub ${QSUB}".qsub"


################################################
###02_06

TASK="02_06_07_ELA_genuslevel"
QSUB="QSUBs/$TASK" 
data=260114
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
export OMP_NUM_THREADS=1

echo "setwd(\"${mainDIR}\")" > "scl_Rscripts/${TASK}$data.R"
cat "/aptmp/mikihito/02_ELA_v060_241113_2/Script/${TASK}.R" >> "scl_Rscripts/${TASK}$data.R"

sed -i "39s/.*/n.core=${ncpu}/" "scl_Rscripts/${TASK}$data.R"
sed -i "40s/.*/ELA_prep_dir='02_01_03_ELA_prep_abundance_threshold_genuslevel'/" "scl_Rscripts/${TASK}$data.R"
sed -i "41s/.*/dir_02_05_03='02_05_03_summary_taxa_select_genuslevel'/" "scl_Rscripts/${TASK}$data.R"

Rscript --vanilla "scl_Rscripts/${TASK}$data.R"

EOF

qsub ${QSUB}".qsub"

###############################################
###03_01

TASK="03_01_03_ELA_withRA_4step_genuslevel"
QSUB="QSUBs/$TASK" 
data=260109
LOGPATH="LOGs_${TASK}_${data}"
mkdir -p LOGs2/$LOGPATH


for f in `seq 1 177`
do

cat << EOF > ${QSUB}$f".qsub"
#!/bin/bash
#PBS -q ${que}
#PBS -N $TASK$f
#PBS -l select=1:ncpus=${ncpu_each}:mem=${nmem}gb
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
export OMP_NUM_THREADS=1

echo "setwd(\"${mainDIR}\")" > "scl_Rscripts/${TASK}$data$f.R"
cat "/aptmp/mikihito/02_ELA_v060_241113_2/Script/${TASK}.R" >> "scl_Rscripts/${TASK}$data$f.R"
sed -i "40s/.*/ELA_prep_dir='02_01_03_ELA_prep_abundance_threshold_genuslevel'/" "scl_Rscripts/${TASK}$data$f.R"
sed -i "41s/.*/n.core=${ncpu}/" "scl_Rscripts/${TASK}$data$f.R"
sed -i "42s/.*/dir_02_05_03='02_05_03_summary_taxa_select_genuslevel_250108'/" "scl_Rscripts/${TASK}$data$f.R"
sed -i "43s/.*/sp=$f/" "scl_Rscripts/${TASK}$data$f.R"


Rscript --vanilla "scl_Rscripts/${TASK}$data$f.R"
EOF

qsub ${QSUB}$f".qsub"

done

###############################################
###03_03

TASK="03_02_04_randELA_withRA_4s_fixP_genuslevel"
QSUB="QSUBs/$TASK" 
data=260109
ncpu=1
LOGPATH="LOGs_${TASK}_${data}"
mkdir -p LOGs2/$LOGPATH
mkdir -p Rscripts/${TASK}


### --1-1000
for f in `seq 1 3000`
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
export OMP_NUM_THREADS=1


echo "setwd(\"${mainDIR}\")" > "scl_Rscripts/${TASK}$data$f.R"
cat "${mainDIR}/Script/${TASK}.R" >> "scl_Rscripts/${TASK}$data$f.R"
sed -i "40s/.*/ELA_prep_dir='02_01_03_ELA_prep_abundance_threshold_genuslevel'/" "scl_Rscripts/${TASK}$data$f.R"
sed -i "41s/.*/n.core=${ncpu}/" "scl_Rscripts/${TASK}$data$f.R"
sed -i "42s/.*/dir_02_05_03='02_05_03_summary_taxa_select_genuslevel_250108'/" "scl_Rscripts/${TASK}$data$f.R"
sed -i "43s/.*/rd=$f/" "scl_Rscripts/${TASK}$data$f.R"

Rscript --vanilla "scl_Rscripts/${TASK}$data$f.R"
EOF

qsub ${QSUB}$f".qsub"

done


###############################################
###03_03

TASK="03_02_04_02_Zconv_ELA_withRA_4step_genuslevel"
QSUB="QSUBs/$TASK" 
data=260114
LOGPATH="LOGs_${TASK}_${data}"
mkdir -p LOGs2/$LOGPATH
mkdir -p Rscripts/${TASK}


### --1-1000
for f in `seq 1 177`
do

cat << EOF > ${QSUB}$f".qsub"
#!/bin/bash
#PBS -q ${que}
#PBS -N $TASK$f
#PBS -l select=1:ncpus=${ncpu_each}:mem=${nmem}gb
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
export OMP_NUM_THREADS=1


echo "setwd(\"${mainDIR}\")" > "scl_Rscripts/${TASK}$data$f.R"
cat "${mainDIR}/Script/${TASK}.R" >> "scl_Rscripts/${TASK}$data$f.R"
sed -i "20s/.*/ELA_prep_dir='02_01_03_ELA_prep_abundance_threshold_genuslevel'/" "scl_Rscripts/${TASK}$data$f.R"
sed -i "21s/.*/n.core=${ncpu}/" "scl_Rscripts/${TASK}$data$f.R"
sed -i "22s/.*/dir_03_01_03='03_01_03_ELA_withRA_4step_genuslevel'/" "scl_Rscripts/${TASK}$data$f.R"
sed -i "23s/.*/dir_03_02_04='03_02_04_randELA_withRA_4s_fixP_genuslevel'/" "scl_Rscripts/${TASK}$data$f.R"
sed -i "24s/.*/sp=$f/" "scl_Rscripts/${TASK}$data$f.R"

Rscript --vanilla "scl_Rscripts/${TASK}$data$f.R"
EOF

qsub ${QSUB}$f".qsub"

done

################################################

TASK="03_02_04_03_graphics_Zconv_ELA_withRA_4step_genuslevel"
QSUB="QSUBs/$TASK" 
data=260114
LOGPATH="LOGs_${TASK}_${data}"
mkdir -p LOGs2/$LOGPATH
mkdir -p Rscripts/${TASK}

cat << EOF > ${QSUB}".qsub"
#!/bin/bash
#PBS -q ${que}
#PBS -N $TASK
#PBS -l select=1:ncpus=${ncpu_each}:mem=${nmem}gb
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
export OMP_NUM_THREADS=1

echo "setwd(\"${mainDIR}\")" > "scl_Rscripts/${TASK}$data.R"
cat "/aptmp/mikihito/02_ELA_v060_241113_2/Script/${TASK}.R" >> "scl_Rscripts/${TASK}$data.R"

sed -i "16s/.*/n.core=${ncpu}/" "scl_Rscripts/${TASK}$data.R"
sed -i "17s/.*/dir_03_02_04_02='03_02_04_02_Zconv_ELA_withRA_4step_genuslevel'/" "scl_Rscripts/${TASK}$data.R"

Rscript --vanilla "scl_Rscripts/${TASK}$data.R"

EOF

qsub ${QSUB}".qsub"

