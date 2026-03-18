#!bin/bash
que=SMALL
ncpu_each=7
ncpu=8
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

TASK="03_02_00_prep_spatial_block"
QSUB="QSUBs/$TASK" 
data=251225
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

sed -i "5s/.*/ELA_prep_dir='02_01_ELA_prep_abundance_threshold'/" "scl_Rscripts/${TASK}$data$f.R"

Rscript --vanilla "scl_Rscripts/${TASK}$data.R"

EOF

qsub ${QSUB}".qsub"


###############################################
###03_03

TASK="03_02_02_randELA_withRA_4s_Spatial_block_withfixP_260104"
QSUB="QSUBs/$TASK" 
data=260111
ncpu=1
LOGPATH="LOGs_${TASK}_${data}"
mkdir -p LOGs2/$LOGPATH
mkdir -p Rscripts/${TASK}


### --1-1000
for f in `seq 301 3000`
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
sed -i "40s/.*/ELA_prep_dir='02_01_ELA_prep_abundance_threshold'/" "scl_Rscripts/${TASK}$data$f.R"
sed -i "41s/.*/n.core=${ncpu}/" "scl_Rscripts/${TASK}$data$f.R"
sed -i "42s/.*/dir_02_05='02_05_summary_taxa_select'/" "scl_Rscripts/${TASK}$data$f.R"
sed -i "43s/.*/dir_03_02_00='03_02_00_prep_spatial_block'/" "scl_Rscripts/${TASK}$data$f.R"
sed -i "44s/.*/rd=$f/" "scl_Rscripts/${TASK}$data$f.R"

Rscript --vanilla "scl_Rscripts/${TASK}$data$f.R"
EOF

qsub ${QSUB}$f".qsub"

done


###############################################
###03_03

TASK="03_02_02_02_Zconv_ELA_withRA_4step_spatialblocking"
QSUB="QSUBs/$TASK" 
data=260114
LOGPATH="LOGs_${TASK}_${data}"
mkdir -p LOGs2/$LOGPATH
mkdir -p Rscripts/${TASK}


### --1-1000
for f in `seq 1 176`
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
sed -i "20s/.*/ELA_prep_dir='02_01_ELA_prep_abundance_threshold'/" "scl_Rscripts/${TASK}$data$f.R"
sed -i "21s/.*/n.core=${ncpu}/" "scl_Rscripts/${TASK}$data$f.R"
sed -i "22s/.*/dir_03_01='03_01_ELA_withRA_4step_250227'/" "scl_Rscripts/${TASK}$data$f.R"
sed -i "23s/.*/dir_03_02_02='03_02_02_randELA_withRA_4s_Spatial_block_withfixP_260104'/" "scl_Rscripts/${TASK}$data$f.R"
sed -i "24s/.*/sp=$f/" "scl_Rscripts/${TASK}$data$f.R"

Rscript --vanilla "scl_Rscripts/${TASK}$data$f.R"
EOF

qsub ${QSUB}$f".qsub"

done

################################################

TASK="03_02_02_03_graphics_Zconv_ELA_withRA_4step_spatialblocking"
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
sed -i "17s/.*/dir_03_02_02_02='03_02_02_02_Zconv_ELA_withRA_4step_spatialblocking'/" "scl_Rscripts/${TASK}$data.R"

Rscript --vanilla "scl_Rscripts/${TASK}$data.R"

EOF

qsub ${QSUB}".qsub"
