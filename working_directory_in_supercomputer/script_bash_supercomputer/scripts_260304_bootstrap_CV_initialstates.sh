#!bin/bash
que=APC
ncpu=2
ncpu_each=2
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
###02_06_02

TASK="02_06_02_ELA_bootstrap_260304_2"
QSUB="QSUBs/$TASK" 
data=260304
LOGPATH="LOGs_${TASK}_${data}"
mkdir -p LOGs2/$LOGPATH

for f in `seq 1 1100`
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

## load app
source /etc/profile.d/modules.sh
module load $MODULE
export OMP_NUM_THREADS=1

echo "setwd(\"${mainDIR}\")" > "scl_Rscripts/${TASK}$data$f.R"
cat "/aptmp/mikihito/02_ELA_v060_241113_2/Script/${TASK}.R" >> "scl_Rscripts/${TASK}$data$f.R"

sed -i "32s/.*/n.core=${ncpu}/" "scl_Rscripts/${TASK}$data$f.R"
sed -i "33s/.*/ELA_prep_dir='02_01_ELA_prep_abundance_threshold'/" "scl_Rscripts/${TASK}$data$f.R"
sed -i "34s/.*/dir_02_05='02_05_summary_taxa_select'/" "scl_Rscripts/${TASK}$data$f.R"
sed -i "35s/.*/dir_02_06='02_06_ELA'/" "scl_Rscripts/${TASK}$data$f.R"
sed -i "36s/.*/boot_id=${f}/" "scl_Rscripts/${TASK}$data$f.R"

Rscript --vanilla "scl_Rscripts/${TASK}$data$f.R"

EOF

qsub ${QSUB}$f".qsub"
done

################################################
###02_06_03

TASK="02_06_03_02_summarize_bootELA_bcom_DG"
QSUB="QSUBs/$TASK" 
data=260304
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

sed -i "3s/.*/n.core=${ncpu}/" "scl_Rscripts/${TASK}$data.R"
sed -i "4s/.*/dir_02_06_02='02_06_02_ELA_bootstrap_260304_2'/" "scl_Rscripts/${TASK}$data.R"
sed -i "5s/.*/dir_02_06='02_06_ELA'/" "scl_Rscripts/${TASK}$data.R"

Rscript --vanilla "scl_Rscripts/${TASK}$data.R"

EOF

qsub ${QSUB}".qsub"


################################################
###02_06_03

TASK="02_06_03_01_summarize_ELA_bootstrap_enedif_260304_2"
QSUB="QSUBs/$TASK" 
data=260304
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

sed -i "24s/.*/n.core=${ncpu}/" "scl_Rscripts/${TASK}$data.R"
sed -i "25s/.*/dir_02_06_02='02_06_02_ELA_bootstrap_260304_2'/" "scl_Rscripts/${TASK}$data.R"
sed -i "26s/.*/dir_02_06='02_06_ELA'/" "scl_Rscripts/${TASK}$data.R"
sed -i "27s/.*/dir_03_10='03_10_graphics_states_flow_flow_Spl_250508'/" "scl_Rscripts/${TASK}$data.R"

Rscript --vanilla "scl_Rscripts/${TASK}$data.R"

EOF

qsub ${QSUB}".qsub"


################################################
###02_06_04

TASK="02_06_04_ELA_initial_state_numebr_260130"
QSUB="QSUBs/$TASK" 
data=260130
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

sed -i "32s/.*/n.core=${ncpu}/" "scl_Rscripts/${TASK}$data.R"
sed -i "33s/.*/ELA_prep_dir='02_01_ELA_prep_abundance_threshold'/" "scl_Rscripts/${TASK}$data.R"
sed -i "34s/.*/dir_02_05='02_05_summary_taxa_select'/" "scl_Rscripts/${TASK}$data.R"
sed -i "35s/.*/dir_02_06='02_06_ELA'/" "scl_Rscripts/${TASK}$data.R"

Rscript --vanilla "scl_Rscripts/${TASK}$data.R"

EOF

qsub ${QSUB}".qsub"

################################################
###02_06_05

TASK="02_06_05_graphics_summarize_ELA_bootstrap"
QSUB="QSUBs/$TASK" 
data=251215
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
echo "set.seed(1234)" > "scl_Rscripts/${TASK}$data.R"
cat "/aptmp/mikihito/02_ELA_v060_241113_2/Script/${TASK}.R" >> "scl_Rscripts/${TASK}$data.R"

sed -i "11s/.*/dir_02_06='02_06_ELA'/" "scl_Rscripts/${TASK}$data.R"
sed -i "12s/.*/dir_02_06_03='02_06_03_summarize_ELA_bootstrap'/" "scl_Rscripts/${TASK}$data.R"
sed -i "13s/.*/dir_03_10_02='03_10_02_SScomposition_recolor'/" "scl_Rscripts/${TASK}$data.R"

Rscript --vanilla "scl_Rscripts/${TASK}$data.R"

EOF

qsub ${QSUB}".qsub"

################################################
###02_06_04

TASK="02_06_08_ELA_100fold_CV"
QSUB="QSUBs/$TASK" 
data=260108
LOGPATH="LOGs_${TASK}_${data}"
mkdir -p LOGs2/$LOGPATH
mkdir -p Rscripts/${TASK}

for f in `seq 1 100`
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

## load app
source /etc/profile.d/modules.sh
module load $MODULE
export OMP_NUM_THREADS=1

echo "setwd(\"${mainDIR}\")" > "scl_Rscripts/${TASK}$data$f.R"
cat "/aptmp/mikihito/02_ELA_v060_241113_2/Script/${TASK}.R" >> "scl_Rscripts/${TASK}$data$f.R"

sed -i "26s/.*/n.core=${ncpu}/" "scl_Rscripts/${TASK}$data$f.R"
sed -i "27s/.*/ELA_prep_dir='02_01_ELA_prep_abundance_threshold'/" "scl_Rscripts/${TASK}$data$f.R"
sed -i "28s/.*/dir_02_05='02_05_summary_taxa_select'/" "scl_Rscripts/${TASK}$data$f.R"
sed -i "29s/.*/dir_02_06='02_06_ELA'/" "scl_Rscripts/${TASK}$data$f.R"
sed -i "30s/.*/nfold=${f}/" "scl_Rscripts/${TASK}$data$f.R"

Rscript --vanilla "scl_Rscripts/${TASK}$data$f.R"

EOF

qsub ${QSUB}$f".qsub"

done

################################################
###02_06_04

TASK="02_06_08_02_graphics_ELA_100fold_CV"
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
echo "set.seed(1234)" > "scl_Rscripts/${TASK}$data.R"
cat "/aptmp/mikihito/02_ELA_v060_241113_2/Script/${TASK}.R" >> "scl_Rscripts/${TASK}$data.R"

sed -i "10s/.*/dir_02_06_08='02_06_08_ELA_100fold_CV'/" "scl_Rscripts/${TASK}$data.R"

Rscript --vanilla "scl_Rscripts/${TASK}$data.R"

EOF

qsub ${QSUB}".qsub"
###############################################
###03_01

TASK="03_01_04_ELA_withRA_4step_boot"
QSUB="QSUBs/$TASK" 
data=260124
LOGPATH="LOGs_${TASK}_${data}"
mkdir -p LOGs2/$LOGPATH

for r in 47 48
do

cat << EOF > ${QSUB}$r".qsub"
#!/bin/bash
#PBS -q ${que}
#PBS -N $TASK$r
#PBS -l select=1:ncpus=${ncpu_each}:mem=${nmem}gb
#PBS -e $mainDIR/LOGs2/$LOGPATH/$TASK$r.log
#PBS -j eo
$walltime

## define variables        
mainDIR=${mainDIR}
cd ${mainDIR}

ncpu=$ncpu; nmem=$nmem
QSUB=$QSUB
data=$data
r=$r

## load app
source /etc/profile.d/modules.sh
module load $MODULE
export OMP_NUM_THREADS=1

echo "setwd(\"${mainDIR}\")" > "scl_Rscripts/${TASK}$data$r.R"
cat "/aptmp/mikihito/02_ELA_v060_241113_2/Script/${TASK}.R" >> "scl_Rscripts/${TASK}$data$r.R"
sed -i "40s/.*/ELA_prep_dir='02_01_ELA_prep_abundance_threshold'/" "scl_Rscripts/${TASK}$data$r.R"
sed -i "41s/.*/n.core=${ncpu}/" "scl_Rscripts/${TASK}$data$r.R"
sed -i "42s/.*/dir_02_05='02_05_summary_taxa_select'/" "scl_Rscripts/${TASK}$data$r.R"
sed -i "44s/.*/boot_id=$r/" "scl_Rscripts/${TASK}$data$r.R"


Rscript --vanilla "scl_Rscripts/${TASK}$data$r.R"
EOF

qsub ${QSUB}$r".qsub"

done

###############################################
###03_01

TASK="03_02_05_randELA_withRA_4step_boot"
QSUB="QSUBs/$TASK" 
data=260127
ncpu=1
LOGPATH="LOGs_${TASK}_${data}"
mkdir -p LOGs2/$LOGPATH

for r in `seq 1 10000`
do

cat << EOF > ${QSUB}$r".qsub"
#!/bin/bash
#PBS -q ${que}
#PBS -N $TASK$r
#PBS -l select=1:ncpus=${ncpu}:mem=${nmem}gb
#PBS -e $mainDIR/LOGs2/$LOGPATH/$TASK$r.log
#PBS -j eo
$walltime

## define variables        
mainDIR=${mainDIR}
cd ${mainDIR}

ncpu=$ncpu; nmem=$nmem
QSUB=$QSUB
data=$data
r=$r

## load app
source /etc/profile.d/modules.sh
module load $MODULE
export OMP_NUM_THREADS=1

echo "setwd(\"${mainDIR}\")" > "scl_Rscripts/${TASK}$data$r.R"
cat "/aptmp/mikihito/02_ELA_v060_241113_2/Script/${TASK}.R" >> "scl_Rscripts/${TASK}$data$r.R"
sed -i "40s/.*/ELA_prep_dir='02_01_ELA_prep_abundance_threshold'/" "scl_Rscripts/${TASK}$data$r.R"
sed -i "41s/.*/n.core=${ncpu}/" "scl_Rscripts/${TASK}$data$r.R"
sed -i "42s/.*/dir_02_05='02_05_summary_taxa_select'/" "scl_Rscripts/${TASK}$data$r.R"
sed -i "43s/.*/seed=$r/" "scl_Rscripts/${TASK}$data$r.R"
sed -i "44s/.*/dir_03_01_04='03_01_04_ELA_withRA_4step_boot'/" "scl_Rscripts/${TASK}$data$r.R"


Rscript --vanilla "scl_Rscripts/${TASK}$data$r.R"
EOF

qsub ${QSUB}$r".qsub"

done


###############################################
###03_01

TASK="03_02_05_02_summarize_randELA_withRA_4step_boot"
QSUB="QSUBs/$TASK" 
data=260130
LOGPATH="LOGs_${TASK}_${data}"
mkdir -p LOGs2/$LOGPATH

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
sed -i "9s/.*/dir_03_02_05='03_02_05_randELA_withRA_4step_boot'/" "scl_Rscripts/${TASK}$data.R"
sed -i "10s/.*/dir_03_01_04='03_01_04_ELA_withRA_4step_boot'/" "scl_Rscripts/${TASK}$data.R"


Rscript --vanilla "scl_Rscripts/${TASK}$data.R"
EOF

qsub ${QSUB}".qsub"


