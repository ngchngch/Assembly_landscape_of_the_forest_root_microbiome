#!bin/bash
que=APC
ncpu=8
ncpu_each=4
nmem=16
data=251213

mainDIR="/aptmp/hiroakif/noguchi_260514"
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
################################################
################################################
dir_05_01='05_01_01_gLV_simulation_Nsp80_nonoise_260514'
dir_05_01_00='05_01_01_00_summarize_gLV_simulation_Nsp80_nonoise_260514'
dir_05_02="05_02_01_ELA_withRA_multistep_eachseed_Nsp80_nonoise_260527"
dir_05_03="05_03_01_randELA_withRA_4step_Nsp80_nonoise_260527"
dir_05_04="05_04_01_summarize_ELA_withRA_4step_Nsp80_nonoise_260514"
dir_05_05="05_05_01_Zconv_ELA_withRA_4step_Nsp80_nonoise_260527"
###05_01

TASK=${dir_05_01}
QSUB="QSUBs/$TASK" 
data=260514
LOGPATH="LOGs_${TASK}_${data}"
mkdir -p LOGs2/$LOGPATH

for f in `seq 1 50`
do
for c in 1 2 3
do

cat << EOF > ${QSUB}$f$c".qsub"
#!/bin/bash
#PBS -q ${que}
#PBS -N $TASK$f$c
#PBS -l select=1:ncpus=${ncpu_each}:mem=${nmem}gb
#PBS -e $mainDIR/LOGs2/$LOGPATH/$TASK$f$c.log
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

echo "setwd(\"${mainDIR}\")" > "scl_Rscripts/${TASK}$data$f$c.R"
cat "/aptmp/hiroakif/noguchi_260514/Script/${TASK}.R" >> "scl_Rscripts/${TASK}$data$f$c.R"

sed -i "25s/.*/seed=${f}/" "scl_Rscripts/${TASK}$data$f$c.R"
sed -i "26s/.*/connect=${c}/" "scl_Rscripts/${TASK}$data$f$c.R"
sed -i "27s/.*/n.core=${ncpu}/" "scl_Rscripts/${TASK}$data$f$c.R"

Rscript --vanilla "scl_Rscripts/${TASK}$data$f$c.R"

EOF

qsub ${QSUB}$f$c".qsub"

done
done

#####
TASK=${dir_05_01_00}
QSUB="QSUBs/$TASK" 
data=260514
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
echo "set.seed(1234)" >> "scl_Rscripts/${TASK}$data.R"
cat "/aptmp/hiroakif/noguchi_260514/Script/${TASK}.R" >> "scl_Rscripts/${TASK}$data.R"

sed -i "8s/.*/dir_05_01='05_01_gLV_simulation_Nsp80_nonoise_260514'/" "scl_Rscripts/${TASK}$data.R"

Rscript --vanilla "scl_Rscripts/${TASK}$data.R"

EOF

qsub ${QSUB}".qsub"

################################################
###05_01

TASK=${dir_05_02}
QSUB="QSUBs/$TASK" 
data=260514
LOGPATH="LOGs_${TASK}_${data}"
mkdir -p LOGs2/$LOGPATH

for f in `seq 1 15`
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
cat "/aptmp/hiroakif/noguchi_260514/Script/${TASK}.R" >> "scl_Rscripts/${TASK}$data$f.R"

sed -i "42s/.*/r=${f}/" "scl_Rscripts/${TASK}$data$f.R"
sed -i "43s/.*/dir_05_01='05_01_gLV_simulation_Nsp80_nonoise_260514'/" "scl_Rscripts/${TASK}$data$f.R"
sed -i "44s/.*/dir_05_01_00='${dir_05_01_00}'/" "scl_Rscripts/${TASK}$data$f.R"
sed -i "45s/.*/n.core=${ncpu}/" "scl_Rscripts/${TASK}$data$f.R"

Rscript --vanilla "scl_Rscripts/${TASK}$data$f.R"

EOF

qsub ${QSUB}$f".qsub"

done

################################################
################################################
###05_03


TASK=${dir_05_03}
QSUB="QSUBs/$TASK" 
data=260514
ncpu=1
LOGPATH="LOGs_${TASK}_${data}"
mkdir -p LOGs2/$LOGPATH

for r in `seq 4 3000`
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

## load app
source /etc/profile.d/modules.sh
module load $MODULE
export OMP_NUM_THREADS=1

echo "setwd(\"${mainDIR}\")" > "scl_Rscripts/${TASK}$data$r.R"
cat "/aptmp/hiroakif/noguchi_260514/Script/${TASK}.R" >> "scl_Rscripts/${TASK}$data$r.R"

sed -i "43s/.*/dir_05_01='05_01_gLV_simulation_Nsp80_nonoise_260514'/" "scl_Rscripts/${TASK}$data$r.R"
sed -i "44s/.*/dir_05_01_00='${dir_05_01_00}'/" "scl_Rscripts/${TASK}$data$r.R"
sed -i "45s/.*/n.core=${ncpu}/" "scl_Rscripts/${TASK}$data$r.R"
sed -i "46s/.*/cur_seed=${r}/" "scl_Rscripts/${TASK}$data$r.R"

Rscript --vanilla "scl_Rscripts/${TASK}$data$r.R"

EOF

qsub ${QSUB}$r".qsub"

done
################################################
#####################################################
#05_04
TASK=${dir_05_04}
QSUB="QSUBs/$TASK" 
data=260514
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
echo "set.seed(1234)" > "scl_Rscripts/${TASK}$data.R"
cat "/aptmp/hiroakif/noguchi_260514/Script/${TASK}.R" >> "scl_Rscripts/${TASK}$data.R"


sed -i "42s/.*/dir_05_01 = '05_01_gLV_simulation_Nsp80_nonoise_260514'/" "scl_Rscripts/${TASK}$data.R"
sed -i "43s/.*/dir_05_02 = '${dir_05_02}'/" "scl_Rscripts/${TASK}$data.R"
sed -i "44s/.*/dir_05_01_00='${dir_05_01_00}'/" "scl_Rscripts/${TASK}$data$f.R"

Rscript --vanilla "scl_Rscripts/${TASK}$data.R"

EOF

qsub ${QSUB}".qsub"

#####################################################
#####################################################
#05_04
TASK=${dir_05_05}
QSUB="QSUBs/$TASK" 
data=260514
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
echo "set.seed(1234)" >> "scl_Rscripts/${TASK}$data.R"
cat "/aptmp/hiroakif/noguchi_260514/Script/${TASK}.R" >> "scl_Rscripts/${TASK}$data.R"

sed -i "20s/.*/dir_05_02 = '${dir_05_02}'/" "scl_Rscripts/${TASK}$data.R"
sed -i "21s/.*/dir_05_03 = '05_03_randELA_withRA_4step_Nsp80_nonoise_260527'/" "scl_Rscripts/${TASK}$data.R"
sed -i "22s/.*/dir_05_01_00='${dir_05_01_00}'/" "scl_Rscripts/${TASK}$data.R"
sed -i "23s/.*/n.core=${ncpu}/" "scl_Rscripts/${TASK}$data.R"

Rscript --vanilla "scl_Rscripts/${TASK}$data.R"

EOF

qsub ${QSUB}".qsub"

################################################
###05_03

#####################################################
#05_04
TASK="05_07_01_merge_result_eachseed_Nsp80_nonoise_260514"
SAVEDIR="05_07_01_merge_result_eachseed_Nsp80_nonoise_260514"
QSUB="QSUBs/$TASK" 
data=260514
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
SAVEDIR=$SAVEDIR

## load app
source /etc/profile.d/modules.sh
module load $MODULE
export OMP_NUM_THREADS=1

echo "setwd(\"${mainDIR}\")" > "scl_Rscripts/${TASK}$data.R"
echo "set.seed(1234)" > "scl_Rscripts/${TASK}$data.R"
cat "/aptmp/hiroakif/noguchi_260514/Script/${TASK}.R" >> "scl_Rscripts/${TASK}$data.R"

sed -i "6s/.*/dir_05_04 = '${dir_05_04}'/" "scl_Rscripts/${TASK}$data.R"
sed -i "6s/.*/dir_05_04 = '${dir_05_04}'/" "scl_Rscripts/${TASK}$data.R"
sed -i "7s/.*/dir_05_05 = '${dir_05_05}'/" "scl_Rscripts/${TASK}$data.R"
sed -i "8s/.*/dir_05_06 = '${dir_05_06}'/" "scl_Rscripts/${TASK}$data.R"
sed -i "9s/.*/dir_05_01_00='${dir_05_01_00}'/" "scl_Rscripts/${TASK}$data.R"
sed -i "10s/.*/n.core=${ncpu}/" "scl_Rscripts/${TASK}$data.R"
sed -i "11s/.*/save.dir='${SAVEDIR}'/" "scl_Rscripts/${TASK}$data.R"

Rscript --vanilla "scl_Rscripts/${TASK}$data.R"

EOF

qsub ${QSUB}".qsub"
