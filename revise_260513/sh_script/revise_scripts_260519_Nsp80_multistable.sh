#!bin/bash
que=APC
ncpu=14
ncpu_each=4
nmem=24
data=260516

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
dir_07_02="Basedata/07_02_evaluate_multistable_dynamics_260519"
dir_07_03="07_03_ELA_withRA_multistable_Nsp80_260519"
dir_07_04="07_04_randELA_withRA_multistable_Nsp80_260519"
dir_07_05="07_05_Zconv_ELA_withRA_4step_Nsp80_multistable_260519"

################################################
###

TASK=${dir_07_03}
QSUB="QSUBs/$TASK" 
data=260516
LOGPATH="LOGs_${TASK}_${data}"
mkdir -p LOGs2/$LOGPATH

for f in `seq 1 3`
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
echo "top=${f}" >>  "scl_Rscripts/${TASK}$data$f.R"
echo "dir_07_02='${dir_07_02}'" >>  "scl_Rscripts/${TASK}$data$f.R"
echo "n.core=${ncpu}" >>  "scl_Rscripts/${TASK}$data$f.R"

cat "/aptmp/hiroakif/noguchi_260514/Script/${TASK}.R" >> "scl_Rscripts/${TASK}$data$f.R"

# sed -i "4s/.*/top=${f}/" "scl_Rscripts/${TASK}$data$f.R"
# sed -i "5s/.*/dir_07_02='${dir_07_02}'/" "scl_Rscripts/${TASK}$data$f.R"
# sed -i "6s/.*/n.core=${ncpu}/" "scl_Rscripts/${TASK}$data$f.R"

Rscript --vanilla "scl_Rscripts/${TASK}$data$f.R"

EOF

qsub ${QSUB}$f".qsub"

done

################################################
################################################
###05_03


TASK=${dir_07_04}
QSUB="QSUBs/$TASK" 
data=260516
LOGPATH="LOGs_${TASK}_${data}"
mkdir -p LOGs2/$LOGPATH

for r in `seq 1 5`
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
echo "dir_07_02='${dir_07_2}'" >> "scl_Rscripts/${TASK}$data$r.R"
echo "n.core=1" >> "scl_Rscripts/${TASK}$data$r.R"
echo "seed=${r}" >> "scl_Rscripts/${TASK}$data$r.R"

cat "/aptmp/hiroakif/noguchi_260514/Script/${TASK}.R" >> "scl_Rscripts/${TASK}$data$r.R"

# sed -i "4s/.*/dir_07_02='${dir_07_2}'/" "scl_Rscripts/${TASK}$data$r.R"
# sed -i "5s/.*/n.core=1/" "scl_Rscripts/${TASK}$data$r.R"
# sed -i "6s/.*/seed=${r}/" "scl_Rscripts/${TASK}$data$r.R"

Rscript --vanilla "scl_Rscripts/${TASK}$data$r.R"

EOF

qsub ${QSUB}$r".qsub"

done
################################################
#####################################################
#
TASK=${dir_07_05}
QSUB="QSUBs/$TASK" 
data=260518
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
echo "n.core=1" >> "scl_Rscripts/${TASK}$data$r.R"
echo "set.seed(1234)" > "scl_Rscripts/${TASK}$data.R"
cat "/aptmp/hiroakif/noguchi_260514/Script/${TASK}.R" >> "scl_Rscripts/${TASK}$data.R"

sed -i "4s/.*/dir_07_03 = '${dir_07_03}'/" "scl_Rscripts/${TASK}$data.R"
sed -i "5s/.*/dir_07_04 = '${dir_07_04}'/" "scl_Rscripts/${TASK}$data.R"

Rscript --vanilla "scl_Rscripts/${TASK}$data.R"

EOF

qsub ${QSUB}".qsub"

#####################################################


#####################################################
#05_04
TASK="07_06_merge_result_eachseed_Nsp80_multistable"
SAVEDIR="07_06_merge_result_eachseed_Nsp80_multistable"
QSUB="QSUBs/$TASK" 
data=260519
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

sed -i "6s|.*|dir_07_02 = '${dir_07_02}'|" "scl_Rscripts/${TASK}$data.R"
sed -i "7s/.*/dir_07_05 = '${dir_07_05}'/" "scl_Rscripts/${TASK}$data.R"
sed -i "10s/.*/n.core=${ncpu}/" "scl_Rscripts/${TASK}$data.R"
sed -i "11s/.*/save.dir='${SAVEDIR}'/" "scl_Rscripts/${TASK}$data.R"

Rscript --vanilla "scl_Rscripts/${TASK}$data.R"

EOF

qsub ${QSUB}".qsub"
