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
###05_01

TASK="06_01_gLV_samplesize_validation_260108"
QSUB="QSUBs/$TASK" 
data=260108
LOGPATH="LOGs_${TASK}_${data}"
mkdir -p LOGs2/$LOGPATH

for s in `seq 9 9`
do
for f in `seq 1 3`
do
for c in `seq 1 4`
do

cat << EOF > ${QSUB}$f$c$s".qsub"
#!/bin/bash
#PBS -q ${que}
#PBS -N $TASK$f$c$s
#PBS -l select=1:ncpus=${ncpu_each}:mem=${nmem}gb
#PBS -e $mainDIR/LOGs2/$LOGPATH/$TASK$f$c$s.log
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

echo "setwd(\"${mainDIR}\")" > "scl_Rscripts/${TASK}$data$f$c$s.R"
cat "/aptmp/mikihito/02_ELA_v060_241113_2/Script/${TASK}.R" >> "scl_Rscripts/${TASK}$data$f$c$s.R"

sed -i "21s/.*/n.core=${ncpu}/" "scl_Rscripts/${TASK}$data$f$c$s.R"
sed -i "22s/.*/seed=${f}/" "scl_Rscripts/${TASK}$data$f$c$s.R"
sed -i "23s/.*/nsp=${c}/" "scl_Rscripts/${TASK}$data$f$c$s.R"
sed -i "24s/.*/nsamp=${s}/" "scl_Rscripts/${TASK}$data$f$c$s.R"

Rscript --vanilla "scl_Rscripts/${TASK}$data$f$c$s.R"

EOF

qsub ${QSUB}$f$c$s".qsub"

done
done
done

#####################################################
#05_04
TASK="06_02_gLV_graphics_samplesize_validation"
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
echo "set.seed(1234)" > "scl_Rscripts/${TASK}$data.R"
cat "/aptmp/mikihito/02_ELA_v060_241113_2/Script/${TASK}.R" >> "scl_Rscripts/${TASK}$data.R"

sed -i "10s/.*/dir_06_01 = '06_01_gLV_samplesize_validation_260108'/" "scl_Rscripts/${TASK}$data.R"

Rscript --vanilla "scl_Rscripts/${TASK}$data.R"

EOF

qsub ${QSUB}".qsub"
