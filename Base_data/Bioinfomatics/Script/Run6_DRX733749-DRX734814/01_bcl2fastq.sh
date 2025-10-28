####################################################################
# Current directory should be the output run folder of the MiSeq run.

experiment_name=Run6

####################################################################

configureBclToFastq.pl --fastq-cluster-count 0 --use-bases-mask Y270n,Y8,Y8,Y30n --input-dir ./Data/Intensities/BaseCalls --output-dir $experiment_name

cd ./$experiment_name
make -j32

####################################################################

