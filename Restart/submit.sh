#!/bin/bash --login

#$ -N NaY.N2.100kpa        # Name of the job
#$ -cwd                 # Job will run from the current directory
#$ -l s_rt=110:00:00
#$ -l avx

# To prevent segmentation fault 
ulimit -s unlimited

# Load the required version
module load libs/intel-17.0/fftw3/serial/double/3.3.8 libs/lapack/3.5.0/gcc-4.8.5
export RASPA_DIR=/mnt/iusers01/ceas01/mbdxqzd2/RASPA2_17122021
$RASPA_DIR/bin/simulate

