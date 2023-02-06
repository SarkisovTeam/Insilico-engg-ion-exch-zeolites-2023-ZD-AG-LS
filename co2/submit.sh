#!/bin/bash --login
#$ -cwd                 # Job will run from the current directory
#$ -l avx

# To prevent segmentation fault 
ulimit -s unlimited

# Load the required version
export RASPA_DIR=/mnt/iusers01/ceas01/mbdxqzd2/RASPA2_17122021
$RASPA_DIR/bin/simulate

