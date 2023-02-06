#!/bin/bash --login

#$ -N ProcOpt       # Name of the job
#$ -cwd             # Job will run from the current directory
#$ -pe smp.pe 24    # run using 32 cores
#$ -l avx       # run using the skylake node 

# load matlab module
module load apps/binapps/matlab/R2020a

# running the compiled MATLAB code
# follow directions in Tutorial_Run_on_CSF3.docx 
./run_FullOptimization.sh $MATLAB_HOME
