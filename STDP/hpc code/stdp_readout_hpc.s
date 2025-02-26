#!/bin/bash
#
#SBATCH --array=1-250
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --mem=4GB
#SBATCH --job-name=k11=0.045
#SBATCH --mail-user=wt2113@nyu.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=slurm..out  # %A is job ID, %a is array index

module purge
module load matlab/2021b

# Navigate to the directory where the MATLAB script is located
cd $HOME/synaptic_plasticity/decrease_k11/pre=1.3/k11=0.045

# Execute MATLAB script without GUI
matlab -nodisplay -nosplash -nodesktop -r "run('$HOME/synaptic_plasticity/decrease_k11/pre=1.3/k11=0.045/stdp_readout_hpc.m'); exit;"