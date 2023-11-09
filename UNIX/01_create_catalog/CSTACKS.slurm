#!/bin/bash
#SBATCH --job-name=ustacks
#SBATCH --output=cstacks_PMI.out
#SBATCH --error=cstacks_PMI.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ewitte@ucsc.edu  #!!change this to your own email
#SBATCH --partition=128x24
#SBATCH --nodes=1
#SBATCH --mem=22GB
#SBATCH --ntasks-per-node=24

# Above are all the specifics SLURM needs to run your specific job.

# --ntasks-per-node will be used in doParallel.R to specify the number
# of cores to use on the machine. Using 24
# now load programs needed

module load stacks2

cstacks -P ./ -M /hb/groups/bernardi_lab/eric/Lionfish/Files/pop_map_PMI_10.txt -n 4 -p 15 
#cstacks creates the catalog with the 10 best samples
