#!/bin/bash
#SBATCH --job-name=gstacks
#SBATCH --output=gstacks_PMI.out
#SBATCH --error=gstacks_PMI.err
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

gstacks -P ./ -M /hb/groups/bernardi_lab/eric/Lionfish/Files/pop_map_PMI_66_RS_Med.txt --model marukilow -t 8
