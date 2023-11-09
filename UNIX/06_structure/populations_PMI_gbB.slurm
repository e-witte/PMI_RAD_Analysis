#!/bin/bash
#SBATCH --job-name=populations
#SBATCH --output=populations_PMI.out
#SBATCH --error=populations_PMI.err
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

populations -P ./ -M /hb/groups/bernardi_lab/eric/Lionfish/Files/pop_map_PMI_66_RS_Med.txt -B /hb/groups/bernardi_lab/eric/Lionfish/Files/PMI_whitelist_GB_184.txt -r 0.8 --write_single_snp --vcf --genepop --structure --fstats -t 15
