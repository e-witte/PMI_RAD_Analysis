#!/bin/bash
#SBATCH --job-name=BL_a3sd_struc
#SBATCH --output=BL_a3sd_struc.out
#SBATCH --error=BL_a3sd_struc.err
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

module load structure

structure -K 2 -L 23382 -N 66 -i /hb/groups/bernardi_lab/eric/Lionfish/PMI_genome/structure_ave3sd/BL/populations.structure_Ba3sd_PMI_66 -o PMI_66_23231_Ba3sd_K2_output1

