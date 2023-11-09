#!/bin/bash
#SBATCH --job-name=BL_a3sd
#SBATCH --output=BL_a3sd.out
#SBATCH --error=BL_a3sd.err
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

populations -P /hb/groups/bernardi_lab/eric/Lionfish/PMI_genome/ref_stacks/ -O /hb/groups/bernardi_lab/eric/Lionfish/PMI_genome/structure_ave3sd/BL/ --popmap /hb/groups/bernardi_lab/eric/Lionfish/Files/pop_map_PMI_66_RS_Med.txt -B /hb/groups/bernardi_lab/eric/Lionfish/Files/ave_3sd_whitelist_ref_gen_PMI_528.txt -r 0.80 -t 24  --write_single_snp --structure --genepop --fstats --vcf

