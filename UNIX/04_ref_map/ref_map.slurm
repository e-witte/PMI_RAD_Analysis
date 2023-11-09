#!/bin/bash
#SBATCH --output=pmi_sam_bam.out
#SBATCH --error=pmi_sam_bam.err
#SBATCH --partition=128x24
#SBATCH --nodes=1
#SBATCH --job-name=pmi_sam_bam
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ewitte@ucsc.edu  #!!change this to your own email
#SBATCH --mem=22GB
#SBATCH --ntasks-per-node=24

# Above are all the specifics SLURM needs to run your specific job. 

# --ntasks-per-node will be used in doParallel.R to specify the number
# of cores to use on the machine. Using 24

# now load programs needed

module load bowtie
module load samtools
module load stacks2

ref_map.pl -T 24 --samples /hb/groups/bernardi_lab/eric/Lionfish/PMI_genome/bam_files --popmap /hb/groups/bernardi_lab/eric/Lionfish/Files/pop_map_PMI_66_RS_Med.txt -o /hb/groups/bernardi_lab/eric/Lionfish/PMI_genome/ref_stacks/
