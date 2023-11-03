#!/bin/bash
#SBATCH --output=pmi_sam.out
#SBATCH --error=pmi_sam.err
#SBATCH --partition=128x24
#SBATCH --nodes=1
#SBATCH --job-name=pmi_sam
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

#index your reference genome

bowtie2-build /hb/groups/bernardi_lab/eric/Lionfish/PMI_genome/genome/GCA_947000775.1.fasta index