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

#Create a variable called "files" with a list of samples for bowtie to convert into SAM files.

files="PMI_MAS_031401
PMI_MAS_031402
PMI_MAL_051901
PMI_MAL_051902
PMI_MAL_051903
PMI_MAL_051904
PMI_MAL_051905_2
PMI_MAL_051906
PMI_MAL_051907_F
PMI_MAL_051908_2
PMI_LEB_0215633
PMI_LEB_0216630
PMI_LEB_0516634
PMI_LEB_635
PMI_LEB_651
PMI_LEB_733
PMI_LEB_734
PMI_LEB_743
PMI_LEB_744
PMI_LEB_746
PMI_CYP_041801
PMI_CYP_041803
PMI_CYP_041804
PMI_CYP_041805
PMI_CYP_041806
PMI_CYP_041808
PMI_CYP_041809
PMI_CYP_041810
PMI_CYP_041811
PMI_CYP_041812
PMI_CYP_041813
PMI_CYP_041814
PMI_CYP_041817
PMI_CYP_041818
PMI_CYP_041819
PMI_CYP_041820
PMI_CYP_041821
PMI_CYP_041826
PMI_CYP_041829
PMI_CYP_041830
PMI_CYP_041833
PMI_CYP_041834
PMI_CYP_041835
PMI_CYP_041836
PMI_CYP_041838
PMI_CYP_041839
PMI_CYP_041840
PMI_CYP_041841
PMI_CYP_041842
PMI_CYP_041843
PMI_CYP_041845
PMI_CYP_041847
PMI_CYP_041849
PMI_CYP_041851
PMI_CYP_061601
PMI_CYP_061602
PMI_CYP_071606
PMI_RHO_041840
PMI_RHO_041841
PMI_RHO_041844
PMI_RHO_041845
PMI_RHO_041846
PMI_RHO_041847
PMI_RHO_041848
PMI_RHO_041850
PMI_RHO_041852"

#for loop to execute the bowtie2 command for each sample in "files"
for sample in $files
do
   bowtie2 -q -p 24 --phred33 --very-sensitive -x index -U /hb/groups/bernardi_lab/eric/Lionfish/RAD_samples_PMI/$sample.fq -S /hb/groups/bernardi_lab/eric/Lionfish/PMI_genome/sam_files/$sample.sam
#sort and convert SAM files to BAM files using samtools   
   samtools view -b /hb/groups/bernardi_lab/eric/Lionfish/PMI_genome/sam_files/$sample.sam | samtools sort > /hb/groups/bernardi_lab/eric/Lionfish/PMI_genome/bam_files/$sample.bam
#remove the large .sam file
rm $sample.sam
done