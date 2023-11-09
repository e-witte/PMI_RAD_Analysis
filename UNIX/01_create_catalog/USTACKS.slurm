#!/bin/bash
#SBATCH --job-name=ustacks
#SBATCH --output=ustacks_PMI.out
#SBATCH --error=ustacks_PMI.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ewitte@ucsc.edu  #!!change this to your own email
#SBATCH --partition=128x24
#SBATCH --nodes=1
#SBATCH --mem=22GB
#SBATCH --ntasks-per-node=24

module load stacks2

FILES=("PMI_MAL_051901"
"PMI_MAL_051902"
"PMI_MAL_051905_2"
"PMI_MAL_051907_F"
"PMI_MAL_051908_2"
"PMI_LEB_635"
"PMI_CYP_041806"
"PMI_CYP_041830"
"PMI_CYP_041835"
"PMI_CYP_041843")

numbers=("1"
"2"
"3"
"4"
"5"
"6"
"7"
"8"
"9"
"10")

length=${#FILES[@]}

for ((i=0;i<$length;i++));
do
ustacks -f /hb/groups/bernardi_lab/eric/Lionfish/RAD_samples_PMI/${FILES[$i]}.fq -i ${numbers[$i]} -o ./ -M 3 -m 3 -p8

done
