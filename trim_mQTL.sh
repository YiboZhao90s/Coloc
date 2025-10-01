#!/bin/bash
#
# SLURM directives:
#SBATCH --job-name=Chr11_mQTL
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=200G
#SBATCH --time=24:00:00
#SBATCH --export=NONE

# commands to run
cd /scratch/gen1/yz735/coloc
zcat Lung.mQTLs.regular.txt.gz | awk '{ if (substr($2,1,5)=="chr11") print }' > Lung.mQTLs.regular_chr11.txt
zcat Lung.mQTLs.conditional.txt.gz | awk '{ if (substr($2,1,5)=="chr11") print }' > Lung.mQTLs.conditional_chr11.txt



# list loaded modules, host we are running on and date
hostname
date


