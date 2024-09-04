#!/bin/bash
#SBATCH -J counts
#SBATCH -p general
#SBATCH -o feature.txt
#SBATCH -e counts.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=syennapu@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=84:00:00
#SBATCH --mem=100G
#SBATCH -A students
# Load modules

featureCounts -a /N/scratch/syennapu/HTP_PROJECT/gencode.v32.chr_patch_hapl_scaff.annotation.gtf -o /N/scratch/syennapu/HTP_PROJECT/counts.txt -T 8 -p *.bam
