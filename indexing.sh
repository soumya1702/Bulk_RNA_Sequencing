#!/bin/bash
#SBATCH -J samples
#SBATCH -p general
#SBATCH -o indexing.txt
#SBATCH -e indexing.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=syennapu@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=84:00:00
#SBATCH --mem=100G
#SBATCH -A students
# Load modules
module load 'star'

GENOMEDIR="/N/scratch/syennapu/HTP_PROJECT/"
mkdir -p $GENOMEDIR/indexingfiles
echo "$NOW1 STAR: Generate Genome Index 1"
echo "======================================"

STAR --runThreadN 12 --runMode genomeGenerate --genomeDir $GENOMEDIR/indexingfiles --genomeFastaFiles $GENOMEDIR/GRCh38.p13.genome.fa --sjdbGTFfile $GENOMEDIR/gencode.v32.chr_patch_hapl_scaff.annotation.gtf --sjdbOverhang 149

echo "======================================"
