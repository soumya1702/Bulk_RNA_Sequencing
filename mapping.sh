#!/bin/bash
#SBATCH -J mapping
#SBATCH -p general
#SBATCH -o mapping.txt
#SBATCH -e mapping.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=syennapu@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=48:00:00
#SBATCH --mem=100G
#SBATCH -A students

#loading star module
module load 'star'
#loading python module
module load 'python/3.9.8'

#creating a new directory for output files
mkdir -p /N/scratch/syennapu/HTP_PROJECT/mappingfiles
chmod +rwx /N/scratch/syennapu/HTP_PROJECT/mappingfiles

GENOMEDIR="/N/scratch/syennapu/HTP_PROJECT/"
FASTQ_FILES="/N/scratch/syennapu/HTP_PROJECT/fastqfiles"

for ((i=176; i<=217; i++)); do
  SAMPLE_NAME="SRR23246${i}"

  # Check if both input files exist
  if [[ ! -f "$FASTQ_FILES/${SAMPLE_NAME}_1.fastq.gz" || ! -f "$FASTQ_FILES/${SAMPLE_NAME}_2.fastq.gz" ]]; then
    echo "sample $SAMPLE_NAME does not exist"
    continue  # Skip the current iteration
  fi

  STAR --runThreadN 12 --genomeDir $GENOMEDIR/indexingfiles --readFilesIn $FASTQ_FILES/${SAMPLE_NAME}_1.fastq.gz $FASTQ_FILES/${SAMPLE_NAME}_2.fastq.gz \
      --readFilesCommand zcat \
      --quantMode GeneCounts \
      --outFileNamePrefix /N/scratch/syennapu/HTP_PROJECT/mappingfiles/${SAMPLE_NAME} --outSAMmapqUnique 60 --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx
done
