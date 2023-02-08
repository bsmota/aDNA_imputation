#!/bin/bash
#SBATCH --job-name=chunk
#SBATCH --array=1-22
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --time=60:00
#SBATCH --mem=4G
#SBATCH --error=chunking_%a.er
#SBATCH --output=chunking_%a.out 

#generate chunks for variants in reference panel for each chromosome (where imputation will act)
#chromosome
CHR=${SLURM_ARRAY_TASK_ID}

#Path to GLIMPSE_chunk 
BIN=/GLIMPSE-v1.1.1/chunk/bin/GLIMPSE_chunk

#Reference panel
REF=chr${CHR}.reference.panel.bcf

#output file, chunk coordinates
OUT=/coordinates/referencel.panel.chr${CHR}.txt


#Generate chunks, with window size of 1Mb with a buffer of 200kb size
$BIN --input $REF --region chr$CHR --window-size 1000000 --buffer-size 200000  --output $OUT

