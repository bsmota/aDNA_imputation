#!/bin/bash
#SBATCH --array=1-22
#SBATCH --job-name ligate
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --export=NONE
#SBATCH --mem 2G
#SBATCH --time=1:00:00
#SBATCH --error=ligate_%a.er
#SBATCH --output=ligate_%a.out

#this script ligates all imputed chunks in one chromosome

#chromosome
CHR=$SLURM_ARRAY_TASK_ID

#software location of GLIMPSE_ligate
BIN=/GLIMPSE-v1.1.1/ligate/bin/GLIMPSE_ligate


LST=/step2_ligate/lists/chr$CHR.list.txt
ls /step1_impute/output/chr$CHR.reg*.bcf > $LST

OUT=/step2_ligate/ligated/chr$CHR.ligated.vcf.gz
OUL=/step2_ligate/ligated/chr$CHR.ligated.log

#Run GLIMPSE_ligate
$BIN --input $LST --output $OUT --log $OUL
#index output file
bcftools index -f $OUT
