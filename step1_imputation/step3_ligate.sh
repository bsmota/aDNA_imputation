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

#software location of GLIMPSE_ligate
BIN=/GLIMPSE-v1.1.1/ligate/bin/GLIMPSE_ligate

#chromosome
CHR=$SLURM_ARRAY_TASK_ID

LST=/step2_ligate/lists/chr$CHR.list.txt

ls /step1_impute/output/1000GP_nygc_umich.chr$CHR.reg*.vcf.gz > $LST


OUT=${PDIR}/step2_ligate/ligated/chr$CHR.ligated.vcf.gz
OUL=${PDIR}/step2_ligate/ligated/chr$CHR.ligated.log

$BIN --input $LST --output $OUT --log $OUL
bcftools index -f $OUT
