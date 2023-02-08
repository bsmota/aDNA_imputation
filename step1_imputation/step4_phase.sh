#!/bin/bash
#SBATCH --job-name phase
#SBATCH --array=1-22
#SBATCH --account amalaspi_popgen
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --export=NONE
#SBATCH --mem 2G
#SBATCH --time=1:00:00
#SBATCH --error=phase_%a.er
#SBATCH --output=phase_%a.out

CHR=${SLURM_ARRAY_TASK_ID}

#location of GLIMPSE phasing tool
BIN=/GLIMPSE-v1.1.1/sample/bin/GLIMPSE_sample_static

#input (ligated imputed data)
VCF=/step2_ligate/ligated/chr$CHR.ligated.vcf.gz
#output
OUT=/step3_phased/chr$CHR.phased.vcf.gz

#phase the data in one chromosome
$BIN --input ${VCF} --solve --output ${OUT}
bcftools index -f $OUT
