#!/bin/bash
#SBATCH --job-name filter
#SBATCH --array=1-N
#SBATCH --account amalaspi_popgen
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --export=NONE
#SBATCH --mem 2G
#SBATCH --time=2:00:00
#SBATCH --error=splitQual_%a.er
#SBATCH --output=splitQual_%a.out

#filter imputed files based on GP and MAF (here GP>0.80 and MAF>5%)

#FILE with sample ids (one per line, N lines)
LST=samples.txt

#extract sample id
SPL=$(cat $LST | head -n ${SLURM_ARRAY_TASK_ID} | tail -n 1)

#Concatenated (all chromosomes in one file) imputed data 
INP=/step2_ligate/ligated/autosomes.ligated.vcf.gz

#output
OUT=/filtered/${SPL}.${GP}.bcf


# split vcf by sample
bcftools view -s $SID $INP -i'INFO/RAF[0]>0.05 && INFO/RAF[0]<0.95' | bcftools view -i'FORMAT/GP[0:*]>0.99' -O b -o $OUT
bcftools index -f $OUT
