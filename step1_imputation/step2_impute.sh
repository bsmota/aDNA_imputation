#!/bin/bash
#SBATCH --array=1-X
#SBATCH --job-name ImpLC
#SBATCH --account=amalaspi_popgen
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --export=NONE
#SBATCH --mem 15G
#SBATCH --output=impute_%a.out
#SBATCH --error=impute_%a.er

#this script allows to impute the X chunks in parallel

#software location
IMP=/GLIMPSE-v1.1.1/phase/bin/GLIMPSE_phase

#chunks coordinates and info
LST=/coordinates/referencel.panel.chr${CHR}.txt

DATA=$(cat $LST | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)

#chunk number
IDG=$(echo $DATA | cut -d" " -f1)
#chromosome
CHR=$(echo $DATA | cut -d" " -f2)
#chunk with buffer coordinates
IREG=$(echo $DATA | cut -d" " -f3)
#chunk without buffer coordinate
OREG=$(echo $DATA | cut -d" " -f4)

#Reference panel
REF=chr${CHR}.referencePanel_genotypes.bcf

#genotype likelihoods
GLS=/calls/chr${CHR}/chr${CHR}.merged.bcf

#genetic map (if build 37 of the reference genome)
MAP=/GMAP_shapeit4/chr${CHR}.b37.gmap.gz

#output file
OUTD=/step1_impute/output/chr$CHR.reg$IDG.bcf

#Run GLIMPSE_phase to impute variants in chunk
$IMP --input $GLS --reference $REF --input-region $IREG --output-region $OREG --map $MAP --output $OUTD

bcftools index -f $OUTD

