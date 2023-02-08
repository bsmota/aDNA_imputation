#!/bin/bash
#SBATCH --job-name=call
#SBATCH --array=1-N*22
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=2:00:00
#SBATCH --output=call_%a.out
#SBATCH --error=call_%a.er

#generate genotype likelihoods for N samples, one sample and one chromosome at a time
#file with chromosomes, sample ids and bam filenames (file has N*22 lines)
LST=input_data.txt

#Chromosome
CHR=$(cat $LST | head -n ${SLURM_ARRAY_TASK_ID} | tail -n 1 | awk '{print $1}')

#Sample id
SPL=$(cat $LST | head -n ${SLURM_ARRAY_TASK_ID} | tail -n 1 | awk '{print $2}')

#Bam filenames
BAM=$(cat $LST | head -n ${SLURM_ARRAY_TASK_ID} | tail -n 1 | awk '{print $3}')

#path to directory with bam files
BAMDIR=/path/toBams

#path to reference genome
FASTA=reference_genome.fa

#List of variants in reference panel 
#(#CHROM\tPOS\tID\tREF\tALT, it can be generated with from genotypes file as bcftools view -G reference.panel.genotypes.bcf)
VPOS=chr${CHR}.referencePanel.vcf.gz

#List of variants in reference panel 
#(CHROM\tPOS\tREF,ALT, no header)
TSV=chr${CHR}.referencePanel.tsv.gz

SPL1=/calls/chr${CHR}/${SPL}.calls.spl
OUT=/calls/chr${CHR}/${SPL}.calls.bcf
TMP=/calls/chr${CHR}/${SPL}.calls.tmp.bcf

#Call genotypes using bcftools
echo ${BAMDIR}/${BAM} $SPL > $SPL1
bcftools mpileup -f $FASTA -I -E -a 'FORMAT/DP' --ignore-RG -T $VPOS -r ${CHR} ${BAMDIR}/${BAM} | bcftools call -Aim -C alleles -T $TSV -Oz -o $OUT
bcftools reheader -s $SPL1 -o $TMP $OUT
mv $TMP $OUT
bcftools index -f $OUT




