#!/bin/bash
#SBATCH --job-name=call
#SBATCH --array=1-N*22
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=8:00:00
#SBATCH --output=val0_%a.out
#SBATCH --error=val0_%a.er

#file with three columns: chromosome, sample id, bam filename
LST=calling_data.txt

#chromosome
CHR=$(cat $LST | head -n ${SLURM_ARRAY_TASK_ID} | tail -n 1 | awk '{print $1}')
#sample id
SPL=$(cat $LST | head -n ${SLURM_ARRAY_TASK_ID} | tail -n 1 | awk '{print $2}')
#bam filename
BAM=$(cat $LST | head -n ${SLURM_ARRAY_TASK_ID} | tail -n 1 | awk '{print $3}')

#location of bam files
BAMDIR=/bams

FASTA=reference.genome.fa #reference genome data was aligned to

#Position files
#List of variants in reference panel 
#(#CHROM\tPOS\tID\tREF\tALT, it can be generated with from genotypes file as bcftools view -G reference.panel.genotypes.bcf)
VPOS=chr${CHR}.referencePanel.vcf.gz

#List of variants in reference panel 
#(CHROM\tPOS\tREF,ALT, no header)
TSV=chr${CHR}.referencePanel.tsv.gz

#Output VCF and LOG files
SPL1=/calls/chr${CHR}/${SPL}.raw.calls.spl
OUT=/calls/chr${CHR}/${SPL}.raw.calls.vcf.gz
TMP=/calls/chr${CHR}/${SPL}.raw.calls.tmp.vcf.gz

#Call genotypes using bcftools
echo ${BAMDIR}/${BAM} $SPL > $SPL1
bcftools mpileup -f $FASTA -I -E -a 'FORMAT/DP' --ignore-RG -T $VPOS  -r $CHR ${BAMDIR}/${BAM} | bcftools call -Aim -C alleles -T $TSV -Oz -o $OUT
bcftools reheader -s $SPL1 -o $TMP $OUT
mv $TMP $OUT
bcftools index -f $OUT
