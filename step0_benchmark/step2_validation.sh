#!/bin/bash
#SBATCH --job-name=val0
#SBATCH --array=1-N*22
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=8:00:00
#SBATCH --output=val0_%a.out
#SBATCH --error=val0_%a.er

#generate validation data for imputation experiments
#1. genotype calling for high-coverage genomes with minimum mapping quality of 20,  base quality of 30 and option -C 50
#2. keep only sites in 1000G accessible genome mask
#3. remove sites in repeat regions
#4. remove sites that are outlier in depth and, but imposing that minimum depth>8

#file with three columns: chr, sample id, bam filename
LST=/calling_data.txt

#chromosome
CHR=$(cat $LST | head -n ${SLURM_ARRAY_TASK_ID} | tail -n 1 | awk '{print $3}')
#sample id
SPL=$(cat $LST | head -n ${SLURM_ARRAY_TASK_ID} | tail -n 1 | awk '{print $1}')
#bam filename
BAM=$(cat $LST | head -n ${SLURM_ARRAY_TASK_ID} | tail -n 1 | awk '{print $2}')

#Pointer to project
FASTA=hs.build37.1.fa

#Position file	
#List of variants in reference panel 
#(#CHROM\tPOS\tID\tREF\tALT, it can be generated with from genotypes file as bcftools view -G reference.panel.genotypes.bcf)
VPOS=chr${CHR}.referencePanel.vcf.gz

#List of variants in reference panel 
#(CHROM\tPOS\tREF,ALT, no header)
TSV=chr${CHR}.referencePanel.tsv.gz

#Accessible genome mask
MSK=/accessible_genome_masks/chr${CHR}.bed.gz
#what to add to vcf/bcf file header: 
#'##INFO=<ID=MASK,Number=1,Type=String,Description="site passing the 1000 Genomes strict mask">'
HEA=/accessible_genome_masks/header.txt

#USCI repeats regions
RPT=/Repeats/by_chr_sorted/chr${CHR}.bed.gz
#what to add to vcf/bcf file header: 
#'##INFO=<ID=REPEATS,Number=1,Type=String,Description="site found in a UCSC repeat sequence">'
HEA=/Repeats/by_chr_sorted/header.txt

#Output1 VCF and LOG files
SPL1=/calls/chr${CHR}/${SPL}.raw.Q20.q30.calls.spl
OUT=/calls/chr${CHR}/${SPL}.raw.Q20.q30.calls.vcf.gz
TMP=/calls/chr${CHR}/${SPL}.raw.Q20.q30.calls.tmp.vcf.gz

#output2
OUT2=/calls/chr${CHR}/${SPL}.raw.Q20.q30.1000G.mask.noRepeats.calls.vcf.gz

#output3
OUT3=/calls/chr${CHR}/${SPL}.raw.Q20.q30.1000G.mask.noRepeats.qual.DP.calls.vcf.gz

#Call genotypes using bcftools
echo ${BAMDIR}/${BAM} $SPL > $SPL1
bcftools mpileup -f $FASTA -I -E -a 'FORMAT/DP' --ignore-RG -T $VPOS -Q 20 -q 30 -C 50 -r $CHR ${BAMDIR}/${BAM} | bcftools call -Aim -C alleles -T $TSV -Oz -o $OUT
bcftools reheader -s $SPL1 -o $TMP $OUT
mv $TMP $OUT
bcftools index -f $OUT

#1000G mask and remove repeats
bcftools annotate -a $MSK -m +MASK=strict -h $HEA -c CHROM,POS,FROM,TO $OUT | bcftools view -i 'INFO/MASK="strict"' | bcftools annotate -a $RPT -m REPEATS=repeats -h $HEA2 -c CHROM,POS,FROM,TO | bcftools view --exclude 'REPEATS="repeats"' -Oz -o $OUT
bcftools index -f $OUT2

#depth and QUAL filter
DOC=$(bcftools query -f '%INFO/DP\n' $OUT |  awk 'BEGIN { s = 0; l=0; } { s+=$1; l++; } END { print s/l;}')

LOW=$(python3 limitsDoC.py $DOC | awk '{print $1}')
UPP=$(python3 limitsDoC.py $DOC | awk '{print $2}')
echo $LOW $UPP

bcftools filter --exclude "FORMAT/DP<$LOW |  FMT/DP>${UPP} & QUAL<30" $OUT2 -Oz -o $OUT3
bcftools index -f $OUT3

rm $OUT
rm $OUT2

