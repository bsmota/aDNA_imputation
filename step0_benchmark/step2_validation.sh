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


#Just SNPs for now
VTYPE=snps_only
REF=1000GP_nygc_umich
#Position file	

MSK=/work/FAC/FBM/DBC/amalaspi/popgen/shared_ressources/accessible_genome_masks/diana_masks/chr${CHR}.bed.gz
HEA=/work/FAC/FBM/DBC/amalaspi/popgen/shared_ressources/accessible_genome_masks/diana_masks/header.txt


OUT2=${PDIR1}/chr${CHR}/${SPL}.raw.Q20.q30.1000G.mask.calls.vcf.gz

mkdir -p ${PDIR1}/chr${CHR}/


#Output VCF and LOG files
SPL1=${PDIR1}/chr${CHR}/${SPL}.raw.Q20.q30.calls.spl
OUT=${PDIR1}/chr${CHR}/${SPL}.raw.Q20.q30.calls.vcf.gz
TMP=${PDIR1}/chr${CHR}/${SPL}.raw.Q20.q30.calls.tmp.vcf.gz
LOG=${PDIR1}/chr${CHR}/${SPL}.raw.Q20.q30.calls.log


#DOC=$(bcftools query -f '%INFO/DP\n' $OUT |  awk 'BEGIN { s = 0; l=0; } { s+=$1; l++; } END { print s/l;}')

echo $DOC; 

#LOW=$(python3 limitsDoC.py $DOC | awk '{print $1}')
#UPP=$(python3 limitsDoC.py $DOC | awk '{print $2}')
echo $LOW $UPP


#Call genotypes using bcftools
echo ${BAMDIR}/${BAM} $SPL > $SPL1
#bcftools mpileup -f $FASTA -I -E -a 'FORMAT/DP' --ignore-RG -T $VPOS -Q 20 -q 30 -C 50 -r $CHR ${BAMDIR}/${BAM} | bcftools call -Aim -C alleles -T $TSV -Oz -o $OUT
#bcftools reheader -s $SPL1 -o $TMP $OUT
mv $TMP $OUT
#bcftools index -f $OUT
echo $CHR $BAM $? >> $LOG


#1000G mask
#bcftools annotate -a $MSK -m MASK=strict -h $HEA -c CHROM,POS,FROM,TO $OUT | bcftools view -i 'INFO/MASK="strict"' -Oz -o $OUT2
#bcftools index -f $OUT2

#remove repeats, filter for DP and QUAL
RPT=/work/FAC/FBM/DBC/amalaspi/popgen/shared_ressources/Repeats/by_chr_sorted/chr${CHR}.bed.gz
HEA=/work/FAC/FBM/DBC/amalaspi/popgen/shared_ressources/Repeats/by_chr_sorted/header.txt

OUT3=${PDIR1}/chr${CHR}/${SPL}.raw.Q20.q30.1000G.mask.noRepeats.qual.calls.vcf.gz

bcftools annotate -a $RPT -m REPEATS=repeats -h $HEA -c CHROM,POS,FROM,TO $OUT2 | bcftools view --exclude 'REPEATS="repeats"' | bcftools filter --exclude "QUAL<30" -Oz -o $OUT3
bcftools index -f $OUT3
