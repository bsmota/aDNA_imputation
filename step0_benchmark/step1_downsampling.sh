#!/bin/bash
#SBATCH --job-name=downsampling
#SBATCH --array 1-N
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=16:00:00
#SBATCH --output=downsampling_%a.out
#SBATCH --error=downsampling_%a.er

#downsample N genomes to coverages in the range 0.1x-2.0x
#file with two columns, sample id and bam filename
LST=sampleInfo.txt

#extract sample id
SPL=$(cat $LST | head -n ${SLURM_ARRAY_ID} | tail -n 1 | awk '{print $1}')
BAM=$(cat $LST | head -n ${SLURM_ARRAY_ID} | tail -n 1 | awk '{print $2}')

#directory where bam files are located
BAMDIR=/bams/

#extract depth per position across all chromosomes for original genome
:> ${SPL}.cov.gz
for CHR in {1..22}; do
	bcftools query -f '%INFO/DP\n' /calls/chr${CHR}/${SPL}.raw.calls.vcf.gz | gzip -c >> /calls/${SPL}.cov.gz
done


echo ${BAMDIR}/${BAM}

#Downsample across different coverages
for COV in 0.1 0.25 0.5 0.75 1.0 2.0; do
	mkdir -p /calls/cov${COV}/bams

	DS_BAM=/calls/cov${COV}/bams/${SPL}.bam

	#Compute fraction of reads to sample from INFO/DP
	FRAC=$(zcat ${SPL}.cov.gz | awk -v c=$COV 'BEGIN { s = 0; l=0; } { s+=$1; l++; } END { print 1+c*l/s; }')	
	echo $FRAC > /calls/cov${COV}/bams/${SPL}.frac
	echo "Fraction=" $FRAC

	#1. DOWNSAMPLING SEQUENCING READS
        echo "Downsampling Reads Coverage = " $COV "x ********************************"
       	echo "DS_BAM = " $DS_BAM  

	samtools view -T $FASTA -s $FRAC -bo $DS_BAM ${BAMDIR}/${BAM}
       	samtools index $DS_BAM

done

