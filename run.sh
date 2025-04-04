#!/bin/bash

export REFGENOME=$1 # The path to the reference hg38 reference genome
export TBTFa=$2 # The path to the twoBitToFa executable
export SAMPLES=$3 # A comma separated list of sample names as they are labelled in the VCF
export VCF=$4 # A VCF with all known variants in the sample
export FASTQ1=$5 # Fastq file 1 of paired end sequencing
export FASTQ2=$6 # Fastq file 2 of paired end sequencing
export SAMPID=$7 # Name for the processed data

# Create directories
mkdir Genome
mkdir VCF
for i in ${SAMPLES//,/ }; do
  mkdir ${i}
done

# Create personal genome
for ((c=1;c<23;c+=1)); do
  python write_genome_and_all_vcf.py -s ${SAMPLES} -c ${c} -v ${VCF} --twoBitToFa ${TBTFa} --twobit ${REFGENOME}
  cat Genome/chr${c}.fasta >> Genome/ref_genome.fasta
  cat VCF/chr${c}.vcf >> VCF/all.vcf
done
gzip VCF/all.vcf

bowtie2-build --threads 32 -p Genome/ref_genome.fasta Genome/ref_genome.bowtie

# Align fastqs
bowtie2 -N --no-mixed --no-discordant --dovetail -X 2000 -t -p 32 -x Genome/ref_genome.bowtie -1 ${FASTQ1} -2 ${FASTQ2} -S ${SAMPID}.sam

# Convert to GRCh38 coordinates
python convert_to_hg38.py -i ./${SAMPID}.sam -o ${SAMPID}.processed.sam
samtools view -bS ${SAMPID}.processed.sam > ${SAMPID}.processed.bam
samtools sort -O bam -T ${SAMPID}.temp ${SAMPID}.processed.bam > ${SAMPID}.processed.sorted.bam
samtools index ${SAMPID}.processed.sorted.bam

# Count alleles
python count_alleles.py -s ${SAMPID}.processed.sorted.bam -o ${SAMPID}.vcf -v VCF/all.vcf.gz
