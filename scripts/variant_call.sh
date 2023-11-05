#!/bin/bash

#variant calling steps
#directories
ref="/home/bryan/variant_call/supporting_files/hg38/hg38.fa"
known_sites="/home/bryan/variant_call/supporting_files/hg38/Homo_sapiens_assembly38.dbsnp138.vcf"
aligned_reads="/home/bryan/variant_call/aligned_reads"
reads="/home/bryan/variant_call/reads"
results="/home/bryan/variant_call/results"
data="/home/bryan/variant_call/data"
trim="/home/bryan/variant_call/supporting_files/TruSeq3-PE.fa"

#STEP 1- QC Run with fastqc
echo "STEP 1- Running QC with fastqc"

fastqc ${reads}/SLGFSK-T_231336_r1_chr5_12_17.fastq.gz -o ${reads}/
fastqc ${reads}/SLGFSK-T_231336_r2_chr5_12_17.fastq.gz -o ${reads}/

#STEP 2- Trim Reads with Trimmomatic
echo "STEP 2- Trim Reads"

java -jar ~/apps/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 8 ${reads}/SLGFSK-T_231336_r1_chr5_12_17.fastq.gz ${reads}/SLGFSK-T_231336_r2_chr5_12_17.fastq.gz \
    trimmed.paired.SLGFSK-T_231336_r1_chr5_12_17.fastq.gz trimmed.unpaired.SLGFSK-T_231336_r1_chr5_12_17.fastq.gz \
    trimmed.paired.SLGFSK-T_231336_r2_chr5_12_17.fastq.gz trimmed.unpaired.SLGFSK-T_231336_r2_chr5_12_17.fastq.gz \
    ILLUMINACLIP:${trim}:2:30:10:8:true HEADCROP:3 TRAILING:10 MINLEN:25

#STEP 3- Rerun QC on trimmed reads
echo "STEP 3- Rerun QC"
fastqc ${reads}/trimmed.paired.SLGFSK-T_231336_r1_chr5_12_17.fastq.gz -o ${reads}/
fastqc ${reads}/trimmed.paired.SLGFSK-T_231336_r2_chr5_12_17.fastq.gz -o ${reads}/
fastqc ${reads}/trimmed.unpaired.SLGFSK-T_231336_r1_chr5_12_17.fastq.gz -o ${reads}/
fastqc ${reads}/trimmed.unpaired.SLGFSK-T_231336_r2_chr5_12_17.fastq.gz -o ${reads}/

#STEP 4- Map to reference using BWA-MEM
echo "STEP 4- Map to reference using BWA-MEM"

#BWA index reference
bwa index ${ref}

#BWA alignment
bwa mem -t 4 -R "@RG\tID:SLGFSK-T_231336\tPL:ILLUMINA\tSM:SLGFSK-T_231336" ${ref} ${reads}/trimmed.paired.SLGFSK-T_231336_r1_chr5_12_17.fastq.gz ${reads}/trimmed.paired.SLGFSK-T_231336_r2_chr5_12_17.fastq.gz >${aligned_reads}/SLGFSK-T_231336.paired.sam

#STEP 5- Mark Duplicates and Sort
echo "STEP 5- Mark Duplicates and Sort with GATK4"

gatk MarkDuplicatesSpark -I ${aligned_reads}/SLGFSK-T_231336.paired.sam -O ${aligned_reads}/SLGFSK-T_231336_sorted_dedup_reads.bam

#STEP 6- Base quality recalibration
echo "STEP 6- Base quality recalibration"

#Build Model
gatk BaseRecalibrator -I ${aligned_reads}/SLGFSK-T_231336_sorted_dedup_reads.bam -R ${ref} --known-sites ${known_sites} -O ${data}/recal_data.table

#Apply the model to adjust the base quality score
gatk ApplyBQSR -I ${aligned_reads}/SLGFSK-T_231336_sorted_dedup_reads.bam -R ${ref} --bqsr-recal-file ${data}/recal_data.table -O ${aligned_reads}/SLGFSK-T_231336_sorted_dedup_bqsr_reads.bam

#STEP 7- Collect Alignment and Insert size Metrics
echo "STEP 7- Collect Alignment and Insert size Metrics"

gatk CollectAlignmentSummaryMetrics -R ${ref} -I ${aligned_reads}/SLGFSK-T_231336_sorted_dedup_bqsr_reads.bam -O ${aligned_reads}/alignment_metrics.txt
gatk CollectInsertSizeMetrics -INPUT ${aligned_reads}/SLGFSK-T_231336_sorted_dedup_bqsr_reads.bam -OUTPUT ${aligned_reads}/insert_size_metrics.txt -Histogram_FILE ${aligned_reads}/insert_size_metrics.pdf

#STEP 8- Call Variants - gatk haplotype caller
echo "STEP 8- Call Variants - gatk haplotype caller"

gatk HaplotypeCaller -R ${ref} -I ${aligned_reads}/SLGFSK-T_231336_sorted_dedup_bqsr_reads.bam -O ${results}/raw_variants.vcf

#Extract SNPs and INDELs
gatk SelectVariants -R ${ref} -V ${results}/raw_variants.vcf --select-type SNP -O ${results}/raw_snps.vcf
gatk SelectVariants -R ${ref} -V ${results}/raw_variants.vcf --select-type INDEL -O ${results}/raw_indels.vcf