#!/bin/bash

#variant calling steps
#directories
ref="/home/bryan/variant_call/supporting_files/hg38/hg38.fa"
results="/home/bryan/variant_call/results"

#Filter SNPs
gatk VariantFiltration \
    -R ${ref} \
    -V ${results}/raw_snps.vcf \
    -O ${results}/filtered_snps.vcf \
    -filter-name "QD_filter" -filter "QD < 2.0" \
    -filter-name "FS-filter" -filter "FS > 60.0" \
    -filter-name "MQ_filter" -filter "MQ < 40.0" \
    -filter-name "SOR_filter" -filter "SOR > 4.0" \
    -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
    -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"\
    -genotype-filter-expression "DP < 10" \
    -genotype-filter-name "DP_filter" \
    -genotype-filter-expression "GQ < 10" \
    -genotype-filter-name "GQ_filter"

#Filter INDELs
gatk VariantFiltration \
    -R ${ref} \
    -V ${results}/raw_indels.vcf \
    -O ${results}/filtered_indels.vcf \
    -filter-name "QD_filter" -filter "QD < 2.0" \
    -filter-name "FS-filter" -filter "FS > 200.0" \
    -filter-name "SOR_filter" -filter "SOR > 10.0" \
    -genotype-filter-expression "DP < 10" \
    -genotype-filter-name "DP_filter" \
    -genotype-filter-expression "GQ < 10" \
    -genotype-filter-name "GQ_filter"

#Select variants that PASS filters
gatk SelectVariants \
    --exclude-filtered \
    -V ${results}/filtered_snps.vcf \
    -O ${results}/analysis_ready_snps.vcf

gatk SelectVariants \
    --exclude-filtered \
    -V ${results}/filtered_indels.vcf \
    -O ${results}/analysis_ready_indels.vcf

#To exclude variants that failed genotype filter
cat analysis_ready_snps.vcf |grep -v -E "DP_filter|GQ_filter" > analysis_ready_snps_filteredGT.vcf
cat analysis_ready_indels.vcf |grep -v -E "DP_filter|GQ_filter" > analysis_ready_indels_filteredGT.vcf

#Annotate Variants with Funcotator
gatk Funcotator \
    --variant ${results}/analysis_ready_snps_filteredGT.vcf \
    --reference ${ref} \
    --ref-version hg38 \
    --data-sources-path /home/bryan/funcotator_dataSources.v1.7.20200521s \
    --output ${results}/analysis_ready_snps_filteredGT_funcotated.vcf \
    --output-file-format VCF 

gatk Funcotator \
    --variant ${results}/analysis_ready_indels_filteredGT.vcf \
    --reference ${ref} \
    --ref-version hg38 \
    --data-sources-path /home/bryan/funcotator_dataSources.v1.7.20200521s \
    --output ${results}/analysis_ready_indels_filteredGT_funcotated.vcf \
    --output-file-format VCF

gatk VariantsToTable \
    -V ${results}/analysis_ready_snps_filteredGT_funcotated.vcf -F AC -F AN -F DP -F AF -F FUNCOTATION \
    -O ${results}/output_snps.table

gatk VariantsToTable \
    -V ${results}/analysis_ready_indels_filteredGT_funcotated.vcf -F AC -F AN -F DP -F AF -F FUNCOTATION \
    -O ${results}/output_indels.table
