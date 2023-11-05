#!/bin/bash

#download reference files
wget -P /home/bryan/variant_call/supporting_files/hg38/ https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip /home/bryan/variant_call/supporting_files/hg38/hg38.fa.gz

#index reference file
samtools faidx /home/bryan/variant_call/supporting_files/hg38/hg38.fa

# ref dictionary
gatk CreateSequenceDictionary -R /home/bryan/variant_call/supporting_files/hg38/hg38.fa -O /home/bryan/variant_call/supporting_files/hg38/hg38.dict

#download known sites files for Base Quality score read Callibration
wget -P /home/bryan/variant_call/supporting_files/hg38/ https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
wget -P /home/bryan/variant_call/supporting_files/hg38/ https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx

#download trim data
wget -c -P /home/bryan/variant_call/supporting_files/ https://raw.githubusercontent.com/timflutre/trimmomatic/master/adapters/TruSeq3-PE.fa