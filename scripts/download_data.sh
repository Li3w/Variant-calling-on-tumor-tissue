#!/bin/bash

#download tumor data 
#wget -c -O ${fastqdir}/normal_r1.fastq.gz https://zenodo.org/record/2582555/files/SLGFSK-N_231335_r1_chr5_12_17.fastq.gz
#wget -c -O  ${fastqdir}/normal_r2.fastq.gz https://zenodo.org/record/2582555/files/SLGFSK-N_231335_r2_chr5_12_17.fastq.gz
wget -P /home/bryan/variant_call/reads https://zenodo.org/record/2582555/files/SLGFSK-T_231336_r1_chr5_12_17.fastq.gz
wget -P /home/bryan/variant_call/reads https://zenodo.org/record/2582555/files/SLGFSK-T_231336_r2_chr5_12_17.fastq.gz



