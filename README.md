# Variant-calling-on-tumor-tissue
1. Download data
   - The Human reference genome hg38 was downloaded from https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
   - The genome of the tumor tissue genome was downloaded from https://zenodo.org/record/2582555/files/SLGFSK-T_231336_r1_chr5_12_17.fastq.gz

2. Evaluate quality of raw data with FastQC
   - Discovered that the raw data had bionomial distrubution of GC content

3. Trim raw reads with Trimmomatic

4. Align the reads with BWA-MEM
   - The reads were aligned with BWA tools, human reference genome 38 (hg38) was used
  
5. Mark Duplicates and sort with GATK4
   - This tool locates and tags duplicate reads in a BAM or SAM file, where duplicate reads are defined as originating from a single fragment of DNA
  
6. Recalibrate base quality
   - the BaseRecalibrator tool builds a model of covariation based on the input data and a set of known variants, producing a recalibration file; then the ApplyBQSR tool adjusts the base quality scores in the data based on the model, producing a new BAM file.

7. Call variants with GATK4

8. Filtered SNPs with the following filter
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

9. Filtered INDELs with the following filter
    -filter-name "QD_filter" -filter "QD < 2.0" \
    -filter-name "FS-filter" -filter "FS > 200.0" \
    -filter-name "SOR_filter" -filter "SOR > 10.0" \
    -genotype-filter-expression "DP < 10" \
    -genotype-filter-name "DP_filter" \
    -genotype-filter-expression "GQ < 10" \
    -genotype-filter-name "GQ_filter"

10. Annotated SNPs and INDELs that pass the filter and recorded the final results in a table
