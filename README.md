# ALU_insertion
Usage example:

BAM=~/GITHUB/ALLELE_INSERTION/NGS625_48_333370_QI_F_VCP1R134Via_Pan4119_S48_R1_001.bam
SAMPLE_ID=$(basename $BAM | grep -oP 'NGS[^_]+_\d+')

docker run --rm \
  -v ~/GITHUB/ALLELE_INSERTION/NGS625_48_333370_QI_F_VCP1R134Via_Pan4119_S48_R1_001.bam:/app/data/${SAMPLE_ID}.bam \
  -v ~/GITHUB/ALLELE_INSERTION/NGS625_48_333370_QI_F_VCP1R134Via_Pan4119_S48_R1_001.bai:/app/data/${SAMPLE_ID}.bai \
  -v ~/GITHUB/ALLELE_INSERTION/hs37d5.fa:/app/data/reference.fa \
  -v ~/GITHUB/ALLELE_INSERTION/hs37d5.fa.fai:/app/data/reference.fa.fai \
  -v ~/GITHUB/ALLELE_INSERTION/hs37d5.fa.nhr:/app/data/reference.fa.nhr \
  -v ~/GITHUB/ALLELE_INSERTION/hs37d5.fa.nin:/app/data/reference.fa.nin \
  -v ~/GITHUB/ALLELE_INSERTION/hs37d5.fa.nsq:/app/data/reference.fa.nsq \
  -v $(pwd)/output:/app/output \
  seglh/scramble:latest \
  --bam /app/data/${SAMPLE_ID}.bam \ 
  --bai /app/data/${SAMPLE_ID}.bai

Description of docker run options:
--bam /app/data/NGS625_48.bam        # required: path to BAM file
--bai /app/data/NGS625_48.bai        # required: path to BAI file
--window 100                          # default: 50  — coverage analysis window size
--polyA_window 20                     # default: 10  — polyA detection window size
--threshold 15.0                      # default: 10.0 — coverage change threshold (%)
--min_polyA_len 8                     # default: 10  — minimum polyA/T stretch length
--merge_gap 5                         # default: 1   — max gap for merging nearby variants
--proximity 20                        # default: 10  — max distance from variant to coverage tract
--bed /app/data/my_regions.bed        # default: None — restrict analysis to bed file regions
--verbose                             # default: off — enable detailed logging
