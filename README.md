# ALU_insertion
## Build instructions

```
make build

# or, in cases where you want to remove the cache and create an entirely new docker image
make cleanbuild
```
## Usage example:

```
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
```

## Docker additional run options:
|Option|Default value|Description|
|------|-------------|-----------|
|--window|50|coverage analysis window size|
|--polyA_window|10|polyA detection window size|
|--threshold|10.0|coverage change threshold (%)|
|--min_polyA_len|10|minimum polyA/T stretch length|
|--merge_gap|1|max gap for merging nearby variants|
|--bed|None|restrict analysis to bed file regions|
|--verbose|off|enable detailed logging|


## To test docker locally
make cleanbuild | bash ./test_command.sh
make build | bash ./test_command.sh

## To test scramble_filtering_vcf_updated_v2.py locally - replace test bam path
python scramble_filtering_vcf_updated_v2.py \
  --vcf "output/NGS625_48_ALU_ins.vcf" \
  --bam "/home/isabeljohnsondavies/GITHUB/ALLELE_INSERTION/NGS625_48_333370_QI_F_VCP1R134Via_Pan4119_S48_R1_001.bam" \ 
  --id "$SAMPLE_ID"