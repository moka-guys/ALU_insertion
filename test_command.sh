# positive controls
#########################################################################################
# NGS625_48 - confirmed ALU at 19:11238723
BAM=~/GITHUB/ALLELE_INSERTION/NGS625_48_333370_QI_F_VCP1R134Via_Pan4119_S48_R1_001.bam
#BAM=~/GITHUB/ALLELE_INSERTION/NGS683_38_361006_PO_U_VCP1R134Via_Pan4119_S38_R1_001.bam
#BAM=~/GITHUB/ALLELE_INSERTION/NGS683_39_361007_PO_U_VCP1R134Via_Pan4119_S39_R1_001.bam
#BAM=~/Downloads/NGS748_24_370072_DB_F_VCP1R134Via_Pan4119_S24_R1_001.bam
SAMPLE_ID=$(basename $BAM | grep -oP 'NGS[^_]+_\d+')
REF_DIR=~/GITHUB/ALLELE_INSERTION/hs37d5


docker run --rm \
  -v $(realpath $BAM):/app/data/${SAMPLE_ID}.bam \
  -v $(realpath ${BAM%.bam}.bai):/app/data/${SAMPLE_ID}.bai \
  -v $(realpath $REF_DIR/hs37d5.fa):/app/data/reference.fa \
  -v $(realpath $REF_DIR/hs37d5.fa.fai):/app/data/reference.fa.fai \
  -v $(realpath $REF_DIR/hs37d5.fa.nhr):/app/data/reference.fa.nhr \
  -v $(realpath $REF_DIR/hs37d5.fa.nin):/app/data/reference.fa.nin \
  -v $(realpath $REF_DIR/hs37d5.fa.nsq):/app/data/reference.fa.nsq \
  -v $(pwd)/output:/app/output \
  seglh/scramble:latest \
  --bam /app/data/${SAMPLE_ID}.bam \
  --bai /app/data/${SAMPLE_ID}.bai \
  --dx_project_id project-J522GFQ0B3vbkvXv7VzgqJXF 

#########################################################################################

# NGS683_38 - confirmed ALU at 

BAM=~/GITHUB/ALLELE_INSERTION/NGS683_38_361006_PO_U_VCP1R134Via_Pan4119_S38_R1_001.bam
SAMPLE_ID=$(basename $BAM | grep -oP 'NGS[^_]+_\d+')

docker run --rm \
  -v $(realpath $BAM):/app/data/${SAMPLE_ID}.bam \
  -v $(realpath ${BAM%.bam}.bai):/app/data/${SAMPLE_ID}.bai \
  -v $(realpath $REF_DIR/hs37d5.fa):/app/data/reference.fa \
  -v $(realpath $REF_DIR/hs37d5.fa.fai):/app/data/reference.fa.fai \
  -v $(realpath $REF_DIR/hs37d5.fa.nhr):/app/data/reference.fa.nhr \
  -v $(realpath $REF_DIR/hs37d5.fa.nin):/app/data/reference.fa.nin \
  -v $(realpath $REF_DIR/hs37d5.fa.nsq):/app/data/reference.fa.nsq \
  -v $(pwd)/output:/app/output \
  seglh/scramble:latest \
  --bam /app/data/${SAMPLE_ID}.bam \
  --bai /app/data/${SAMPLE_ID}.bai \
  --dx_project_id project-J522GFQ0B3vbkvXv7VzgqJXF 

######################################################################################

# NGS683_39 - confirmed ALU at 

BAM=~/GITHUB/ALLELE_INSERTION/NGS683_39_361007_PO_U_VCP1R134Via_Pan4119_S39_R1_001.bam
SAMPLE_ID=$(basename $BAM | grep -oP 'NGS[^_]+_\d+')

docker run --rm \
  -v $(realpath $BAM):/app/data/${SAMPLE_ID}.bam \
  -v $(realpath ${BAM%.bam}.bai):/app/data/${SAMPLE_ID}.bai \
  -v $(realpath $REF_DIR/hs37d5.fa):/app/data/reference.fa \
  -v $(realpath $REF_DIR/hs37d5.fa.fai):/app/data/reference.fa.fai \
  -v $(realpath $REF_DIR/hs37d5.fa.nhr):/app/data/reference.fa.nhr \
  -v $(realpath $REF_DIR/hs37d5.fa.nin):/app/data/reference.fa.nin \
  -v $(realpath $REF_DIR/hs37d5.fa.nsq):/app/data/reference.fa.nsq \
  -v $(pwd)/output:/app/output \
  seglh/scramble:latest \
  --bam /app/data/${SAMPLE_ID}.bam \
  --bai /app/data/${SAMPLE_ID}.bai \
  --dx_project_id project-J522GFQ0B3vbkvXv7VzgqJXF 

###########################################################################################

# NGS748_24 - test sample

BAM=~/Downloads/NGS748_24_370072_DB_F_VCP1R134Via_Pan4119_S24_R1_001.bam
SAMPLE_ID=$(basename $BAM | grep -oP 'NGS[^_]+_\d+')

docker run --rm \
  -v $(realpath $BAM):/app/data/${SAMPLE_ID}.bam \
  -v $(realpath ${BAM%.bam}.bai):/app/data/${SAMPLE_ID}.bai \
  -v $(realpath $REF_DIR/hs37d5.fa):/app/data/reference.fa \
  -v $(realpath $REF_DIR/hs37d5.fa.fai):/app/data/reference.fa.fai \
  -v $(realpath $REF_DIR/hs37d5.fa.nhr):/app/data/reference.fa.nhr \
  -v $(realpath $REF_DIR/hs37d5.fa.nin):/app/data/reference.fa.nin \
  -v $(realpath $REF_DIR/hs37d5.fa.nsq):/app/data/reference.fa.nsq \
  -v $(pwd)/output:/app/output \
  seglh/scramble:latest \
  --bam /app/data/${SAMPLE_ID}.bam \
  --bai /app/data/${SAMPLE_ID}.bai \
  --dx_project_id project-J522GFQ0B3vbkvXv7VzgqJXF 

# test command for local
# python scramble_filtering_vcf_updated_v2.py \
#     --vcf "output/NGS625_48_ALU_ins.vcf" \
#     --bam "/home/isabeljohnsondavies/GITHUB/ALLELE_INSERTION/NGS625_48_333370_QI_F_VCP1R134Via_Pan4119_S48_R1_001.bam" \
#     --id "$SAMPLE_ID"