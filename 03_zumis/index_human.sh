#!/bin/bash
##
## Build STAR indices for human hg38
##
## Izaskun Mallona
## 20 Sept 2021

FA=GRCh38.p13.genome.fa
GTF=gencode.v38.basic.annotation.gtf
NTHREADS=20
ID=GRCh38_gencode_38
STAR=/home/imallona/soft/star/bin/STAR

mkdir -p ~/giulia/indices
cd "$_"

cd ~/giulia/indices
mkdir -p "$ID"


wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/"$FA".gz
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/"$GTF".gz

## installs

pigz --decompress -p "$NTHREADS" "$GTF".gz
pigz --decompress -p "$NTHREADS" "$FA".gz

$STAR --runThreadN "$NTHREADS" \
        --runMode genomeGenerate \
        --genomeDir "$ID" \
        --genomeFastaFiles "$FA" \
        --sjdbGTFfile "$GTF" \
        --sjdbOverhang 61

pigz -p "$NTHREADS" "$FA"
pigz -p "$NTHREADS" "$GTF"
