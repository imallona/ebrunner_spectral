#!/bin/bash
##
## Build STAR indices for mouse mm38
##
## Izaskun Mallona
## 20 Sept 2021

FA=GRCm38.p6.genome.fa
GTF=gencode.vM25.annotation.gtf
NTHREADS=20
ID=GRCm38_gencode_M25


mkdir -p ~/giulia/indices
cd "$_"


wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/"$FA".gz
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/"$GTF".gz


cd ~/giulia/indices
  
mkdir -p "$ID"

pigz --decompress -p "$NTHREADS" "$GTF".gz
pigz --decompress -p "$NTHREADS" "$FA".gz

STAR --runThreadN "$NTHREADS" \
        --runMode genomeGenerate \
        --genomeDir "$ID" \
        --genomeFastaFiles "$FA" \
        --sjdbGTFfile "$GTF" \
        --sjdbOverhang 61

pigz -p "$NTHREADS" "$FA" ;
pigz -p "$NTHREADS" "$GTF"
