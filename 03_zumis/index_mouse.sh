#!/bin/bash
##
## Build STAR indices for mouse mm38
##
## Izaskun Mallona
## 20 Sept 2021

FA=Mus_musculus.GRCm38.dna.primary_assembly.fa
GTF=gencode.vM25.annotation.gtf
NTHREADS=20
ID=GRCm38_gencode_M25

wget http://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/dna/"$FA".gz
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/"$GTF".gz

mkdir -p ~/indices
cd "$_"

## installs
cd ~/indices
  
mkdir -p "$ID"

pigz --decompress -p "$NTHREADS" "$GTF".gz
pigz --decompress -p "$NTHREADS" "$FA".gz

sed -i 's/>/>chr/g' "$FA" 

STAR --runThreadN "$NTHREADS" \
        --runMode genomeGenerate \
        --genomeDir "$ID" \
        --genomeFastaFiles "$FA" \
        --sjdbGTFfile "$GTF" \
        --sjdbOverhang 61

pigz -p "$NTHREADS" "$FA" ;
pigz -p "$NTHREADS" "$GTF"
