#!/bin/bash
##
## Build STAR indices for human hg38 + mouse m38 + alien_WPRE
## Requirements: run human, mouse, alien indexing first (so FASTAs/GTFs are retrieved by them)
##
## Izaskun Mallona
## 3rd Dec 2021


NTHREADS=20
export PATH=/home/imallona/soft/star/STAR-2.7.10b/source:$PATH


ID=human_mouse_alien_WPRE_combined_test

mkdir -p /home/imallona/indices
cd /home/imallona/indices

mkdir -p "$ID"
cd "$ID"


NAMES=(human mouse alien)

declare -A FAS
FAS=([human]="/home/gmoro/indices/GRCh38.p13.genome.fa" [mouse]="/home/gmoro/indices//GRCm38.p6.genome.fa" [alien]="/home/gmoro/indices//alien_WPRE.fa")

declare -A GTFS
GTFS=([human]="/home/gmoro/indices//gencode.v38.basic.annotation.gtf" [mouse]="/home/gmoro/indices//gencode.vM25.annotation.gtf" [alien]="/home/gmoro/indices//alien_genes_WPRE.gtf")



# check files exist/their sizes
for ref in ${NAMES[@]}
do
    echo "$ref"
    
    echo ${GTFS[$ref]}
    ls -lh ${GTFS[$ref]}

    echo ${FAS[$ref]}
    ls -lh ${FAS[$ref]}

    echo "-----------------"
done

## prepend the species/assembly to each scaffold/chromosome

for ref in ${NAMES[@]}
do
    echo "$ref"
    sed "s/^>/>${ref}_/g" ${FAS[$ref]} > $(basename ${FAS[$ref]}.prepend)

    grep -v "#" ${GTFS[$ref]} | sed "s/^/${ref}_/g" > $(basename ${GTFS[$ref]}.prepend)
done

cat *fa.prepend > combined.fa
cat *gtf.prepend > combined.gtf

rm *prepend

STAR --runMode genomeGenerate \
      --runThreadN "$NTHREADS" \
      --genomeDir "$ID" \
      --genomeFastaFiles combined.fa

pigz -p "$NTHREADS" combined.fa
pigz -p "$NTHREADS" combined.gtf

