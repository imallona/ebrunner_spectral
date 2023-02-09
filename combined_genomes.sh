#!/bin/bash
##
## Build STAR indices for human hg38 + mouse m38 + alien_WPRE
## Requirements: run human, mouse, alien indexing first (so FASTAs/GTFs are retrieved by them)
##
## Izaskun Mallona
## 3rd Dec 2021


NAMES=(human mouse alien_WPRE)

declare -A FAS
FAS=([human]="../GRCh38.p13.genome.fa" [mouse]="../GRCm38.p6.genome.fa" [alien]="../alien_WPRE.fa")

declare -A GTFS
GTFS=([human]="../gencode.v38.basic.annotation.gtf" [mouse]="../gencode.vM25.annotation.gtf" [alien]="../alien_genes_WPRE.gtf")

NTHREADS=20
ID=human_mouse_alien_WPRE_combined
STAR=/home/imallona/soft/star/bin/STAR

mkdir -p /home/gmoro/indices
cd "$_"

cd /home/gmoro/indices
mkdir -p "$ID"
cd "$ID"

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

$STAR --runMode genomeGenerate \
      --runThreadN "$NTHREADS" \
      --genomeDir "$ID" \
      --genomeFastaFiles combined.fa

pigz -p "$NTHREADS" combined.fa
pigz -p "$NTHREADS" combined.gtf

