#!/bin/bash
##
## Build STAR indices for human hg38 + mouse m38 + alien
## Requirements: run human, mouse, alien indexing first (so FASTAs/GTFs are retrieved by them)
##
## Izaskun Mallona
## 3rd Dec 2021


NAMES=(human mouse alien)

declare -A FAS
FAS=([human]="../GRCh38.p13.genome.fa" [mouse]="../GRCm38.p6.genome.fa" [alien]="../alien.fa")

declare -A GTFS
GTFS=([human]="../gencode.v38.basic.annotation.gtf" [mouse]="../gencode.vM25.annotation.gtf" [alien]="../alien.gtf")

NTHREADS=20
ID=human_mouse_alien_combined
STAR=/home/imallona/soft/star/bin/STAR

mkdir -p ~/giulia/indices
cd "$_"

cd ~/giulia/indices
mkdir -p "$ID"
cd "$ID"

## prepend the species/assembly to each
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

