#!/bin/bash
##
## Build STAR indices for tdTomato and eGFP
##
## Izaskun Mallona
## 20 Sept 2021

cd /home/gmoro/indices
STAR=/home/imallona/soft/star/bin/STAR
GENOME=alien_WPRE.fa
GTF=alien_genes_WPRE.gtf
NTHREADS=10
ID=alien_WPRE

## installs
cd /home/gmoro/indices

mkdir -p "$ID"

# $STAR --runThreadN "$NTHREADS" \
#         --runMode genomeGenerate \
#         --genomeDir "$ID" \
#         --genomeFastaFiles "$GENOME" \
#         --sjdbGTFfile "$GTF" \
#         --sjdbOverhang 59 \
#         --genomeSAindexNbases 7

$STAR --runMode genomeGenerate \
      --runThreadN "$NTHREADS" \
      --genomeDir "$ID" \
      --genomeFastaFiles "$GENOME" \
      --genomeSAindexNbases 7

pigz -p "$NTHREADS" "$GENOME"
pigz -p "$NTHREADS" "$GTF"

