#!/bin/bash
##

r1=/home/gmoro/fastqs_PDGFRA/single_clean/o303001_1-Unmodified_S1_R1_001.fastq.gz
r2=/home/gmoro/fastqs_PDGFRA/single_clean/o303001_1-Unmodified_S1_R2_001.fastq.gz
SRC=/home/imallona/src/ebrunner_spectral

export PATH=/home/imallona/soft/star/STAR-2.7.10b/source:$PATH

#output dir
rm -r /home/gmoro/star_solo_PDGFRA/unmod_single_clean/unmod_single_clean_wta_output
mkdir /home/gmoro/star_solo_PDGFRA/unmod_single_clean/unmod_single_clean_wta_output
cd /home/gmoro/star_solo_PDGFRA/unmod_single_clean/unmod_single_clean_wta_output

nice -n 19 STAR --runThreadN 5 \
     --genomeDir /home/imallona/giulia/indices/GRCm38_gencode_M25_gfp_star27b/ \
     --readFilesCommand zcat \
     --outFileNamePrefix unmod_single_clean_wta/ \
     --readFilesIn "$r2" "$r1" \
     --soloType CB_UMI_Complex \
     --soloAdapterSequence GTGANNNNNNNNNGACA \
     --soloCBposition 2_-9_2_-1 2_4_2_12 2_17_2_25 \
     --soloUMIposition 3_10_3_17 \
     --soloCBwhitelist "$SRC"/07_barcodes_translation_sbg/data/CLS1.txt "$SRC"/07_barcodes_translation_sbg/data/CLS2.txt "$SRC"/07_barcodes_translation_sbg/data/CLS3.txt \
     --soloCBmatchWLtype 1MM \
     --soloCellFilter EmptyDrops_CR \
     --outSAMattributes CB UB gx gn sS CR CY UR UY\
     --outSAMtype BAM SortedByCoordinate \
     --quantMode GeneCounts \
     --sjdbGTFfile /home/imallona/giulia/indices/combined.gtf



