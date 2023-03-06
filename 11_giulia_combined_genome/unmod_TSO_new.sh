#!/bin/bash
##

r1=/home/gmoro/fastqs_seventh_scRNAseq/o307161_1-Unmodified_S4_R1_001.fastq.gz
r2=/home/gmoro/fastqs_seventh_scRNAseq/o307161_1-Unmodified_S4_R2_001.fastq.gz
SRC=/home/imallona/src/ebrunner_spectral

export PATH=/home/imallona/soft/star/STAR-2.7.10b/source:$PATH

#output dir

mkdir /home/gmoro/star_solo_cell_lines/unmod/TSO_new
cd /home/gmoro/star_solo_cell_lines/unmod/TSO_new

nice -n 19 STAR --runThreadN 5 \
     --genomeDir /home/gmoro/indices/correct_human_mouse_alien_WPRE_combined/ \
     --readFilesCommand zcat \
     --outFileNamePrefix unmod_TSO_new/ \
     --readFilesIn "$r2" "$r1" \
     --soloType CB_UMI_Complex \
     --soloAdapterSequence AATGNNNNNNNNNCCAC \
     --soloCBposition 2_-9_2_-1 2_4_2_12 2_17_2_25 \
     --soloUMIposition 3_10_3_17 \
     --soloCBwhitelist "$SRC"/07_barcodes_translation_sbg/data/CLS1.txt "$SRC"/07_barcodes_translation_sbg/data/CLS2.txt "$SRC"/07_barcodes_translation_sbg/data/CLS3.txt \
     --soloCBmatchWLtype 1MM \
     --soloCellFilter None \
     --soloCellReadStats Standard \
     --sjdbOverhang 61 \
     --outSAMattributes CB UB gx gn sS CR CY UR UY\
     --outSAMtype BAM SortedByCoordinate \
     --quantMode GeneCounts \
     --sjdbGTFfile /home/gmoro/indices/combined_correct.gtf
