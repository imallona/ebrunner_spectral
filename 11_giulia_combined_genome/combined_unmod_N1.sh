#!/bin/bash
##

## parameters

r1=/home/gmoro/fastqs_seventh_scRNAseq/o307161_2-Unmodified_N1_S1_R1_001.fastq.gz
r2=/home/gmoro/fastqs_seventh_scRNAseq/o307161_2-Unmodified_N1_S1_R2_001.fastq.gz
SRC=/home/imallona/src/ebrunner_spectral
SAMPLE=unmod_N1
WD=/home/gmoro/star_solo_cell_lines/featurecounts/"$SAMPLE"

MODE1=WTA
MODE2=TSO
MODE3=featurecounts

FEATURECOUNTS=/home/imallona/soft/subread-2.0.4-source/bin/featureCounts
GENOME_FA=/home/gmoro/indices/combined_correct.fa
NTHREADS=20

export PATH=/home/imallona/soft/star/STAR-2.7.10b/source:$PATH

## star_solo_WTA

#output dir

mkdir "$WD"/"${SAMPLE}${MODE1}"
cd "$WD"/"${SAMPLE}${MODE1}"

nice -n 19 STAR --runThreadN 5 \
     --genomeDir /home/gmoro/indices/correct_human_mouse_alien_WPRE_combined/ \
     --readFilesCommand zcat \
     --outFileNamePrefix "${SAMPLE}${MODE1}"/ \
     --readFilesIn "$r2" "$r1" \
     --soloType CB_UMI_Complex \
     --soloAdapterSequence GTGANNNNNNNNNGACA \
     --soloCBposition 2_-9_2_-1 2_4_2_12 2_17_2_25 \
     --soloUMIposition 3_10_3_17 \
     --soloCellReadStats Standard \
     --sjdbOverhang 61 \
     --soloCBwhitelist "$SRC"/07_barcodes_translation_sbg/data/CLS1.txt "$SRC"/07_barcodes_translation_sbg/data/CLS2.txt "$SRC"/07_barcodes_translation_sbg/data/CLS3.txt \
     --soloCBmatchWLtype 1MM \
     --soloCellFilter EmptyDrops_CR 10000 0.99 10 45000 90000 500 0.01 20000 0.05 10000 \
     --outSAMattributes NH HI AS nM NM MD jM jI MC ch CB UB gx gn sS CR CY UR UY\
     --outSAMtype BAM SortedByCoordinate \
     --quantMode GeneCounts \
     --sjdbGTFfile /home/gmoro/indices/combined_correct.gtf

echo "starsolo WTA finished"





## star_solo_TSO

#output dir

mkdir "$WD"/"${SAMPLE}${MODE2}"
cd "$WD"/"${SAMPLE}${MODE2}"

nice -n 19 STAR --runThreadN 5 \
     --genomeDir /home/gmoro/indices/correct_human_mouse_alien_WPRE_combined/ \
     --readFilesCommand zcat \
     --outFileNamePrefix "${SAMPLE}${MODE2}"/ \
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
     --outSAMattributes NH HI AS nM NM MD jM jI MC ch CB UB gx gn sS CR CY UR UY\
     --outSAMtype BAM SortedByCoordinate \
     --quantMode GeneCounts \
     --sjdbGTFfile /home/gmoro/indices/combined_correct.gtf

echo "starsolo TSO finished"






## featurecounts alien sequences

TSO="$WD"/"${SAMPLE}${MODE2}"/"${SAMPLE}${MODE2}"/Aligned.sortedByCoord.out.bam # paths to starsolo bam files
WTA="$WD"/"${SAMPLE}${MODE1}"/"${SAMPLE}${MODE1}"/Aligned.sortedByCoord.out.bam #paths to starsolo bam files

mkdir "$WD"/"${SAMPLE}${MODE3}"
cd "$WD"/"${SAMPLE}${MODE3}" #directory where we want to save the files

## Linking the .bam file and giving it a new name

ln -s $TSO tso.bam #this only needs to be run once otherwise error
ln -s $WTA wta.bam #this only needs to be run once otherwise error

## Indexing .bam file

samtools index tso.bam
samtools index wta.bam

## Extracting alien chromosomes

samtools view -h tso.bam alien_egfp_WPRE alien_tdtomato_WPRE | samtools view -Sb > alien_tso.bam
samtools view -h wta.bam alien_egfp_WPRE alien_tdtomato_WPRE | samtools view -Sb > alien_wta.bam

## Merging .bam file

samtools merge -@ $NTHREADS - -r alien_merged.bam alien_tso.bam alien_wta.bam > alien_merged.bam
echo "finished merging"

## Defining parameters

TAB="$(printf '\t')" #defining what tab is
valid_rgs="$WD"/"${SAMPLE}${MODE1}"/"${SAMPLE}${MODE1}"/Solo.out/Gene/filtered/barcodes.tsv #these are the barcodes defined by starsolo, need to change every time

## Extracing wta barcodes

while read -r barcode
do
    echo "@RG${TAB}ID:alien_tso_${barcode}"
    echo "@RG${TAB}ID:alien_wta_${barcode}"

done < "$valid_rgs" > header_barcodes

sed 's/@RG\tID://g' header_barcodes > wta_tso_filtered_barcodes #extracting the barcodes from the header

echo "finished extracting barcodes"

## Extracting header

samtools view -H alien_merged.bam > curr_header.txt

## Selecting CBs with barcode and adding @RG to header

samtools view alien_merged.bam | \
    grep -vP "CB:Z:-\t" | \
    awk '{OFS=FS="\t"; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,"RG"substr($29,3 )"_"substr($27, 6)}' |  uniq -f 27 | \
    cat curr_header.txt header_barcodes - | \
    samtools view -Sb -@ $NTHREADS | \
    samtools sort -@ $NTHREADS  > alien_merged_header.bam

echo "finished adding RG"

## Indexing new .bam file

samtools index -@ $NTHREADS alien_merged_header.bam

echo "finished indexing new file"

## Filtering the barcodes extracted before in the header

samtools view -h alien_merged_header.bam -R wta_tso_filtered_barcodes | \
    samtools view -Sb >  alien_merged_filtered.bam

echo "finished filtering barcodes"

## Getting .gtf file

grep alien /home/gmoro/indices/combined_correct.gtf > alien.gtf

## Removing files we don't need

rm alien_tso.bam
rm alien_wta.bam
rm alien_merged.bam
rm header_barcodes
rm wta_tso_filtered_barcodes
rm curr_header.txt
rm alien_merged_header.bam

# Featurecounts

$FEATURECOUNTS -a alien.gtf \
               -o alien_only \
               alien_merged_filtered.bam \
               -F GTF \
               -t exon \
               -g gene_id \
               -f \
               -O \
               -M  \
               -T "$NTHREADS" \
               --byReadGroup #\
               #-J \
               # -G "$GENOME_FA"

