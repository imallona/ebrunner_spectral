#!/bin/bash
##
## Merges TSO/WTA bamfiles, adds a readgroup to them, updates the CB to add CB+RG
##

FEATURECOUNTS=/home/imallona/soft/subread-2.0.4-source/bin/featureCounts
GENOME_FA=/home/gmoro/indices/combined_correct.fa
NTHREADS=20

echo "samtools merge TSO /WTA"

cd /home/gmoro/featurecounts_cell_line_experiment/RoCK_ROI #directory where we want to save the files
TSO=/home/gmoro/star_solo_cell_lines/RoCK_ROI/TSO_new/RoCK_ROI_TSO_new/Aligned.sortedByCoord.out.bam # paths to starsolo bam files
WTA=/home/gmoro/star_solo_cell_lines/RoCK_ROI/WTA_new/RoCK_ROI_WTA_new/Aligned.sortedByCoord.out.bam #paths to starsolo bam files

## Linking the .bam file and giving it a new name

ln -s $TSO tso.bam #this only needs to be run once otherwise error
ln -s $WTA wta.bam #this only needs to be run once otherwise error

## Checking if .bam files linked

samtools view -h tso.bam  | head -10
samtools view -h wta.bam  | head -10

## test with tiny bam file

#samtools view -h tso.bam  | head -1000 | samtools view -Sb > tiny_tso.bam
#samtools view -h wta.bam | head -1000 | samtools view -Sb > tiny_wta.bam

## Merging the CB with the RG in the MO column, filtering the lines with the CB not empty and deduplicating UMIs

samtools merge -@ $NTHREADS - -r tso.bam wta.bam  > merged.bam
samtools view -H merged.bam  > header_merged.sam

samtools view -h merged.bam | grep -vP "CR:Z:\t" | \
awk '{OFS=FS="\t"; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,$29,"MO"substr($29,3 )"_"substr($27, 6)}' | uniq -f 27 | cat header_merged.sam - |  \
samtools view -Sb -@ $NTHREADS | \
samtools sort -@ $NTHREADS > filtered_merged.bam

rm merged.bam
rm tiny_merged.sam


$FEATURECOUNTS -a /home/gmoro/indices/combined_correct.gtf \
               -o RoCK_ROI \
               filtered_merged.bam \
               -F GTF \
               -t exon \
               -g transcript_id \
               -f \
               -O \
               -M  \
               -J \
               -G "$GENOME_FA" \
               -T "$NTHREADS" \
               --byReadGroup


