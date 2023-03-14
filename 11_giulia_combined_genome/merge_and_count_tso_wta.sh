#!/bin/bash
##
## Merges TSO/WTA bamfiles, adds a readgroup to them, updates the CB to add CB+RG
##
##  6th March 2023

## featurecounts some features, to tackle the multimapping problem

## subread install
## from https://sourceforge.net/projects/subread/
# mkdir -p /home/imallona/soft/subread
# tar xzvf *tar.gz
# cd subread-2.0.4-source/src
# make -f Makefile.Linux

FEATURECOUNTS=/home/imallona/soft/subread-2.0.4-source/bin/featureCounts
WD=/home/imallona/another_test
GENOME_FA=/home/gmoro/indices/combined_correct.fa
NTHREADS=20

mkdir -p $WD
cd $WD

## Giulia does not give read access to anyone else (?)

sudo chown -R gmoro:robinsonlab /home/gmoro/star_solo_cell_lines/RoCK/TSO_new/RoCK_TSO_new/
sudo chmod -R g+rx /home/gmoro/star_solo_cell_lines/RoCK/TSO_new/RoCK_TSO_new/
ln -s /home/gmoro/star_solo_cell_lines/RoCK/TSO_new/RoCK_TSO_new/Aligned.sortedByCoord.out.bam tso.bam

sudo chown -R gmoro:robinsonlab /home/gmoro/star_solo_cell_lines/RoCK/WTA_new/RoCK_WTA_new/
sudo chmod -R g+rx /home/gmoro/star_solo_cell_lines/RoCK/WTA_new/RoCK_WTA_new/
ln -s /home/gmoro/star_solo_cell_lines/RoCK/WTA_new/RoCK_WTA_new/Aligned.sortedByCoord.out.bam wta.bam

## index
samtools index tso.bam
samtools index wta.bam

## extract alien features, only
samtools view -h tso.bam alien_egfp_WPRE alien_tdtomato_WPRE | samtools view -Sb > small_tso.bam
samtools view -h wta.bam alien_egfp_WPRE alien_tdtomato_WPRE | samtools view -Sb > small_wta.bam

## merge and add RGs (RG are bam filenames)
samtools merge -@ $NTHREADS - -r small_tso.bam small_wta.bam  > merged_temp.bam

## these are the filtered-in WTA CBs
valid_rgs=/home/gmoro/star_solo_cell_lines/RoCK/WTA_new/RoCK_WTA_new/Solo.out/Gene/filtered/barcodes.tsv
TAB="$(printf '\t')"

## this generate a SAM header with RG entries from the 'filtered in' barcodes for both TSO and WTA
while read -r barcode
do
    
    echo "@RG${TAB}ID:small_tso_${barcode}"
    echo "@RG${TAB}ID:small_wta_${barcode}"

done < "$valid_rgs" > header_barcodes

## this strips the RG flags and keeps the filtered in barcodes, plain
sed 's/@RG\tID://g' header_barcodes > wta_tso_filteredin_barcodes

## extract current bam header
samtools view -H merged_temp.bam > curr_header.txt

## rename CB to RG, reheader, sort
samtools view merged_temp.bam | \
    grep -vP "CB:Z:-\t" | \
    awk '{OFS=FS="\t"; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,"RG"substr($29,3 )"_"substr($27, 6)}' |  uniq -f 27 | \
    cat curr_header.txt header_barcodes - | \
    samtools view -Sb -@ $NTHREADS | \
    samtools sort -@ $NTHREADS  > merged_untested_header.bam


rm merged_temp.bam

## index the reheaded bam
samtools index -@ $NTHREADS merged_untested_header.bam

## filter in only those CBs (RGs) that were filtered in in WTA
samtools view -h merged_untested_header.bam -R wta_tso_filteredin_barcodes | \
    samtools view -Sb >  merge_untested_filtered.bam

## get a small GTF, alien-only
grep alien /home/gmoro/indices/combined_correct.gtf > alien.gtf

## count multimapping, multioverlap included
$FEATURECOUNTS -a alien.gtf \
               -o alien_only \
               merged_untested_header_filtered.bam \
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

rm merged_untested_header.bam curr_header.txt wta_tso_filteredin_barcodes
