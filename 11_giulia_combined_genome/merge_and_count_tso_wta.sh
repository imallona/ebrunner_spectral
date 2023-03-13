#!/bin/bash
##
## Merges TSO/WTA bamfiles, adds a readgroup to them, updates the CB to add CB+RG
##
##  6th March 2023

## featurecounts some features, to tackle the multimapping problem

## subread install
## from https://sourceforge.net/projects/subread/
mkdir -p /home/imallona/soft/subread

tar xzvf *tar.gz

cd subread-2.0.4-source/src

make -f Makefile.Linux

FEATURECOUNTS=/home/imallona/soft/subread-2.0.4-source/bin/featureCounts
WD=/home/gmoro/star_solo_cell_lines/unmod/TSO_new/unmod_TSO_new/
GENOME_FA=/home/gmoro/indices/combined_correct.fa
NTHREADS=20

chown -R gmoro:robinsonlab "$WD"
chmod -R g+rx "$WD"/..



echo "samtools merge TSO /WTA"

cd /home/imallona/giulia/merging_sandbox
TSO=/home/gmoro/star_solo_cell_lines/RoCK_ROI/TSO_new/RoCK_ROI_TSO_new/Aligned.sortedByCoord.out.bam
WTA=/home/gmoro/star_solo_cell_lines/RoCK_ROI/WTA_new/RoCK_ROI_WTA_new/Aligned.sortedByCoord.out.bam

chown -R gmoro:robinsonlab /home/gmoro/star_solo_cell_lines/
chmod -R g+x /home/gmoro/star_solo_cell_lines/

chmod -R g+r /home/gmoro/star_solo_cell_lines/

ln -s $TSO tso.bam
ln -s $WTA wta.bam


samtools view -h tso.bam  | head -1000 | samtools view -Sb > tiny_tso.bam
samtools view -h wta.bam | head -1000 | samtools view -Sb > tiny_wta.bam

echo "mind we need to update the column indices to get the CR and the RG as RG!"
echo "mind we need to deduplicate umis!"
echo ' for that, just print UMI and CB (corrected) last and then uniq with unix'

# samtools merge -@ $NTHREADS - -r tiny_tso.bam tiny_wta.bam | samtools view -h  | grep -vP "CR:Z:\t" | \
#     awk '{OFS=FS="\t"; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,"MO"substr($21,3 ),$21"_"substr($15, 6), $20}' | uniq -f 20 |\
#     samtools view -Sb -@ $NTHREADS | \
#     samtools sort -@ $NTHREADS  > merged.bam

## mind that the uniq does umi dedup using the new RG (with CB and MO) and the UB
samtools merge -@ $NTHREADS - -r tiny_tso.bam tiny_wta.bam | samtools view -h  | grep -vP "CB:Z:-\t" | \
    awk '{OFS=FS="\t"; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,"MO"substr($21,3 ),$21"_"substr($19, 6), $20}' | uniq -f 20 |\
    samtools view -Sb -@ $NTHREADS | \
    samtools sort -@ $NTHREADS  > merged.bam


echo "think of another bamfile/another flag for non canonical, no-CR containing alignments"
echo untested

sed -i 's/"ROI_3_tdtomato/"ROI_3_tdtomato"/g' test.gtf

tail -100 test.gtf  > tail.gtf


echo "beware the gtf is exon-only, does not have genes nor tx"
$FEATURECOUNTS -a tail.gtf \
               -o test_counting \
               merged.bam \
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
               

## how would it look like for a real set of BAM files?
echo 'todo append to the new scripts'

valid_rgs=~/star_solo_cell_lines/RoCK_ROI/WTA_new/RoCK_ROI_WTA_new/Solo.out/Gene/filtered/barcodes.tsv

# adding RGs

TAB="$(printf '\t')"

# while read -r barcode
# do
    
#     echo "@RG${TAB}ID=tiny_wta_${barcode}"
#     echo "@RG${TAB}ID=tiny_tso_${barcode}"

# done < "$valid_rgs" > header_barcodes

## mind that the uniq does umi dedup using the new RG (with CB and MO) and the UB
samtools merge -@ $NTHREADS - -r tiny_tso.bam tiny_wta.bam  > merged_temp.bam

samtools view -H merged_temp.bam > curr_header.txt

samtools view merged_temp.bam | cut -f 27,29 | grep -vP "CB:Z:-" | awk '{FS=OFS="\t"; print "@RG","ID="substr($2,6)"_"substr($1,6)}' > \
                                                                       header_barcodes

samtools view merged_temp.bam | \
    grep -vP "CB:Z:-\t" | \
    awk '{OFS=FS="\t"; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,"RG"substr($29,3 )"_"substr($27, 6)}' |  uniq -f 27 | \
    cat curr_header.txt header_barcodes - | \
    samtools view -Sb -@ $NTHREADS | \
    samtools sort -@ $NTHREADS  > merged_untested_header.bam


echo "beware the gtf is exon-only, does not have genes nor tx"
$FEATURECOUNTS -a /home/gmoro/indices/combined_correct.gtf \
               -o test_flagged \
               merged_untested_header.bam \
               -F GTF \
               -t exon \
               -g gene_id \
               -f \
               -O \
               -M  \
               -T "$NTHREADS" \
               --byReadGroup



## RGs as reported by the WTA - full run, rock+roi


samtools index tso.bam
samtools index wta.bam

samtools view -h tso.bam alien_egfp_WPRE alien_tdtomato_WPRE | samtools view -Sb > small_tso.bam
samtools view -h wta.bam alien_egfp_WPRE alien_tdtomato_WPRE | samtools view -Sb > small_wta.bam



samtools merge -@ $NTHREADS - -r small_tso.bam small_wta.bam  > merged_temp.bam

valid_rgs=~/star_solo_cell_lines/RoCK_ROI/WTA_new/RoCK_ROI_WTA_new/Solo.out/Gene/filtered/barcodes.tsv
TAB="$(printf '\t')"

while read -r barcode
do
    
    echo "@RG${TAB}ID:small_tso_${barcode}"
    echo "@RG${TAB}ID:small_wta_${barcode}"

done < "$valid_rgs" > header_barcodes

sed 's/@RG\tID://g' header_barcodes > wta_tso_filteredin_barcodes


samtools view -H merged_temp.bam > curr_header.txt

# samtools view merged_temp.bam | cut -f 27,29 | grep -vP "CB:Z:-" | awk '{FS=OFS="\t"; print "@RG","ID="substr($2,6)"_"substr($1,6)}' > \
#                                                                        header_barcodes

samtools view merged_temp.bam | \
    grep -vP "CB:Z:-\t" | \
    awk '{OFS=FS="\t"; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28,"RG"substr($29,3 )"_"substr($27, 6)}' |  uniq -f 27 | \
    cat curr_header.txt header_barcodes - | \
    samtools view -Sb -@ $NTHREADS | \
    samtools sort -@ $NTHREADS  > merged_untested_header.bam

samtools index -@ $NTHREADS merged_untested_header.bam
samtools view -h merged_untested_header.bam -R wta_tso_filteredin_barcodes | \
    samtools view -Sb >  merge_untested_filtered.bam

echo "beware the gtf is exon-only, does not have genes nor tx"
grep alien /home/gmoro/indices/combined_correct.gtf > alien.gtf

$FEATURECOUNTS -a alien.gtf \
               -o alien_only \
               merged_untested_header.bam \
               -F GTF \
               -t exon \
               -g gene_id \
               -f \
               -O \
               -M  \
               -T "$NTHREADS" \
               --byReadGroup

