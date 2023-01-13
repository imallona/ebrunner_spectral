#!/bin/bash
##
## Preliminary index build + STARsolo run
##
## Izaskun Mallona
## 05 Jan 2023

# https://github.com/alexdobin/STAR/issues/1607#issuecomment-1310028944

SRC=/home/imallona/src/ebrunner_spectral
r1=/home/gmoro/fastqs_PDGFRA/single_clean/o303001_2-RoCKseq_S2_R1_001.fastq.gz
r2=/home/gmoro/fastqs_PDGFRA/single_clean/o303001_2-RoCKseq_S2_R2_001.fastq.gz

export PATH=/home/imallona/soft/star/STAR-2.7.10b/source:$PATH

mkdir -p ~/ebrunner_spectral/star_solo
cd $_

## index start


FA=~/giulia/indices/GRCm38.p6.genome.fa
GTF=~/giulia/indices/gencode.vM25.annotation.gtf
NTHREADS=20
ID=GRCm38_gencode_M25_gfp_star27b

cd ~/giulia/indices/


cat << EOF > alien.fa
>gfp
ATGCCAGAGCCAGCGAAGTCTGCTCCCGCCCCGAAAAAGGGCTCCAAGAAGGCGGTGACTAAGGCGCAGA
AGAAAGGCGGCAAGAAGCGCAAGCGCAGCCGCAAGGAGAGCTATTCCATCTATGTGTACAAGGTTCTGAA
GCAGGTCCACCCTGACACCGGCATTTCGTCCAAGGCCATGGGCATCATGAATTCGTTTGTGAACGACATT
TTCGAGCGCATCGCAGGTGAGGCTTCCCGCCTGGCGCATTACAACAAGCGCTCGACCATCACCTCCAGGG
AGATCCAGACGGCCGTGCGCCTGCTGCTGCCTGGGGAGTTGGCCAAGCACGCCGTGTCCGAGGGTACTAA
GGCCATCACCAAGTACACCAGCGCTAAGGATCCACCGGTCGCCACCATGGTGAGCAAGGGCGAGGAGCTG
TTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCG
GCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCC
CGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCAC
ATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCA
AGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGA
GCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGC
CACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACA
TCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCT
GCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCAC
ATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAAGCGG
CCAATTCCAGCTGAGCGCCGGTCGCTACCATTACCAGTTGGTCTGGTGTCAGGGGATCCCCCCTCGCTGA
TCAGCCTCGACTGTGCCTTCTAGTTGCCAGCCATCTGTTGTTTGCCCCTCCCCCGTGCCTTCCTTGACCC
TGGAAGGTGCCACTCCCACTGTCCTTTCCTAATAAAATGAAGAAATTGCATCGCATTGTCTGAGTAGGTG
TCATTCTATTCTGGGGGGTGGGGTGGGGCAGGACAGCAAGGGGGAGGATTGGGAAGACAATAGCAGGCAT
GCTGGGGATGCGGTGGGCTCTATGG

EOF


echo -e 'gfp\tALIEN\texon\t1\t1425\t.\t+\t.\tgene_id "gfp"; transcript_id "gfpt";' > alien.gtf

cat $FA alien.fa > combined.fa
cat $GTF alien.gtf > combined.gtf



nice -n19 STAR --runMode genomeGenerate \
      --runThreadN 50 \
      --genomeDir "$ID" \
      --genomeFastaFiles combined.fa

# rm combined.fa

## indexing end



cd ~/ebrunner_spectral/star_solo

nice -n 19 STAR --runThreadN 5 \
     --genomeDir /home/imallona/giulia/indices/GRCm38_gencode_M25_gfp_star27b/ \
     --readFilesCommand zcat \
     --outFileNamePrefix star_outs_wta_test/ \
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
     --sjdbGTFfile ~/giulia/indices/combined.gtf \
     --outTmpDir ~/tmp/tmp

## this is to estimate 10k cells, not 3k (hardcoded by default)
STAR --runMode soloCellFiltering  star_outs_wta_test/Solo.out/Gene/raw \
     star_outs_wta_test_10k/  \
     --soloCellFilter EmptyDrops_CR  10000 0.99 10 45000 90000 500 0.01 20000 0.05 10000

nice -n19 STAR --runThreadN 5 \
     --genomeDir /home/imallona/giulia/indices/GRCm38_gencode_M25_gfp_star27b/ \
     --readFilesCommand zcat \
     --outFileNamePrefix star_outs_tso_test/ \
     --readFilesIn "$r2" "$r1" \
     --soloType CB_UMI_Complex \
     --soloAdapterSequence AATGNNNNNNNNNCCAC \
     --soloCBposition 2_-9_2_-1 2_4_2_12 2_17_2_25 \
     --soloUMIposition 3_10_3_17 \
     --soloCBwhitelist "$SRC"/07_barcodes_translation_sbg/data/CLS1.txt "$SRC"/07_barcodes_translation_sbg/data/CLS2.txt "$SRC"/07_barcodes_translation_sbg/data/CLS3.txt \
     --soloCBmatchWLtype 1MM \
     --soloCellFilter None \
     --outSAMattributes CB UB gx gn sS CR CY UR UY\
     --outSAMtype BAM SortedByCoordinate \
     --quantMode GeneCounts \
     --sjdbGTFfile ~/giulia/indices/combined.gtf \
     --outTmpDir ~/tmp/another


# Extract the CBs-assigned alns for TSO and WTA bam files

cd ~/ebrunner_spectral/star_solo/star_outs_wta_test

samtools view -h Aligned.sortedByCoord.out.bam | grep -vP "CR:Z:\t" | samtools view  -Sb > with_CR.bam
samtools index -@ 10 with_CR.bam


cd ~/ebrunner_spectral/star_solo/star_outs_tso_test

samtools view -h Aligned.sortedByCoord.out.bam | grep -vP "CR:Z:\t" | samtools view  -Sb > with_CR.bam
samtools index -@ 10 with_CR.bam

# awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k4,4n -k5,5n"}' "$ID".gtf > out_sorted.gtf


## Strand-specific mapping

nice -n19 STAR --runThreadN 5 \
     --genomeDir /home/imallona/giulia/indices/GRCm38_gencode_M25_gfp_star27b/ \
     --readFilesCommand zcat \
     --outFileNamePrefix star_outs_tso_strand_forward_more_flags \
     --readFilesIn "$r2" "$r1" \
     --soloType CB_UMI_Complex \
     --soloAdapterSequence AATGNNNNNNNNNCCAC \
     --soloCBposition 2_-9_2_-1 2_4_2_12 2_17_2_25 \
     --soloUMIposition 3_10_3_17 \
     --soloCBwhitelist "$SRC"/07_barcodes_translation_sbg/data/CLS1.txt "$SRC"/07_barcodes_translation_sbg/data/CLS2.txt "$SRC"/07_barcodes_translation_sbg/data/CLS3.txt \
     --soloCellFilter None \
     --outSAMattributes CB UB gx gn sM sS sQ CR CY UR UY\
     --outSAMtype BAM SortedByCoordinate \
     --quantMode GeneCounts \
     --sjdbGTFfile ~/giulia/indices/combined.gtf \
     --outTmpDir ~/tmp/another \
     --soloStrand Forward \
     --soloUMIdedup 1MM_All \
     --soloAdapterMismatchesNmax 1 \
     --soloCBmatchWLtype 1MM
     
