#!/bin/bash

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

rm combined.fa


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
     --outTmpDir ~/tmp/thisone


nice -n19 STAR --runThreadN 5 \
     --genomeDir /home/imallona/giulia/indices/GRCm38_gencode_M25_gfp_star27b/ \
     --readFilesCommand zcat \
     --outFileNamePrefix star_outs_tso_test/ \
     --readFilesIn "$r2" "$r1" \
     --soloType CB_UMI_Complex \
     --soloAdapterSequence AATGANNNNNNNNNCCAC \
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
