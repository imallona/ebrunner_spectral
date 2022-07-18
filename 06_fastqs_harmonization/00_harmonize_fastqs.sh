#!/bin/bash
##
## New rhapsody beads parsing/harmonizing (removal of prepended seqs)
##
## Izaskun Mallona
##
## 04 July 2022


# reads are a mixture of:
# empty - seq1 9nt - linker1 4nt -  seq2 9nt - linker2 4nt - seq3 9nt - UMI 8 nt - T 25nt
# empty - seq1 9nt - linker1 4nt -  seq2 9nt - linker2 4nt - seq3 9nt - UMI 8 nt - TATGCGTAGTAGGTATG
# A     - seq1 9nt - linker1 4nt -  seq2 9nt - linker2 4nt - seq3 9nt - UMI 8 nt - T 25nt
# A     - seq1 9nt - linker1 4nt -  seq2 9nt - linker2 4nt - seq3 9nt - UMI 8 nt - TATGCGTAGTAGGTATG
# GT    - seq1 9nt - linker1 4nt -  seq2 9nt - linker2 4nt - seq3 9nt - UMI 8 nt - T 25nt
# GT    - seq1 9nt - linker1 4nt -  seq2 9nt - linker2 4nt - seq3 9nt - UMI 8 nt - TATGCGTAGTAGGTATG
# TCA   - seq1 9nt - linker1 4nt -  seq2 9nt - linker2 4nt - seq3 9nt - UMI 8 nt - T 25nt
# TCGA  - seq1 9nt - linker1 4nt -  seq2 9nt - linker2 4nt - seq3 9nt - UMI 8 nt - TATGCGTAGTAGGTATG

# A regex might look like
# '([ACTGN]{0,3})([ACTGN]{9})(GTGA|AATG{1})([ACTGN]{9})(GACA|CCAC){1}([ACTGN]{9})([ACTGN]{8})(.*)'

WD=/home/imallona/tmp/regex_more_buckets
DATA=/home/gmoro/fastqs_six_scRNAseq

mkdir -p $WD; cd $WD
# zcat "$DATA"/20220629.B-o2875511-SPECTRAL_1_R1.fastq.gz | head -100 > test.fq


# what about chunking each R1/R2 file into a multiple-of-4 number of lines, and then running the parsing script?
NLINES=10000000
zcat "$DATA"/20220629.B-o2875511-SPECTRAL_1_R1.fastq.gz | split - -l "$NLINES" --filter='gzip > $FILE.r1.gz' part.
zcat "$DATA"/20220629.B-o2875511-SPECTRAL_1_R2.fastq.gz  | split - -l "$NLINES" --filter='gzip > $FILE.r2.gz' part.

# for r1 in $(find . -name "part.*.r1.gz")
# do
#     echo $r1
#     curr=$(basename $r1 .r1.gz)
#     r2="$curr".r2.gz
#     echo $r2

#     Rscript ~/src/ebrunner_spectral/06_fastqs_harmonization/01_harmonize_fastqs.R \
#             -r1 "$r1" \
#             -r2 "$r2" \
#             -o delete_me_chunked
# done

## ok, now parallelize this - 30 cores

N=32

(
    for r1 in $(find . -name "part.*.r1.gz" | xargs -n"$N") 
    do 
        ((i=i%N)); ((i++==0)) && wait
        echo $r1
        curr=$(basename $r1 .r1.gz)
        r2="$curr".r2.gz
        echo $r2
    
        echo $i
        
        nice -n 19 /usr/local/R/R-4.1.0/bin/Rscript ~/src/ebrunner_spectral/06_fastqs_harmonization/01_harmonize_fastqs.R \
            -r1 "$r1" \
            -r2 "$r2" \
            -o "$WD"/delete_me_chunked/"$curr" &

    done
)


## ran till here

# mkdir -p output

# for item in non_matching AATG_CCAC_linkers GTGA_GACA_linkers
# do
#     find delete_me_chunked -name "$item*R2.fastq.gz" | sort | xargs zcat | \
#         gzip -c > output/grouped_"$item"_R2.fastq.gz
#     find delete_me_chunked -name "$item*R1.fastq.gz" | sort | xargs zcat | \
#         gzip -c > output/grouped_"$item"_R1.fastq.gz
# done

# rm *part*gz
# rm -rf delete_me_chunked
