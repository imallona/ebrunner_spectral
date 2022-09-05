#!/bin/bash
##
## zUMI-processes v3 BDR beads exploring the effects of external seqs (TSO or not);
##   CB binning (Hamming dist), TSO/dT bucket (i.e. linkers, tails etc)
##
## Two pass STAR, counts multimappers
##
## Downsampled to 1M reads
##
## Izaskun Mallona
## GPLv3
## 15 Aug 2022

WD=/home/imallona/ebrunner_spectral/zumis_systematic
DATA=/home/imallona/ebrunner_spectral/harmonize_fastqs
OUT="$WD"/out

set +e

mkdir -p $WD
cd $WD

mkdir -p data

cat << EOF > data/gfp_and_others.fa
>egfp
ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAA
ACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTT
CATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAG
TGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACG
TCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGG
CGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCAC
AAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGG
TGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACAC
CCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAA
GACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCA
TGGACGAGCTGTACAAGTAA
>tso_oligomer
TATGCGTAGTAGGTATGTATGCGTAGTAGGTATGTATGCGTAGTAGGTATGTATGCGTAGTAGGTATGTATGCGTAGTAGGTATGTATGCGTAGTAGGTATGTATGCGTAGTAGGTATGTATGCGTAGTAGGTATGTATGCGTAGTAGGTATGTATGCGTAGTAGGTATGTATGCGTAGTAGGTATGTATGCGTAGTAGGTATG
>n1_oligomer
ACAGGAAACTCATGGTGCGTACAGGAAACTCATGGTGCGTACAGGAAACTCATGGTGCGTACAGGAAACTCATGGTGCGTACAGGAAACTCATGGTGCGTACAGGAAACTCATGGTGCGTACAGGAAACTCATGGTGCGTACAGGAAACTCATGGTGCGTACAGGAAACTCATGGTGCGTACAGGAAACTCATGGTGCGTACAGGAAACTCATGGTGCGTACAGGAAACTCATGGTGCGT
EOF


cat << EOF > data/gfp.fa
>egfp
ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAA
ACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTT
CATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAG
TGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACG
TCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGG
CGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCAC
AAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGG
TGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACAC
CCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAA
GACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCA
TGGACGAGCTGTACAAGTAA
EOF
 
# downsampling, how much?


source /home/imallona/virtenvs/zumis/bin/activate 

export R_LIBS=/home/imallona/R/x86_64-pc-linux-gnu-library/3.6/


for run in $(find $DATA -maxdepth 1 -mindepth 1 -type d)
do
    echo "$run start"

    run_id=$(basename $run)

    # for bucket in grouped_AATG_CCAC_else \
    #                   grouped_AATG_CCAC_Tn \
    #                   grouped_AATG_CCAC_TSO \
    #                   grouped_GTGA_GACA_else \
    #                   grouped_GTGA_GACA_Tn \
    #                   grouped_GTGA_GACA_TSO \
    #                   grouped_non_matching \
    #                   grouped_other_linkers_else \
    #                   grouped_other_linkers_Tn \
    #                   grouped_other_linkers_TSO

    for bucket in grouped_AATG_CCAC_TSO \
                      grouped_GTGA_GACA_Tn
    do

        for external_seqs in gfp_and_others gfp
        do
            
            for cb_hamming in 0 4
            do
                for primary_hit in yes no
                do
                    curr_outdir="$OUT"/"$run_id"/"$bucket"/"$external_seqs"/cb_hamming_"$cb_hamming"/primary_hit_"$primary_hit"
                    mkdir -p "$curr_outdir"
                    cd "$curr_outdir"

                    curr="$run_id"_"$bucket"_"$external_seqs"_cb_hamming_"$cb_hamming"_primary_hit_"$primary_hit"
                    
                    echo "$run $bucket $external_seq $cb_hamming $primary_hit start"
              
                cat <<EOF > "$curr".yaml
project: $curr
sequence_files:
  file1:
    name: $run/output/${bucket}_R1.fastq.gz
    base_definition:
    - BC(1-9,14-22,27-35)
    - UMI(36-43)
  file2:
    name: $run/output/${bucket}_R2.fastq.gz
    base_definition: cDNA(1-62)
reference:
  STAR_index: /home/gmoro/imallona2giulia/indices/GRCm38_gencode_M25
  GTF_file: /home/gmoro/imallona2giulia/indices/gencode.vM25.annotation.gtf
  additional_STAR_params: ~
  additional_files: 
  - $WD/data/$external_seqs.fa
out_dir: $curr_outdir
num_threads: 10
mem_limit: 50
filter_cutoffs:
  BC_filter:
    num_bases: 1
    phred: 20
  UMI_filter:
    num_bases: 1
    phred: 20
barcodes:
  barcode_num: ~
  barcode_file: /home/imallona/zumis_finetuning_test/whitelist_v3_by_gmoro.txt
  automatic: yes             ## whitelist added, but automatic as well
  BarcodeBinning: $cb_hamming
  nReadsperCell: 100         ## at least 100 reads per cell
counting_opts:
  introns: yes
  downsampling: 1000000
  strand: 0
  Ham_Dist: 0                ## no Hamming for UMIs
  velocyto: no
  primaryHit: $primary_hit             ## get multimappers
  twoPass: yes               ## two pass mapping
make_stats: yes
which_Stage: Filtering
Rscript_exec: /usr/local/R/R-3.6.0/bin/Rscript
STAR_exec: /home/imallona/soft/star/bin/STAR
pigz_exec: pigz
samtools_exec: /home/imallona/soft/samtools/samtools-1.14/bin/samtools
zUMIs_directory: /home/imallona/soft/zumis_alt/zUMIs
EOF

                echo $curr
                echo 'run me'
                nice -n 19 bash /home/imallona/soft/zumis_alt/zUMIs/zUMIs.sh \
                     -y "$curr".yaml | tee > "$curr_outdir"/zumis_run.log
                
                echo 'immediately remove bam files to save space'
                cd $WD
                done
            done
        done
    done       
done


# test


run=/home/imallona/ebrunner_spectral/harmonize_fastqs/20220629.B-o2875514-SPECTRAL_4
run_id=$(basename $run)

bucket=grouped_AATG_CCAC_TSO
external_seqs=gfp_and_others
cb_hamming=0
primary_hit=no

curr="$run_id"_"$bucket"_"$external_seqs"_cb_hamming_"$cb_hamming"_primary_hit_"$primary_hit"
                  
curr_outdir="$OUT"/"$run_id"/"$bucket"/"$external_seqs"/cb_hamming_"$cb_hamming"/primary_hit_"$primary_hit"
mkdir -p "$curr_outdir"
cd "$curr_outdir"

nice -n 19 bash /home/imallona/soft/zumis_alt/zUMIs/zUMIs.sh \
     -y "$curr".yaml | tee > "$curr_outdir"/zumis_run.log
