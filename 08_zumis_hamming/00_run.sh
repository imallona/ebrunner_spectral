#!/bin/bash

mkdir -p /home/imallona/zumis_finetuning_test
cd $_

cat <<EOF > 00_run.yaml
project: 00_run
sequence_files:
  file1:
    name: /home/imallona/ebrunner_spectral/harmonize_fastqs/20220629.B-o2875513-SPECTRAL_3/output/grouped_AATG_CCAC_TSO_R1.fastq.gz
    base_definition:
    - BC(1-9,14-22,27-35)
    - UMI(36-43)
  file2:
    name: /home/imallona/ebrunner_spectral/harmonize_fastqs/20220629.B-o2875513-SPECTRAL_3/output/grouped_AATG_CCAC_TSO_R2.fastq.gz
    base_definition: cDNA(1-62)
reference:
  STAR_index: /home/gmoro/imallona2giulia/indices/GRCm38_gencode_M25
  GTF_file: /home/gmoro/imallona2giulia/indices/gencode.vM25.annotation.gtf
  additional_STAR_params: ~
  additional_files: 
  - /home/gmoro/indices/egfp_no_cap.fa
out_dir: /home/imallona/zumis_finetuning_test/zumis
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
  BarcodeBinning: 4          ## with Hamming = 4
  nReadsperCell: 100         ## at least 100 reads per cell
counting_opts:
  introns: yes
  downsampling: '0'
  strand: 0
  Ham_Dist: 0                ## no Hamming for UMIs
  velocyto: no
  primaryHit: no             ## get multimappers
  twoPass: yes               ## two pass mapping
make_stats: yes
which_Stage: Filtering
Rscript_exec: /usr/local/R/R-3.6.0/bin/Rscript
STAR_exec: /home/imallona/soft/star/bin/STAR
pigz_exec: pigz
samtools_exec: /home/imallona/soft/samtools/samtools-1.14/bin/samtools
zUMIs_directory: /home/imallona/soft/zumis_alt/zUMIs
EOF

source /home/imallona/virtenvs/zumis/bin/activate 

export R_LIBS=/home/imallona/R/x86_64-pc-linux-gnu-library/3.6/ 

# rm -rf /home/gmoro/mapping_sixth_scRNAseq/10_100/TSO/zUMIs  ## to avoid the error you're having
mkdir -p /home/imallona/zumis_finetuning_test/zumis

cd $_

nice -n 19 bash /home/imallona/soft/zumis_alt/zUMIs/zUMIs.sh \
     -y ../00_run.yaml | tee > /home/imallona/zumis_finetuning_test/00_run.log


