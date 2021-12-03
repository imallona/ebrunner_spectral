#!/bin/bash
##
## Runs zUMIs on Spectral data, including the index generation steps
##
## 19th Oct 2021
## GPLv3
## Izaskun Mallona

NTHREADS=5

echo 'Run on sherborne, user imallona'

echo 'Index human'

bash index_human.sh

echo 'Index mouse'

bash index_mouse.sh

echo 'Index alien'

bash index_alien.sh

echo 'Run one of the datasets / downsampled'

# mkdir -p /home/imallona/giulia/zumi_runs/20210923.B-o26015_1_3-SPECTRAL_unmod

# export R_LIBS=/home/imallona/R/x86_64-pc-linux-gnu-library/4.1/
# bash /home/imallona/soft/zUMIs/zUMIs.sh \
#      -y /home/imallona/src/ebrunner_spectral/03_zumis/prototype_sherborne_2.yaml


# /usr/local/R/R-4.1.0/bin/Rscript /home/imallona/soft/zUMIs/checkyaml.R \
#                                  /home/imallona/src/ebrunner_spectral/03_zumis/prototype_sherborne_2.yaml

cd /home/imallona/giulia/fastqs

mkdir -p /home/imallona/giulia/fastqs/downsampled

for fn in $(find . -maxdepth 1 -name "*SPECTRAL_unmod*fastq.gz" )
do
    echo $fn
    curr=$(basename $fn .fastq.gz)
    pigz --decompress -p $NTHREADS --keep --stdout $fn | \
        head -40000000 | pigz -p $NTHREADS --stdout > \
                              /home/imallona/giulia/fastqs/downsampled/"$curr"_downsampled.fastq.gz
done


## try the example data run

# mkdir -p /home/imallona/giulia/zumis_example
# cd $_

# wget https://github.com/sdparekh/zUMIs/raw/zUMIs-version1/ExampleData/barcoderead_HEK.1mio.fq.gz
# wget https://github.com/sdparekh/zUMIs/raw/zUMIs-version1/ExampleData/cDNAread_HEK.1mio.fq.gz
# source /home/imallona/virtenvs/zumis/bin/activate
# export R_LIBS=/home/imallona/R/x86_64-pc-linux-gnu-library/3.6/

# bash /home/imallona/soft/zUMIs/zUMIs.sh \
#      -y /home/imallona/src/ebrunner_spectral/03_zumis/yaml/runExample_local.yaml

# # it worked


## run prototype / downsampled on mouse

## store temporary files here, and not in /tmp
# mkdir /home/imallona/tmp
# export TMPDIR=/home/imallona/tmp
# export TMP=/home/imallona/tmp

## point to the R libraries
mkdir -p /home/imallona/giulia/zumi_runs/downsampled/alien/

## zUMIs expects uncompressed GTFs
gunzip /home/imallona/giulia/indices/alien.gtf.gz

source /home/imallona/virtenvs/zumis/bin/activate
export R_LIBS=/home/imallona/R/x86_64-pc-linux-gnu-library/3.6/

mkdir -p /home/imallona/giulia/logs
rm -rf /home/imallona/giulia/zumi_runs/downsampled/mouse
mkdir -p /home/imallona/giulia/zumi_runs/downsampled/mouse/

bash /home/imallona/soft/zUMIs/zUMIs.sh \
     -y /home/imallona/src/ebrunner_spectral/03_zumis/yaml/prototype_downsampled_mouse.yaml \
    &> /home/imallona/giulia/logs/zumis_mapping_mouse_r36.log


## looping through the real datasets

# fastqs_path="/home/imallona/giulia/fastqs/"
# indices_path="/home/imallona/giulia/indices/"
# indices="GRCh38_gencode_38 GRCm38_gencode_M25 alien"
# treatments="20210923.B-o26015_1_3-SPECTRAL_unmod 20210923.B-o26015_1_4-SPECTRAL_mod"

## relative to the `03_zumis` code repository path
yaml_paths=./yaml/batch_1/no_barcodes_hm_0

mkdir -p /home/imallona/giulia/zumi_runs/batch_1/no_barcodes_hm_0
mkdir -p /home/imallona/giulia/zumi_runs/batch_1/logs

source /home/imallona/virtenvs/zumis/bin/activate
export R_LIBS=/home/imallona/R/x86_64-pc-linux-gnu-library/3.6/

for yaml in $(find $yaml_paths -name "*yaml")
do
    echo "$yaml"

    bash /home/imallona/soft/zUMIs/zUMIs.sh \
         -y "$yaml" \
        &> /home/imallona/giulia/zumi_runs/batch_1/logs/"$(basename $yaml .yaml)".log

done

## extra yamls for combined runs

for yaml in $(find $yaml_paths -name "combined*mod.yaml")
do
    echo "$yaml"
    bash /home/imallona/soft/zUMIs/zUMIs.sh \
         -y "$yaml" \
        &> /home/imallona/giulia/zumi_runs/batch_1/logs/"$(basename $yaml .yaml)".log
    
done
