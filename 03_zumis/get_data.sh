#!/bin/bash
##
## Giulia Moro's Spectral PoC
##
## Izaskun Mallona
## 24 Sept 2021

## in portmac
WD=/home/imallona/giulia/fastqs
THREADS=20

mkdir -p $WD
cd $_

cat << EOF > data_sources.conf
https://fgcz-gstore.uzh.ch/projects/p3195/NovaSeq_20210923_NOV939_o26015_DataDelivery/20210923.B-o26015_1_1-Crohn_P393_WTA_R1.fastq.gz
https://fgcz-gstore.uzh.ch/projects/p3195/NovaSeq_20210923_NOV939_o26015_DataDelivery/20210923.B-o26015_1_1-Crohn_P393_WTA_R2.fastq.gz
https://fgcz-gstore.uzh.ch/projects/p3195/NovaSeq_20210923_NOV939_o26015_DataDelivery/20210923.B-o26015_1_2-Crohn_P393_MS_R1.fastq.gz
https://fgcz-gstore.uzh.ch/projects/p3195/NovaSeq_20210923_NOV939_o26015_DataDelivery/20210923.B-o26015_1_2-Crohn_P393_MS_R2.fastq.gz
https://fgcz-gstore.uzh.ch/projects/p3195/NovaSeq_20210923_NOV939_o26015_DataDelivery/20210923.B-o26015_1_3-SPECTRAL_unmod_R1.fastq.gz
https://fgcz-gstore.uzh.ch/projects/p3195/NovaSeq_20210923_NOV939_o26015_DataDelivery/20210923.B-o26015_1_3-SPECTRAL_unmod_R2.fastq.gz
https://fgcz-gstore.uzh.ch/projects/p3195/NovaSeq_20210923_NOV939_o26015_DataDelivery/20210923.B-o26015_1_4-SPECTRAL_mod_R1.fastq.gz
https://fgcz-gstore.uzh.ch/projects/p3195/NovaSeq_20210923_NOV939_o26015_DataDelivery/20210923.B-o26015_1_4-SPECTRAL_mod_R2.fastq.gz
EOF

cat data_sources.conf | wget --user imallona -e robots=off --ask-password --no-parent -nH \
                               -i -

# downsampling for Giulia, around 1 GB each

for fn in $(find . -name "*SPECTRAL*fastq.gz")
do
    echo $fn
    curr=$(basename $fn fastq.gz)
    pigz --decompress -p $NTHREADS --keep --stdout $fn | \
        head -120000000 | pigz -p $NTHREADS --stdout > \
                              "$curr"_downsampled.fastq.gz
done


mv 20210923.B-o26015_1_4-SPECTRAL_mod_R2._downsampled.fastq.gz 20210923.B-o26015_1_4-SPECTRAL_mod_downsampled_R2.fastq.gz
mv 20210923.B-o26015_1_4-SPECTRAL_mod_R1._downsampled.fastq.gz 20210923.B-o26015_1_4-SPECTRAL_mod_downsampled_R1.fastq.gz

mv 20210923.B-o26015_1_3-SPECTRAL_unmod_R2._downsampled.fastq.gz 20210923.B-o26015_1_3-SPECTRAL_unmod_downsampled_R2.fastq.gz
mv 20210923.B-o26015_1_3-SPECTRAL_unmod_R1._downsampled.fastq.gz 20210923.B-o26015_1_3-SPECTRAL_unmod_downsampled_R1.fastq.gz

## rsync to remote server

