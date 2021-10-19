#!/bin/bash
##
## Runs zUMIs on Spectral data, including the index generation steps
##
## 19th Oct 2021
## GPLv3
## Izaskun Mallona

echo 'Run on sherborne, user imallona'

echo 'Index human'

bash index_human.sh

echo 'Index mouse'

bash index_mouse.sh

echo 'Index alien'

bash index_alien.sh

echo 'Run one of the datasets'

mkdir -p /home/imallona/giulia/zumi_runs/20210923.B-o26015_1_3-SPECTRAL_unmod

export R_LIBS=~/R/x86_64-pc-linux-gnu-library/4.1/
bash /home/imallona/soft/zUMIs/zUMIs.sh \
     -y /home/imallona/src/ebrunner_spectral/03_zumis/prototype_sherborne_2.yaml


/usr/local/R/R-4.1.0/bin/Rscript /home/imallona/soft/zUMIs/checkyaml.R \
                                 /home/imallona/src/ebrunner_spectral/03_zumis/prototype_sherborne_2.yaml
