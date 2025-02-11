#!/bin/bash
CONDA_PATH=$(conda info | grep -i 'base environment' | awk '{print $4}')
source $CONDA_PATH/etc/profile.d/conda.sh
#Activate conda environment
conda activate irmaenv2
IRMA FLU PATH/TO/reads/avian/AV-m-EX1_R1_001.fastq.gz	PATH/TO/reads/avian/AV-m-EX1_R2_001.fastq.gz	AV-m-EX1 &
wait
IRMA FLU PATH/TO/reads/avian/AV-m-EX2_R1_001.fastq.gz	PATH/TO/reads/avian/AV-m-EX2_R2_001.fastq.gz	AV-m-EX2 &
wait
IRMA FLU PATH/TO/reads/avian/AV-s-EX1_R1_001.fastq.gz	PATH/TO/reads/avian/AV-s-EX1_R2_001.fastq.gz	AV-s-EX1 &
wait
IRMA PATH/TO/reads/avian/AV-s-EX2_R1_001.fastq.gz	PATH/TO/reads/avian/AV-s-EX2_R2_001.fastq.gz	AV-s-EX2 &
wait
IRMA FLU PATH/TO/reads/avian/AV-w-EX1_R1_001.fastq.gz	PATH/TO/reads/avian/AV-w-EX1_R2_001.fastq.gz	AV-w-EX1 &
wait
IRMA FLU PATH/TO/reads/avian/AV-w-EX2_R1_001.fastq.gz	PATH/TO/reads/avian/AV-w-EX2_R2_001.fastq.gz	AV-w-EX2 &
wait
. /PATH/TO/IRMAscript.sh
