#!/bin/bash
set -e

cd Fastq2Data

bash script01_GettingFASTQ.sh
bash script02_Cutadaptor.sh
bash script03_FilterTrimming.sh
bash script04_Denoising.sh

Rscript script05_mergeSequenceTable.R