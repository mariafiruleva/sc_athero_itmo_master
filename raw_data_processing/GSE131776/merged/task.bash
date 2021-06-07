#!/bin/bash

cd /mnt/tank/scratch/mfiruleva/scn/data/GSE131776/merged


snakemake -j 4 --use-singularity --use-conda --conda-prefix /mnt/tank/scratch/mfiruleva/scn/config --singularity-prefix /mnt/tank/scratch/mfiruleva/scn/config --singularity-args '--bind /mnt/tank/scratch/mfiruleva/scn/config/5d179225:/mnt/tank/scratch/mfiruleva/scn/config/5d179225 --bind /mnt/tank/scratch/mfiruleva/scn/data/GSE131776:/mnt/tank/scratch/mfiruleva/scn/data/GSE131776 --bind /mnt/tank/scratch/mfiruleva/scn/stats/summary.csv:/mnt/tank/scratch/mfiruleva/scn/stats/summary.csv' --verbose
