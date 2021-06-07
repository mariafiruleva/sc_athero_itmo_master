#!/bin/bash

cd /mnt/tank/scratch/mfiruleva/scn/data/GSE155513/SRS7124076


snakemake -j 4 --use-singularity --use-conda --conda-prefix /mnt/tank/scratch/mfiruleva/scn/config --singularity-prefix /mnt/tank/scratch/mfiruleva/scn/config --singularity-args '--bind /scratch/mfiruleva/winter/genomes/mus:/home --bind /scratch/mfiruleva/GSEs_processing/wget_and_process:/files --bind /mnt/tank/scratch/mfiruleva/scn/config/65ef9760:/mnt/tank/scratch/mfiruleva/scn/config/65ef9760 --bind /mnt/tank/scratch/mfiruleva/scn/stats/summary.csv:/mnt/tank/scratch/mfiruleva/scn/stats/summary.csv' --verbose
