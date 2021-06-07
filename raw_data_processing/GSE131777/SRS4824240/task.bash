#!/bin/bash

cd /mnt/tank/scratch/mfiruleva/scn/10_fq/GSE131777/SRS4824240

snakemake --snakefile Snakefile -F -j 4
