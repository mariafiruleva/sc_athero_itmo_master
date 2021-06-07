rule get_cells:
    input: rda='data/{gse_id}.RData', annot='annotations/{gse_id}.tsv'
    output: rda='target/{gse_id}_target.RData'
    benchmark: "benchmarks/get_cells/{gse_id}.txt"
    log: 'logs/get_cells/{gse_id}.txt'
    params: res=lambda wildcards: data[wildcards.gse_id]['res'],
            clusters=lambda wildcards: data[wildcards.gse_id]['clusters']
    message: 'Start get_cells rule, dataset: {wildcards.gse_id}, ident = {params.res}, clusters = {params.clusters}'
    shell: "/scratch/opt/R/3.6.0/bin/Rscript scripts/get_cells.R --in_rda {input.rda} --out_rda {output.rda} --ident {params.res} --clusters {params.clusters} --annot {input.annot} 2> {log}"
    
rule integrate:
    input: rda=expand('target/{gse_id}_target.RData', gse_id=data.keys())
    output: rda='out/athero_merged.RData', h5='out/athero_merged.h5ad'
    benchmark: "benchmarks/integrate.txt"
    log: 'logs/integrate.txt'
    message: 'Start integration with target samples'
    conda: "/mnt/tank/scratch/mfiruleva/scn/config/5d179225.yaml"
    singularity: "docker://mfiruleva/scn:latest"
    shell: "Rscript scripts/merge_samples.R 2> {log}"