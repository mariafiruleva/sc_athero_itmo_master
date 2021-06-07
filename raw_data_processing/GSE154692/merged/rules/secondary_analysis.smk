rule run_analysis:
    """
    Run Seurat processing using count matrix from the get_count_matrix rule.
    """
    input: ['/mnt/tank/scratch/mfiruleva/scn/data/GSE154692/SRS7041373/counts.RData', '/mnt/tank/scratch/mfiruleva/scn/data/GSE154692/SRS7041374/counts.RData', '/mnt/tank/scratch/mfiruleva/scn/data/GSE154692/SRS7041375/counts.RData', '/mnt/tank/scratch/mfiruleva/scn/data/GSE154692/SRS7041376/counts.RData']

    output: h5ad=temp("GSE154692.h5ad"), rda="GSE154692.RData"
    benchmark: "benchmarks/analysis.txt"
    log: "logs/seurat.log"
    conda: "/mnt/tank/scratch/mfiruleva/scn/config/5d179225.yaml"
    singularity: "docker://mfiruleva/scn:latest"
    shell: "Rscript scripts/analysis.R 2> {log}"