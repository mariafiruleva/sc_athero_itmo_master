rule run_analysis:
    """
    Run Seurat processing using count matrix from the get_count_matrix rule.
    """
    input: ['/mnt/tank/scratch/mfiruleva/scn/data/GSE154817/SRS7052767/counts.RData', '/mnt/tank/scratch/mfiruleva/scn/data/GSE154817/SRS7052768/counts.RData']

    output: h5ad=temp("GSE154817.h5ad"), rda="GSE154817.RData"
    benchmark: "benchmarks/analysis.txt"
    log: "logs/seurat.log"
    conda: "/mnt/tank/scratch/mfiruleva/scn/config/5d179225.yaml"
    singularity: "docker://mfiruleva/scn:latest"
    shell: "Rscript scripts/analysis.R 2> {log}"