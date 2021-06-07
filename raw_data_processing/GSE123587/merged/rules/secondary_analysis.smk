rule run_analysis:
    """
    Run Seurat processing using count matrix from the get_count_matrix rule.
    """
    input: ['/mnt/tank/scratch/mfiruleva/scn/data/GSE123587/SRS4137616/counts.RData', '/mnt/tank/scratch/mfiruleva/scn/data/GSE123587/SRS4137620/counts.RData']

    output: h5ad=temp("GSE123587.h5ad"), rda="GSE123587.RData"
    benchmark: "benchmarks/analysis.txt"
    log: "logs/seurat.log"
    conda: "/mnt/tank/scratch/mfiruleva/scn/config/5d179225.yaml"
    singularity: "docker://mfiruleva/scn:latest"
    shell: "Rscript scripts/analysis.R 2> {log}"