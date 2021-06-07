rule run_analysis:
    """
    Run Seurat processing using count matrix from the get_count_matrix rule.
    """
    input: ['/mnt/tank/scratch/mfiruleva/scn/data/GSE149069/SRS6514323/counts.RData', '/mnt/tank/scratch/mfiruleva/scn/data/GSE149069/SRS6514322/counts.RData']

    output: h5ad=temp("GSE149069.h5ad"), rda="GSE149069.RData"
    benchmark: "benchmarks/analysis.txt"
    log: "logs/seurat.log"
    conda: "/mnt/tank/scratch/mfiruleva/scn/config/65ef9760.yaml"
    singularity: "docker://mfiruleva/scn:latest"
    shell: "Rscript scripts/analysis.R 2> {log}"