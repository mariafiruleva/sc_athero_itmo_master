rule run_analysis:
    """
    Run Seurat processing using count matrix from the get_count_matrix rule.
    """
    input: ['/mnt/tank/scratch/mfiruleva/scn/10_fq/GSE131777/SRS4824239/counts.RData', '/mnt/tank/scratch/mfiruleva/scn/10_fq/GSE131777/SRS4824240/counts.RData', '/mnt/tank/scratch/mfiruleva/scn/10_fq/GSE131777/SRS4824241/counts.RData', '/mnt/tank/scratch/mfiruleva/scn/10_fq/GSE131777/SRS4824242/counts.RData']

    output: h5ad=temp("GSE131777.h5ad"), rda="GSE131777.RData"
    benchmark: "benchmarks/analysis.txt"
    log: "logs/seurat.log"
    conda: "/mnt/tank/scratch/mfiruleva/scn/config/5d179225.yaml"
    singularity: "docker://mfiruleva/scn:latest"
    shell: "Rscript scripts/analysis.R 2> {log}"