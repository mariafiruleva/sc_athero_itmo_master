rule run_analysis:
    """
    Run Seurat processing using count matrix from the get_count_matrix rule.
    """
    input: rules.get_count_matrix.output.count_mat

    output: h5ad=temp("SRS4824197.h5ad"), rda="SRS4824197.RData"
    benchmark: "benchmarks/analysis.txt"
    log: "logs/seurat.log"
    conda: "/mnt/tank/scratch/mfiruleva/scn/config/65ef9760.yaml"
    singularity: "docker://mfiruleva/scn:latest"
    shell: "Rscript scripts/analysis.R 2> {log}"