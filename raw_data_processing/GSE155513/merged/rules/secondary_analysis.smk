rule run_analysis:
    """
    Run Seurat processing using count matrix from the get_count_matrix rule.
    """
    input: ['/mnt/tank/scratch/mfiruleva/scn/data/GSE155513/SRS7124065/counts.RData', '/mnt/tank/scratch/mfiruleva/scn/data/GSE155513/SRS7124066/counts.RData', '/mnt/tank/scratch/mfiruleva/scn/data/GSE155513/SRS7124067/counts.RData', '/mnt/tank/scratch/mfiruleva/scn/data/GSE155513/SRS7124068/counts.RData', '/mnt/tank/scratch/mfiruleva/scn/data/GSE155513/SRS7124069/counts.RData', '/mnt/tank/scratch/mfiruleva/scn/data/GSE155513/SRS7124070/counts.RData', '/mnt/tank/scratch/mfiruleva/scn/data/GSE155513/SRS7124071/counts.RData', '/mnt/tank/scratch/mfiruleva/scn/data/GSE155513/SRS7124072/counts.RData', '/mnt/tank/scratch/mfiruleva/scn/data/GSE155513/SRS7124073/counts.RData', '/mnt/tank/scratch/mfiruleva/scn/data/GSE155513/SRS7124074/counts.RData', '/mnt/tank/scratch/mfiruleva/scn/data/GSE155513/SRS7124075/counts.RData', '/mnt/tank/scratch/mfiruleva/scn/data/GSE155513/SRS7124076/counts.RData', '/mnt/tank/scratch/mfiruleva/scn/data/GSE155513/SRS7124077/counts.RData', '/mnt/tank/scratch/mfiruleva/scn/data/GSE155513/SRS7124078/counts.RData']

    output: h5ad=temp("GSE155513.h5ad"), rda="GSE155513.RData"
    benchmark: "benchmarks/analysis.txt"
    log: "logs/seurat.log"
    conda: "/mnt/tank/scratch/mfiruleva/scn/config/5d179225.yaml"
    singularity: "docker://mfiruleva/scn:latest"
    shell: "Rscript scripts/analysis.R 2> {log}"