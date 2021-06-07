rule run_analysis:
    """
    Run Seurat processing using count matrix from the get_count_matrix rule.
    """
    input: ['/mnt/tank/scratch/mfiruleva/scn/data/GSE131776/SRS4824181/counts.RData', '/mnt/tank/scratch/mfiruleva/scn/data/GSE131776/SRS4824182/counts.RData', '/mnt/tank/scratch/mfiruleva/scn/data/GSE131776/SRS4824183/counts.RData', '/mnt/tank/scratch/mfiruleva/scn/data/GSE131776/SRS4824184/counts.RData', '/mnt/tank/scratch/mfiruleva/scn/data/GSE131776/SRS4824185/counts.RData', '/mnt/tank/scratch/mfiruleva/scn/data/GSE131776/SRS4824186/counts.RData', '/mnt/tank/scratch/mfiruleva/scn/data/GSE131776/SRS4824187/counts.RData', '/mnt/tank/scratch/mfiruleva/scn/data/GSE131776/SRS4824188/counts.RData', '/mnt/tank/scratch/mfiruleva/scn/data/GSE131776/SRS4824189/counts.RData', '/mnt/tank/scratch/mfiruleva/scn/data/GSE131776/SRS4824190/counts.RData', '/mnt/tank/scratch/mfiruleva/scn/data/GSE131776/SRS4824191/counts.RData', '/mnt/tank/scratch/mfiruleva/scn/data/GSE131776/SRS4824192/counts.RData', '/mnt/tank/scratch/mfiruleva/scn/data/GSE131776/SRS4824193/counts.RData', '/mnt/tank/scratch/mfiruleva/scn/data/GSE131776/SRS4824194/counts.RData', '/mnt/tank/scratch/mfiruleva/scn/data/GSE131776/SRS4824195/counts.RData', '/mnt/tank/scratch/mfiruleva/scn/data/GSE131776/SRS4824196/counts.RData', '/mnt/tank/scratch/mfiruleva/scn/data/GSE131776/SRS4824197/counts.RData', '/mnt/tank/scratch/mfiruleva/scn/data/GSE131776/SRS4824198/counts.RData']

    output: h5ad=temp("GSE131776.h5ad"), rda="GSE131776.RData"
    benchmark: "benchmarks/analysis.txt"
    log: "logs/seurat.log"
    conda: "/mnt/tank/scratch/mfiruleva/scn/config/5d179225.yaml"
    singularity: "docker://mfiruleva/scn:latest"
    shell: "Rscript scripts/analysis.R 2> {log}"