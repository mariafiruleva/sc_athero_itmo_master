rule run_analysis:
    """
    Run Seurat processing using count matrix from the get_count_matrix rule.
    """
    input: ['/mnt/tank/scratch/mfiruleva/scn/data/GSE116240/SRS3458842/counts.RData',
            '/mnt/tank/scratch/mfiruleva/scn/data/GSE116240/SRS3458845/counts.RData',
            '/mnt/tank/scratch/mfiruleva/scn/data/GSE149069/SRS6514323/counts.RData',
            '/mnt/tank/scratch/mfiruleva/scn/data/GSE149069/SRS6514322/counts.RData',
            '/mnt/tank/scratch/mfiruleva/scn/data/GSE131776/SRS4824181/counts.RData',
            '/mnt/tank/scratch/mfiruleva/scn/data/GSE131776/SRS4824182/counts.RData',
            '/mnt/tank/scratch/mfiruleva/scn/data/GSE131776/SRS4824183/counts.RData',
            '/mnt/tank/scratch/mfiruleva/scn/data/GSE131776/SRS4824184/counts.RData',
            '/mnt/tank/scratch/mfiruleva/scn/data/GSE131776/SRS4824185/counts.RData',
            '/mnt/tank/scratch/mfiruleva/scn/data/GSE131776/SRS4824186/counts.RData',
            '/mnt/tank/scratch/mfiruleva/scn/data/GSE131776/SRS4824187/counts.RData',
            '/mnt/tank/scratch/mfiruleva/scn/data/GSE131776/SRS4824188/counts.RData',
            '/mnt/tank/scratch/mfiruleva/scn/data/GSE131776/SRS4824189/counts.RData',
            '/mnt/tank/scratch/mfiruleva/scn/data/GSE131776/SRS4824190/counts.RData',
            '/mnt/tank/scratch/mfiruleva/scn/data/GSE131776/SRS4824191/counts.RData',
            '/mnt/tank/scratch/mfiruleva/scn/data/GSE131776/SRS4824192/counts.RData',
            '/mnt/tank/scratch/mfiruleva/scn/data/GSE131776/SRS4824193/counts.RData',
            '/mnt/tank/scratch/mfiruleva/scn/data/GSE131776/SRS4824194/counts.RData',
            '/mnt/tank/scratch/mfiruleva/scn/data/GSE131776/SRS4824195/counts.RData',
            '/mnt/tank/scratch/mfiruleva/scn/data/GSE131776/SRS4824196/counts.RData',
            '/mnt/tank/scratch/mfiruleva/scn/data/GSE131776/SRS4824197/counts.RData',
            '/mnt/tank/scratch/mfiruleva/scn/data/GSE131776/SRS4824198/counts.RData',
            '/mnt/tank/scratch/mfiruleva/scn/10_fq/GSE131777/SRS4824239/counts.RData',
            '/mnt/tank/scratch/mfiruleva/scn/10_fq/GSE131777/SRS4824240/counts.RData',
            '/mnt/tank/scratch/mfiruleva/scn/10_fq/GSE131777/SRS4824241/counts.RData',
            '/mnt/tank/scratch/mfiruleva/scn/10_fq/GSE131777/SRS4824242/counts.RData',
            '/mnt/tank/scratch/mfiruleva/scn/data/GSE123587/SRS4137616/counts.RData',
            '/mnt/tank/scratch/mfiruleva/scn/data/GSE123587/SRS4137620/counts.RData',
            '/mnt/tank/scratch/mfiruleva/scn/data/GSE155513/SRS7124065/counts.RData',
            '/mnt/tank/scratch/mfiruleva/scn/data/GSE155513/SRS7124066/counts.RData',
            '/mnt/tank/scratch/mfiruleva/scn/data/GSE155513/SRS7124067/counts.RData',
            '/mnt/tank/scratch/mfiruleva/scn/data/GSE155513/SRS7124068/counts.RData',
            '/mnt/tank/scratch/mfiruleva/scn/data/GSE155513/SRS7124069/counts.RData',
            '/mnt/tank/scratch/mfiruleva/scn/data/GSE155513/SRS7124070/counts.RData',
            '/mnt/tank/scratch/mfiruleva/scn/data/GSE155513/SRS7124071/counts.RData',
            '/mnt/tank/scratch/mfiruleva/scn/data/GSE155513/SRS7124072/counts.RData',
            '/mnt/tank/scratch/mfiruleva/scn/data/GSE155513/SRS7124073/counts.RData',
            '/mnt/tank/scratch/mfiruleva/scn/data/GSE155513/SRS7124074/counts.RData',
            '/mnt/tank/scratch/mfiruleva/scn/data/GSE155513/SRS7124075/counts.RData',
            '/mnt/tank/scratch/mfiruleva/scn/data/GSE155513/SRS7124076/counts.RData',
            '/mnt/tank/scratch/mfiruleva/scn/data/GSE155513/SRS7124077/counts.RData',
            '/mnt/tank/scratch/mfiruleva/scn/data/GSE155513/SRS7124078/counts.RData',
            '/mnt/tank/scratch/mfiruleva/scn/data/GSE154692/SRS7041373/counts.RData',
            '/mnt/tank/scratch/mfiruleva/scn/data/GSE154692/SRS7041374/counts.RData',
            '/mnt/tank/scratch/mfiruleva/scn/data/GSE154692/SRS7041375/counts.RData',
            '/mnt/tank/scratch/mfiruleva/scn/data/GSE154692/SRS7041376/counts.RData',
            '/mnt/tank/scratch/mfiruleva/scn/data/GSE154817/SRS7052767/counts.RData',
            '/mnt/tank/scratch/mfiruleva/scn/data/GSE154817/SRS7052768/counts.RData']

    output: h5ad=temp("athero_merged.h5ad"), rda="athero_merged.RData"
    benchmark: "benchmarks/analysis.txt"
    log: "logs/seurat.log"
    conda: "/mnt/tank/scratch/mfiruleva/scn/config/5d179225.yaml"
    singularity: "docker://mfiruleva/scn:latest"
    shell: "Rscript scripts/analysis.R 2> {log}"
