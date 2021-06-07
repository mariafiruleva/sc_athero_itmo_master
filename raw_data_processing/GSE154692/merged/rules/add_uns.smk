rule add_uns_to_h5ad:
    """
    Add uns information to h5ad object after Seurat processing
    """
    input: rules.run_analysis.output.h5ad
    output: "GSE154692_with_uns.h5ad"
    conda: "/mnt/tank/scratch/mfiruleva/scn/config/5d179225.yaml"
    params: kallisto="/mnt/tank/scratch/mfiruleva/scn/data/GSE154692/SRS7041373/kallisto.sh", s_d="/mnt/tank/scratch/mfiruleva/scn/data/GSE154692/SRS7041373/sample_description.csv", summary="/mnt/tank/scratch/mfiruleva/scn/stats/summary.csv"

    shell:
         "python scripts/add_uns.py --h5 {input} --h5_out {output} --kallisto_script {params.kallisto} \
         --s_d {params.s_d} --summary_file {params.summary}"
