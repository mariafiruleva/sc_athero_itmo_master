rule add_uns_to_h5ad:
    """
    Add uns information to h5ad object after Seurat processing
    """
    input: rules.run_analysis.output.h5ad
    output: "athero_merged_with_uns.h5ad"
    conda: "/mnt/tank/scratch/mfiruleva/scn/config/5d179225.yaml"
    shell:
         "python scripts/add_uns.py --h5 {input} --h5_out {output}"
