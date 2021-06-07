configfile: "/mnt/tank/scratch/mfiruleva/scn/config/config_mus.yaml"


index = config['index']
transcripts_to_genes = config['transcripts_to_genes']


rule define_version:
    """
    Define version of 10x chemistry used for particular dataset using whitelists and CB:Z tag in the
    first 1e3 reads from bam file.
    In the output BAM file, the original uncorrected barcode is encoded in the CR tag, and the corrected
    barcode sequence is encoded in the CB tag.
    Reads that are not able to be assigned a corrected barcode will not have a CB tag.
    Source:
    https://kb.10xgenomics.com/hc/en-us/articles/115003822406-How-does-Cell-Ranger-correct-barcode-sequencing-errors-
    """
    input: bams=expand("{accession}_tmp.bam", accession=bam_files.keys()), s_d="sample_description.csv"
    output: kallisto_script="kallisto.sh"
    conda: "/mnt/tank/scratch/mfiruleva/scn/config/65ef9760.yaml"
    log: "logs/define_version.log"
    singularity: "docker://mfiruleva/scn:latest"
    benchmark: "benchmarks/define_version.txt"
    params: white_10xv1="/files/10xv1_whitelist.txt", white_10xv2="/files/10xv2_whitelist.txt",
            white_10xv3="/files/10xv3_whitelist.txt", threads = 4
    shell:
        """
        python scripts/define_version.py --s_d {input.s_d} --tmp_bam {input.bams[0]} \
        --threads {params.threads} --index {index} --transcripts_to_genes {transcripts_to_genes} \
        --white_10xv1 {params.white_10xv1} --white_10xv2 {params.white_10xv2} --white_10xv3 {params.white_10xv3}
        """

