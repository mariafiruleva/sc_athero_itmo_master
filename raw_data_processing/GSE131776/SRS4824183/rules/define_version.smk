configfile: "/mnt/tank/scratch/mfiruleva/scn/config/config_mus.yaml"


index = config['index']
transcripts_to_genes = config['transcripts_to_genes']


rule define_version:
    """
    Define version of 10x chemistry used for particular dataset. If dataset corresponds
    to 10xv1, then downstream analysis won't be started because kallisto needs
    index file for 10xv1 platform.
    If length of the reads doesn't correspond to any 10x versions, then analysis also
    won't be pushed, and dataset must be validated manually.
    If length of the reads corresponds to particular 10x version, then defined version (and
    its whitelist) will be used in pseudoalignment step.
    """
    input: reads=expand("{accession}_1.fastq.gz", accession=barcodes.keys())
    output: kallisto_script="kallisto.sh"
    conda: "/mnt/tank/scratch/mfiruleva/scn/config/65ef9760.yaml"
    log: "logs/define_version.log"
    benchmark: "benchmarks/define_version.txt"
    params: white_10xv2="/files/10xv2_whitelist.txt", white_10xv3="/files/10xv3_whitelist.txt"
    threads: 4
    shell: "python scripts/define_version.py --fq_r1 {input.reads[0]} --threads {threads} --index {index} --transcripts_to_genes {transcripts_to_genes} --white_10xv2 {params.white_10xv2} --white_10xv3 {params.white_10xv3}"

