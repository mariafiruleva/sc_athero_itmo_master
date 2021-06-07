
rule download_fq_header:
    """
    Download the first read for forward fastq file
    and gzipped it in to R1.gz file. Pipeline uses this file for definition of 10x version:
    """
    output: temp("{accession}_1.fastq.gz")
    conda: "/mnt/tank/scratch/mfiruleva/scn/config/65ef9760.yaml"
    log: "logs/prepare_loading_{accession}.log"
    benchmark: "benchmarks/prepare_loading_{accession}.txt"
    params: link=lambda wildcards, output: barcodes[f'{wildcards.accession}']
    shell: "wget -q -O - ftp://{params.link}| zcat | head -400000 | gzip > {output} || true"

