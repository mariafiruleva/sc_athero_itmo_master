configfile: "/mnt/tank/scratch/mfiruleva/scn/config/config_mus.yaml"


index = config['index']
transcripts_to_genes = config['transcripts_to_genes']


rule define_version:
    input: "{accession}_sra.txt"
    output: temp("{accession}_1.fastq.gz"), temp("{accession}_2.fastq.gz"), sra=temp("sra/{accession}.sra")
    log: fq_dump="logs/fq_dump/{accession}_dump.log"
    conda: "/mnt/tank/scratch/mfiruleva/scn/config/65ef9760.yaml"
    params: run_id=lambda wildcards: wildcards.accession,
            white_10xv2="/files/10xv2_whitelist.txt", white_10xv3="/files/10xv3_whitelist.txt", max_size='50G'
    threads: 4
    shell:
         """
         python scripts/define_version.py --sra_file {input} --threads {threads} --index {index} --transcripts_to_genes {transcripts_to_genes} --white_10xv2 {params.white_10xv2} --white_10xv3 {params.white_10xv3}
         prefetch {params.run_id} -o {output.sra} --max-size {params.max_size}
         parallel-fastq-dump -s {output.sra} --split-files --threads {threads} -O . --tmpdir . --gzip
         """

