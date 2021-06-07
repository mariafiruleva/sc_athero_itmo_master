

rule get_bam_header:
    output: temp("{accession}_tmp.bam")
    conda: "/mnt/tank/scratch/mfiruleva/scn/config/65ef9760.yaml"
    singularity: "docker://mfiruleva/scn:latest"
    log: "logs/get_bam_{accession}.log"
    benchmark: "benchmarks/get_bam_{accession}.txt"
    params: link=lambda wildcards, output: bam_files[f'{wildcards.accession}']
    shell: "wget -q -O - {params.link} | head -1000 > {output} || true"

