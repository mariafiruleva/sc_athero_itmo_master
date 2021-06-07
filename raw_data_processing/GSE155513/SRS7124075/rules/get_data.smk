
rule fastq_dump:
    output: "{accession}_sra.txt"
    params: run_id=lambda wildcards: wildcards.accession
    conda: "/mnt/tank/scratch/mfiruleva/scn/config/65ef9760.yaml"
    threads: 1
    shell:
         """
         esearch -db sra -query {params.run_id} | efetch -format metadata | grep -Po 'average="(\d+)"' | awk '!x[$0]++' | sed 's/average=//g' | sed 's/"//g' > {params.run_id}_sra.txt
         """

