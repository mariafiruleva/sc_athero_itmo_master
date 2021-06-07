



rule kallisto:
    input: expand("{accession}_1.fastq.gz", accession=runs), expand("{accession}_2.fastq.gz", accession=runs)
    output: temp("R1.gz"), temp("R2.gz"), out=temp("bus_out/correct_output.bus"), uncor_out=temp("bus_out/output.bus")
    log: "logs/kallisto.log"
    benchmark: "benchmarks/kallisto.txt"
    conda: "/mnt/tank/scratch/mfiruleva/scn/config/65ef9760.yaml"
    singularity: "docker://mfiruleva/scn:latest"
    threads: 4
    shell: "bash kallisto.sh 2> {log}"
