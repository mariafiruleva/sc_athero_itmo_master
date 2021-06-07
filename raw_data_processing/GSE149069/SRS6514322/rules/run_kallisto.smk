
rule kallisto:
    input: rules.define_version.output.kallisto_script
    output: temp("R1.gz"), temp("R2.gz"), out=temp("bus_out/correct_output.bus"), uncor_out=temp("bus_out/output.bus")
    log: "logs/kallisto.log"
    benchmark: "benchmarks/kallisto.txt"
    conda: "/mnt/tank/scratch/mfiruleva/scn/config/65ef9760.yaml"
    shell: "bash kallisto.sh 2> {log}"