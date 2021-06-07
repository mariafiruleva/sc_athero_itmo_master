
def is_merged():
    if len(bam_files.keys()) > 1:
        return temp('merged.bam')
    else:
        return []

def get_tmp_bams():
    return [temp(f'file{idx}.bam') for idx in range(0, len(bam_files.keys()))]





rule kallisto:
    input: rules.define_version.output.kallisto_script
    output: is_merged(), get_tmp_bams(),
            out=temp("bus_out/correct_output.bus"), uncor_out=temp("bus_out/output.bus")
    log: "logs/kallisto.log"
    benchmark: "benchmarks/kallisto.txt"
    conda: "/mnt/tank/scratch/mfiruleva/scn/config/65ef9760.yaml"
    singularity: "docker://mfiruleva/scn:latest"
    threads: 4
    shell: "bash kallisto.sh 2> {log}"
