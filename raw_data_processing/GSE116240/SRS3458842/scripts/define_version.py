import argparse
import os.path
import pandas as pd
import re


import argparse
import re
import warnings

warnings.simplefilter("ignore", UserWarning)
import numpy as np
import operator
import pysam
import pandas as pd


def get_info(s_d: str, tmp_bam: str, index: str, thread: str, transcripts_to_genes: str,
             white_10xv1: str, white_10xv2: str, white_10xv3: str, kallisto_script='kallisto.sh') -> None:
    whitelist_pathes = {"10xv1": white_10xv1,
                        "10xv2": white_10xv2,
                        "10xv3": white_10xv3}
    with open(white_10xv1) as f:
        whitelist_v1 = np.array(f.read().splitlines())
    with open(white_10xv2) as f:
        whitelist_v2 = np.array(f.read().splitlines())
    with open(white_10xv3) as f:
        whitelist_v3 = np.array(f.read().splitlines())
    with pysam.AlignmentFile(tmp_bam, "rb", ignore_truncation=True) as file:
        result = {"10xv1": 0, "10xv2": 0, "10xv3": 0}
        for idx, line in enumerate(file):
            try:
                if idx > 1e3:
                    break
                barcode = re.sub('-[0-9]*', '', line.get_tag('CB'))
                if barcode.strip() in whitelist_v1:
                    result["10xv1"] += 1
                if barcode.strip() in whitelist_v2:
                    result["10xv2"] += 1
                if barcode.strip() in whitelist_v3:
                    result["10xv3"] += 1
            except:
                continue
    max_n = max(result.items(), key=operator.itemgetter(1))[1]
    res = [x[0] for x in result.items() if x[1] == max_n]
    if not max_n or len(res) > 1:
        raise Exception("Technology wasn't defined")
    technology = res[0]
    whitelist = whitelist_pathes[technology]
    description = pd.read_csv(s_d).reset_index().to_dict("records")
    with open(kallisto_script, 'w') as out_file:
        bam = [f"ftp://{bam_file['submitted_ftp'].split(';')[0]}" for bam_file in description]
        file_names = [f'file{idx}.bam' for idx in range(0, len(bam))]
        if len(file_names) > 30:
            raise Exception(f"Too many files (there are {len(file_names)} files), break the process")
        out_file.write(f"mkfifo {' '.join(file_names)}\n")
        if len(file_names) > 1:
            out_file.write("mkfifo merged.bam\n")
        for idx, file in enumerate(file_names):
            out_file.write(f"curl -Ls {bam[idx]} > {file} &\n")
        if technology == "10xv1":
            if len(file_names) > 1:
                out_file.write(f"samtools cat -o merged.bam file* &\n")
                out_file.write("/usr/bin/bamtofastq-1.2.0 --reads-per-fastq=100000000000000 merged.bam out_bam\n")
            else:
                out_file.write(
                    f"/usr/bin/bamtofastq-1.2.0 --reads-per-fastq=100000000000000 {' '.join(file_names)} out_bam\n")
            out_file.write(f"cat out_bam/*/*R1*.gz > R1.gz\n")
            out_file.write(f"cat out_bam/*/*R2*.gz > R2.gz\n")
            out_file.write(f"cat out_bam/*/*R3*.gz > R3.gz\n")
            out_file.write(
                f'kallisto bus -i {index} -x 2,0,14:1,0,10:0,0,0 -t {thread} -o bus_out/ R1.gz R3.gz R2.gz\n')
        else:
            if len(file_names) > 1:
                out_file.write(f"samtools cat -o merged.bam file* &\n")
                out_file.write(
                    f'kallisto bus --bam -i {index} -x {technology} -t {thread} -o bus_out/ merged.bam\n')
            else:
                out_file.write(
                    f'kallisto bus --bam -i {index} -x {technology} -t {thread} -o bus_out/ {" ".join(file_names)}\n')
        out_file.write('mkdir bus_out/tmp\n')
        out_file.write(f'bustools correct -w {whitelist} -o bus_out/correct_output.bus bus_out/output.bus\n')
        out_file.write(
            f'bustools sort -t {thread} -T bus_out/tmp/ -p bus_out/correct_output.bus | bustools count -o bus_out/genes -g {transcripts_to_genes} -e bus_out/matrix.ec -t bus_out/transcripts.txt --genecounts -\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Define 10x version")
    parser.add_argument('--s_d', type=str, required=True,
                        help='Path to sample description')
    parser.add_argument('--tmp_bam', type=str, required=True,
                        help='Downloaded header of bam file')
    parser.add_argument('--threads', type=int, required=True,
                        help='Number of threads for kallisto')
    parser.add_argument('--index', type=str, required=True,
                        help='Path to kallisto index')
    parser.add_argument('--transcripts_to_genes', type=str, required=True,
                        help='Path to transcripts_to_genes file')
    parser.add_argument('--white_10xv1', type=str, required=True,
                        help='Path to whitelist for 10xv1')
    parser.add_argument('--white_10xv2', type=str, required=True,
                        help='Path to whitelist for 10xv2')
    parser.add_argument('--white_10xv3', type=str, required=True,
                        help='Path to whitelist for 10xv3')
    args = parser.parse_args()
    get_info(args.s_d, args.tmp_bam , args.index, args.threads, args.transcripts_to_genes,
             args.white_10xv1, args.white_10xv2, args.white_10xv3)
