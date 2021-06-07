import argparse
import os.path
import pandas as pd
import re


import subprocess

read_info = dict()
read_info['technology'] = 'NA'


def def_type(read_len: int, postfix: str, whitelist_10xv2: str, whitelist_10xv3: str) -> None:
    if read_len < 12:
        read_info['index'] = postfix
    if read_len == 26:
        read_info['barcode_len'] = read_len
        read_info['barcode'] = postfix
        read_info['technology'] = '10xv2'
        read_info['whitelist'] = whitelist_10xv2
    if read_len == 28:
        read_info['barcode_len'] = read_len
        read_info['barcode'] = postfix
        read_info['whitelist'] = whitelist_10xv3
        read_info['technology'] = '10xv3'
    if read_len > 50:
        read_info['bio'] = postfix


def get_info(sra_file, threads: int, index: str, transcripts_to_genes: str,
             whitelist_10xv2: str, whitelist_10xv3: str, kallisto_script='kallisto.sh') -> None:
    number_of_files = int(subprocess.getoutput(f"wc -l {sra_file}")[0])
    sra = pd.read_csv(f'{sra_file}', names=['length']).transpose()
    sra.columns = range(1, number_of_files + 1)
    for postfix in range(1, number_of_files + 1):
        def_type(sra[postfix][0], postfix=postfix, whitelist_10xv2=whitelist_10xv2,
                 whitelist_10xv3=whitelist_10xv3)
    if read_info['technology'] == 'NA':
        raise Exception("Technology wasn't defined")
    with open(kallisto_script, 'w') as out_file:
        out_file.write("mkfifo R1.gz R2.gz\n")
        out_file.write(f"cat *{read_info['barcode']}.fastq.gz > R1.gz &\n")
        out_file.write(f"cat *{read_info['bio']}.fastq.gz > R2.gz &\n")
        out_file.write(
            f'kallisto bus -i {index} -x {read_info["technology"]} -t {threads} -o bus_out/ R1.gz R2.gz \n')
        out_file.write('mkdir bus_out/tmp\n')
        out_file.write(
            f'bustools correct -w {read_info["whitelist"]} -o bus_out/correct_output.bus bus_out/output.bus\n')
        out_file.write(
            f'bustools sort -t {threads} -T bus_out/tmp/ -p bus_out/correct_output.bus | bustools count -o bus_out/genes -g {transcripts_to_genes} -e bus_out/matrix.ec -t bus_out/transcripts.txt --genecounts -')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Define 10x version")
    parser.add_argument('--sra_file', type=str, required=True,
                        help='File with average read length per each fastq file for corresponding run')
    parser.add_argument('--threads', type=int, required=True,
                        help='Number of threads for parallel-fastq-dump program')
    parser.add_argument('--index', type=str, required=True,
                        help='Path to kallisto index')
    parser.add_argument('--transcripts_to_genes', type=str, required=True,
                        help='Path to transcripts_to_genes file')
    parser.add_argument('--white_10xv2', type=str, required=True,
                        help='Path to whitelist for 10xv2')
    parser.add_argument('--white_10xv3', type=str, required=True,
                        help='Path to whitelist for 10xv3')
    args = parser.parse_args()
    get_info(args.sra_file, args.threads, args.index, args.transcripts_to_genes, args.white_10xv2,
             args.white_10xv3)

