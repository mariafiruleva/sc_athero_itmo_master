import argparse
import os.path
import pandas as pd
import re
import gzip

def tech_version(fq_r1: str, whielist_10xv2: str, whielist_10xv3: str) -> dict:
    res = {'barcode_len': '', 'technology': 'NA', 'whitelist': 'NA'}
    with gzip.open(fq_r1, 'r') as in_file:
        for idx, line in enumerate(in_file, start=1):
            if idx == 2:
                barcode_len = len(line.decode('ascii').strip())
                res['barcode_len'] = barcode_len
                if barcode_len == 26:
                    res['technology'] = "10xv2"
                    res['whitelist'] = whielist_10xv2
                elif barcode_len == 28:
                    res['technology'] = "10xv3"
                    res['whitelist'] = whielist_10xv3
                else:
                    res['technology'] = "NA"
    if res['technology'] == "10xv3":
        tt_count = 0
        nn_count = 0
        with gzip.open(fq_r1, 'r') as in_file:
            for idx, line in enumerate(in_file, start=1):
                if idx > 1e5:
                    break
                if idx and not idx % 2 and idx % 4:
                    if line.decode('ascii')[26:28] == 'TT':
                        tt_count += 1
                    if line.decode('ascii')[26:28] == 'NN':
                        nn_count += 1
        if ((tt_count + nn_count) / 1e5) > 0.5:  # 4e5 / 4 = 1e5
            res['technology'] = "NA"
    return res


def get_info(fq_r1: str, threads: int, index: str, transcripts_to_genes: str,
             whielist_10xv2: str, whielist_10xv3: str, kallisto_script='kallisto.sh'):
    description = pd.read_csv('sample_description.csv').reset_index().to_dict("records")
    r1 = [re.sub('ftp://', '', x['fastq_ftp'].split(';')[0]) for x in description]
    r2 = [re.sub('ftp://', '', x['fastq_ftp'].split(';')[1]) for x in description]
    read_info = tech_version(fq_r1, whielist_10xv2, whielist_10xv3)
    if read_info['technology'] == 'NA':
        raise Exception("Technology wasn't defined")
    if not os.path.isfile(kallisto_script):
        with open(kallisto_script, 'w') as out_file:
            out_file.write("mkfifo R1.gz R2.gz\n")
            out_file.write("curl -Ls ftp://{" + ','.join(r1) + "} > R1.gz &\n")
            out_file.write("curl -Ls ftp://{" + ','.join(r2) + "} > R2.gz &\n")
            out_file.write(
                f'kallisto bus -i {index} -x {read_info["technology"]} -t {threads} -o bus_out/ R1.gz R2.gz \n')
            out_file.write('mkdir bus_out/tmp\n')
            out_file.write(
                f'bustools correct -w {read_info["whitelist"]} -o bus_out/correct_output.bus bus_out/output.bus\n')
            out_file.write(
                f'bustools sort -t {threads} -T bus_out/tmp/ -p bus_out/correct_output.bus | bustools count -o bus_out/genes -g {transcripts_to_genes} -e bus_out/matrix.ec -t bus_out/transcripts.txt --genecounts -')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Define 10x version")
    parser.add_argument('--fq_r1', type=str, required=True,
                        help='fastq file header')
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
    get_info(args.fq_r1, args.threads, args.index, args.transcripts_to_genes, args.white_10xv2, args.white_10xv3)

