import argparse
import re
from urllib.request import urlopen, Request
from xml.etree.ElementTree import parse

import numpy as np
import pandas as pd
import scanpy

def get_markers(markers_file: str) -> dict:
    markers = pd.read_csv(markers_file, sep='\t')
    return {k: np.array(list(v.values())) for k, v in markers.to_dict().items()}

path = "/mnt/tank/scratch/mfiruleva/scn/data/GSE131776/SRS4824187"


def add_uns(h5: str, h5_out: str, s_d: str, summary_file: str, kallisto_script=None, technology=None) -> None:
    file = scanpy.read_h5ad(h5)

    
    description = pd.read_csv(s_d).reset_index().to_dict("records")[0]
    file.uns["expType"] = "counts"
    file.uns["public"] = True
    file.uns["curated"] = False
    file.uns["gse"] = description['GSE']
    file.uns["geo"] = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={description['GSE']}"
    file.uns["study_accession"] = description['study_accession']
    file.uns["species"] = description['scientific_name']
    if isinstance(file.uns["species"], list):
        file.uns["species"] = file.uns["species"][0]
    if technology:
        file.uns['technology'] = technology
    else:
        if description['technology'] != "10x":
            file.uns["technology"] = description['technology']
        else:
            with open(kallisto_script, 'r') as run_file:
                data = run_file.read().replace('\n', '')
            file.uns["technology"] = re.findall('10xv[0-9]*', data)[0]

    link = Request(f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={description["GSE"]}',
                   headers={'User-Agent': 'Mozilla/5.0'})
    link = urlopen(link)
    article = Request(link.url, headers={'User-Agent': 'Mozilla/5.0'})
    response = urlopen(article).read()
    acc_ids = {'SRP': re.findall('SRP\d*', response.decode('utf-8'))[0],
               'PRJNA': re.findall('SRP\d*', response.decode('utf-8'))[0]
               }
    if acc_ids['SRP']:
        var_url = urlopen(
            f'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runtable&term={acc_ids["SRP"]}')
    else:
        var_url = urlopen(
            f'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runtable&term={acc_ids["PRJNA"]}')
    xmldoc = parse(var_url)
    file.uns["title"] = xmldoc.findall('EXPERIMENT_PACKAGE/STUDY/DESCRIPTOR/STUDY_TITLE')[0].text
    study_des = xmldoc.findall('EXPERIMENT_PACKAGE/STUDY/DESCRIPTOR/STUDY_ABSTRACT')[0].text
    file.uns["description"] = re.sub('Overall design:\s*', '', study_des)
    file.uns["design"] = re.sub('Overall design:\s*', '', re.findall('Overall design:.*', study_des)[0])


    
    file.uns["token"] = description['secondary_sample_accession']
    file.uns["sra"] = f"https://www.ncbi.nlm.nih.gov/sra/{description['secondary_sample_accession']}"

    
    
    file.uns['processed_from_panglao'] = False

    
    meta = {'dataset': file.uns['gse'], 'sample': file.uns['token'], 'organism': file.uns['species'],
            'technology': file.uns['technology'], 'path': path}
    pd.DataFrame.from_dict(meta, orient='index').T.to_csv(summary_file, mode='a', header=False, index=False)


    
    file.uns['markers'] = dict()
    resolutions = re.sub('\s', '', "0.2, 0.4, 0.6, 0.8, 1").split(',')
    for res in resolutions:
        file.uns['markers'][f'markers{res}'] = get_markers(f'markers/SCT_snn_res.{res}/markers.tsv')

    file.write_h5ad(h5_out, compression='gzip')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Define 10x version")
    parser.add_argument('--h5', type=str, required=True,
                        help='h5 input filename without uns after Seurat processing')
    parser.add_argument('--h5_out', type=str, required=True,
                        help='h5 output filename with filled uns')
    parser.add_argument('--kallisto_script', type=str, required=False, default=None,
                        help='Path to kallisto script')
    parser.add_argument('--s_d', type=str, required=True,
                        help='Path to sample description file')
    parser.add_argument('--summary_file', type=str, required=True,
                        help='Path to the summary file')
    parser.add_argument('--technology', type=str, required=False, default=None,
                        help='Name of used technology; this argument specified in case of panglao db')
    args = parser.parse_args()

    add_uns(h5=args.h5, h5_out=args.h5_out, kallisto_script=args.kallisto_script,
            s_d=args.s_d, summary_file=args.summary_file, technology=args.technology)