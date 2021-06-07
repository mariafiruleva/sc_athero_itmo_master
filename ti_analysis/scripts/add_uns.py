import argparse
import os.path
import re
from urllib.request import urlopen, Request
from xml.etree.ElementTree import parse

import numpy as np
import pandas as pd
import scanpy

def get_markers(markers_file: str) -> dict:
    markers = pd.read_csv(markers_file, sep='\t')
    return {k: np.array(list(v.values())) for k, v in markers.to_dict().items()}

def add_uns(h5: str, h5_out: str) -> None:
    file = scanpy.read_h5ad(h5)
    description = pd.read_csv(s_d).reset_index().to_dict("records")[0]
    file.uns["expType"] = "counts"
    file.uns["public"] = False
    file.uns["curated"] = False
    file.uns["token"] = 'athero_merged_icohv1Oo'
    file.uns['processed_from_panglao'] = False
    file.uns['markers'] = dict()
    resolutions = re.sub('\s', '', "0.2, 0.4, 0.6, 0.8, 1").split(',')
    for res in resolutions:
        if os.path.exists(f'markers/integrated_snn_res.{res}/markers.tsv'):
            file.uns['markers'][f'markers{res}'] = get_markers(f'markers/integrated_snn_res.{res}/markers.tsv')
    file.write_h5ad(h5_out, compression='gzip')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Define 10x version")
    parser.add_argument('--h5', type=str, required=True,
                        help='h5 input filename without uns after Seurat processing')
    parser.add_argument('--h5_out', type=str, required=True,
                        help='h5 output filename with filled uns')
    args = parser.parse_args()
    add_uns(h5=args.h5, h5_out=args.h5_out)
