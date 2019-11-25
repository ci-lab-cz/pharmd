#!/usr/bin/env python3
# author          : Pavel Polishchuk

import argparse
import os
import shutil
import numpy as np
import pandas as pd
from pmapper.pharmacophore import Pharmacophore as P

from scipy.spatial.distance import pdist


def create_parser():
    parser = argparse.ArgumentParser(description='Search for xyz files in a directory, calc 3D pharmacophore hashes '
                                                 'and save models with distinct hashes to an output directory.')
    parser.add_argument('-i', '--input', metavar='DIR_NAME', required=True, type=str,
                        help='directory with xyz file of pharmacophore models.')
    parser.add_argument('-s', '--bin_step', metavar='NUMERIC', required=False, default=1, type=float,
                        help='bin step. Default: 1.')
    parser.add_argument('-o', '--output', metavar='DIR_NAME', required=True, type=str,
                        help='xyz files having distinct 3D pharmacophore hashes will be stored to the specified '
                             'directory.')
    parser.add_argument('-d', '--output_hashes', metavar='hashes.txt', required=False, type=str, default=None,
                        help='text file with 3D pharmacophore hashes and auxiliary information for all input '
                             'pharmacophores. If not specified the file named hashes.txt will be created in the '
                             'input directory.')
    return parser


def entry_point():
    parser = create_parser()
    args = parser.parse_args()

    if args.output_hashes is None:
        args.output_hashes = os.path.join(os.path.dirname(args.input), 'hashes.txt')

    files = [f for f in os.listdir(args.input) if f.endswith(".xyz")]

    data = []
    for fname in sorted(files):
        p = P(bin_step=args.bin_step)
        p.load_from_xyz(os.path.join(args.input, fname))
        h = p.get_signature_md5()
        features = p.get_feature_coords()
        fstr = ''.join(sorted(item[0] for item in features))
        fcount = len(features)
        fcount_uniq = len(set(item[1] for item in features))
        dist = pdist(np.array(list(item[1] for item in features)), metric='euclidean')
        try:
            max_dist = round(max(dist), 1)
        except ValueError:
            max_dist = 0
        data.append([fname.replace('.xyz', ''), h, fcount, fcount_uniq, fstr, max_dist])

    df = pd.DataFrame(data, columns=['filename', 'hash', 'count', 'ucount', 'features', 'max_dist'])
    df.to_csv(args.output_hashes, sep='\t', index=False)

    df.drop_duplicates('hash', inplace=True)
    os.makedirs(args.output, exist_ok=True)
    for fname in df['filename']:
        shutil.copyfile(os.path.join(args.input, fname + '.xyz'), os.path.join(args.output, fname + '.xyz'))


if __name__ == '__main__':
    entry_point()
