#!/usr/bin/env python3
# author          : Pavel Polishchuk

import argparse
import sys
import os
import pandas as pd
from psearch.database import DB


def get_conf_count(db_fname):
    d = DB(db_fname)
    mol_names = d.get_mol_names()
    data = []
    for mol_name in mol_names:
        n = d.get_conf_count(mol_name)
        data.append((mol_name, n))
    data = pd.DataFrame(data, columns=['mol_id', 'total_conf_count'])
    return data


def create_parser():
    parser = argparse.ArgumentParser(description='Calculate scores from multiple hit lists returned by pharmacophore '
                                                 'model screening.')
    parser.add_argument('-i', '--input', metavar='DIRNAME', required=True, type=str,
                        help='directory with text files containing hit lists for scrrening with individual models. '
                             'All files having txt extension will be considered as hit lists.')
    parser.add_argument('-o', '--output', metavar='output.txt', required=True, type=str,
                        help='text file with calculated scoring.')
    parser.add_argument('-s', '--score', metavar='SCORE', required=False, type=str, default='cca',
                        help='scheme using for calculation of scores. Can be "cca" (conformers coverage appraoch) or '
                             '"cha" (common hits approach). Default: cca.')
    parser.add_argument('-d', '--db_name', metavar='conformers.db', required=False, type=str, default=None,
                        help='DB with precalculated conformers used for screening. It is only required to calculate '
                             '"cca" scores.')
    return parser


def entry_point():
    parser = create_parser()
    args = parser.parse_args()

    args.score = args.score.lower()
    if args.score == 'cca' and args.db_name is None:
        sys.stderr.write('Database should be supplied to calculate cca scores.')
        exit()

    files = [f for f in os.listdir(args.input) if os.path.isfile(os.path.join(args.input, f)) and f.endswith('.txt')]
    if not files:
        sys.stderr.write('There are no txt-files in the input directory.')
        exit()

    if args.score == 'cca':

        conf_count = get_conf_count(args.db_name)

        df = pd.read_csv(os.path.join(args.input, files[0]), sep='\t', header=None, names=['mol_id', 'stereo_id', 'conf_id'])
        for f in files[1:]:
            d = pd.read_csv(os.path.join(args.input, f), sep='\t', header=None, names=['mol_id', 'stereo_id', 'conf_id'])
            df = df.append(d, ignore_index=True)
            df.drop_duplicates(inplace=True)
        df = df.groupby('mol_id').agg('count').reset_index()
        df = df.merge(conf_count, how='left', on='mol_id')
        df['cca_score'] = round(df['conf_id'] / df['total_conf_count'], 3)
        df.to_csv(args.output, sep='\t', columns=['mol_id', 'cca_score'], index=False)

    elif args.score == 'cha':

        df = pd.read_csv(os.path.join(args.input, files[0]), sep='\t', header=None, names=['mol_id', 'stereo_id', 'conf_id'], usecols=['mol_id'])
        df.drop_duplicates(inplace=True)
        for i, f in enumerate(files[1:], 1):
            d = pd.read_csv(os.path.join(args.input, f), sep='\t', header=None, names=['mol_id', 'stereo_id', 'conf_id'], usecols=['mol_id'])
            d.drop_duplicates(inplace=True)
            df = df.append(d, ignore_index=True)
        df = df.reset_index().groupby('mol_id').agg('count').reset_index()
        df['cha_score'] = round(df['index'] / len(files), 3)
        df.to_csv(args.output, sep='\t', columns=['mol_id', 'cha_score'], index=False)


if __name__ == '__main__':
    entry_point()
