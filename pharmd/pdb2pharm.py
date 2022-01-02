#!/usr/bin/env python3

import os
import argparse
from pharmd.md2pharm import get_pharmacophores


def create_parser():
    parser = argparse.ArgumentParser(description='Extract pharmacophore models from a PDB files of a '
                                                 'protein-ligand complex.')
    parser.add_argument('-p', '--pdbs', metavar='FILENAME(S) or DIRNAME(S)', required=True, type=str, nargs='+',
                        help='PDB file')
    parser.add_argument('-g', '--lig_id', metavar='STRING', required=True, type=str,
                        help='three-letter ligand ID')
    parser.add_argument('-o', '--output', metavar='output.pdb', required=False, default=None,
                        help='output PDB file with all extracted frames. Solvent will be omitted.')
    parser.add_argument('-c', '--ncpu', metavar='INTEGER', required=False, type=int, default=1,
                        help='number of CPU to generate pharmacophores. Default: 1.')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='print progress to STDERR.')
    return parser


def entry_point():
    parser = create_parser()
    args = parser.parse_args()
    args.lig_id = args.lig_id.upper()

    pdb_input = args.pdbs
    output = args.output if args.output else os.path.join(os.path.split(pdb_input)[0], 'pharmacophores')
    os.makedirs(output, exist_ok=True)

    if all(os.path.isdir(item) for item in pdb_input):
        for dname in pdb_input:
            dname = os.path.abspath(dname)
            for fpdb in os.listdir(dname):
                for p in get_pharmacophores(os.path.join(dname, fpdb), args.lig_id):
                    p.save_to_xyz(os.path.join(output, f'{os.path.splitext(fpdb)[0]}.xyz'))
    elif all(os.path.isfile(item) for item in pdb_input):
        for fpdb in pdb_input:
            for p in get_pharmacophores(fpdb, args.lig_id):
                p.save_to_xyz(os.path.join(f'{os.path.splitext(fpdb)[0]}.xyz'))
    else:
        raise ValueError('Input queries should be all either files or directories not a mix.')


if __name__ == '__main__':
    entry_point()
