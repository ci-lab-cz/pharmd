#!/usr/bin/env python3

import argparse
import os
import sys
import mdtraj as md
from multiprocessing import Pool, cpu_count
from functools import partial
from mdtraj.core.residue_names import _SOLVENT_TYPES, _WATER_RESIDUES

from pmapper.customize import load_smarts
from pmapper.pharmacophore import Pharmacophore
from io import StringIO
from collections import defaultdict
from CGRtools import PDBRead, to_rdkit_molecule
from CGRtools.algorithms.pharmacophore import (PiStack, PiCation, MetalComplex, Hydrophobic, CationPi, Salts,
                                               HydrogenAcceptor, HydrogenDonor, HalogenDonor, HalogenAcceptor, distance)


def create_parser():
    parser = argparse.ArgumentParser(description='Extract pharmacophore models from an MD trajectory of a '
                                                 'protein-ligand complex.')
    parser.add_argument('-i', '--input', metavar='input.xtc', required=False, type=str,
                        help='Input file with MD trajectory. Formats are the same as MDTraj supports.')
    parser.add_argument('-t', '--topology', metavar='input.pdb', required=False, type=str,
                        help='PDB file with topology')
    parser.add_argument('-f', '--first', metavar='INTEGER', required=False, type=int,
                        help='Staring frame number.')
    parser.add_argument('-l', '--last', metavar='INTEGER', required=False, type=int,
                        help='Last frame number.')
    parser.add_argument('-s', '--stride', metavar='INTEGER', required=False, type=int,
                        help='Step using to extract frames.')
    parser.add_argument('-o', '--output', metavar='output.pdb', required=False,
                        help='output PDB file with all extracted frames. Solvent will be omitted.')
    parser.add_argument('-p', '--pdb_input', metavar='frames.pdb', required=False, type=str, default=None,
                        help='PDB file with multiple frames of a trajectory - an alternative input to MD trajectory. '
                             'Output pharmacophore models will be stored in the same directory; '
                             'output argument would be ignored as all other arguments related to extraction of frames '
                             'from MD trajectory..')
    parser.add_argument('-g', '--lig_id', metavar='STRING', required=True, type=str,
                        help='three-letter ligand ID')
    parser.add_argument('-c', '--ncpu', metavar='INTEGER', required=False, type=int, default=1,
                        help='number of CPU to generate pharmacophores. Default: 1.')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='print progress to STDERR.')
    return parser


legend = {'a': 'aromatic_centers', 'D': 'hydrogen_donor_centers', 'A': 'hydrogen_acceptor_centers',
          'P': 'positive_charged_centers', 'N': 'negative_charged_centers', 'H': 'hydrophobic_centers',
          'M': 'metal_ligands_centers'}
reverse_legend = {Hydrophobic: 'H', PiCation: 'a', PiStack: 'a', CationPi: 'P', MetalComplex: 'M',
                  HydrogenDonor: 'D', HydrogenAcceptor: 'A', HalogenAcceptor: None, HalogenDonor: None}
smarts = load_smarts(load_metal_chelators=True)


def fix_disulphide(mol, dist=4):
    xyz = mol._conformers[0]
    bonds = mol._bonds
    neighbors = mol._neighbors
    charges = mol._charges

    sulph = []
    for n, a in mol.atoms():
        if a.atomic_number == 16 and charges[n] == -1 and neighbors[n] == 1:
            sulph.append(n)

    found = 0
    while len(sulph) > 1:
        n = sulph.pop()
        nc = xyz[n]
        try:
            m = next(m for m in sulph if distance(xyz[m], nc) < dist)
        except StopIteration:
            continue

        sulph.remove(m)
        bonds[n][m] = bonds[m][n] = Bond(1)
        neighbors[n] = neighbors[m] = 2
        charges[n] = charges[m] = 0
        found += 2

    if found:
        # fix C[S+]C
        for n, a in mol.atoms():
            if a.atomic_number == 16 and charges[n] == 1 and len(bonds[n]) == 2:
                charges[n] = 0
                found -= 1

        return not found
    return True


def get_pharmacophore(pdb, ligand, *, radius_multiplier=1.35, fix_distance=4):
    with PDBRead(StringIO(pdb), radius_multiplier=radius_multiplier, ignore=True, element_name_priority=True) as f:
        cmol = next(f)

    if not fix_disulphide(cmol, fix_distance):
        print('Charge not balanced')

    cmol.thiele()  # aromatize

    lig = [n for n, r in cmol.meta['RESIDUE'].items() if r == ligand]
    prt = [n for n, r in cmol.meta['RESIDUE'].items() if r != ligand]

    lig = cmol.substructure(lig)
    prt = cmol.substructure(prt)

    rdkit_lig = to_rdkit_molecule(lig)

    features = Pharmacophore._get_features_atom_ids(rdkit_lig, smarts)

    # apply pmapper ff to ligand
    mapping = dict(enumerate(lig))
    reverse_mapping = {v: k for k, v in mapping.items()}
    for t, ids in features.items():
        if t == 'D':
            lig.__dict__[legend[t]] = tuple((mapping[x], n) for x in ids for x in x
                                            for n in lig._bonds[mapping[x]] if lig._atoms[n].atomic_number == 1)
        elif t == 'a':
            lig.__dict__[legend[t]] = tuple(tuple(mapping[x] for x in x) for x in ids)
        else:
            lig.__dict__[legend[t]] = tuple(mapping[x] for x in ids for x in x)
    lig.__dict__['halogen_acceptor_centers'] = lig.__dict__['halogen_donor_centers'] = ()

    # mapping of atoms into groups of ff
    features_mapping = {}
    for k, v in features.items():
        if k != 'a':
            features_mapping[k] = tmp = defaultdict(list)
            for x in v:
                x = tuple(mapping[x] for x in x)
                for y in x:
                    tmp[y].append(x)

    # find contacts
    found_contacts = lig.find_contacts(prt)

    # extract ligand atoms and groups
    contacts = {'a': set(), 'D': set(), 'A': set(), 'P': set(), 'N': set(), 'H': set(), 'M': set(), None: set()}
    for c in found_contacts:
        if isinstance(c, Salts):
            if c.source in lig.positive_charged_centers:
                contacts['P'].add(c.source)
            else:
                contacts['N'].add(c.source)
        else:
            contacts[reverse_legend[type(c)]].add(c.source)

    del contacts[None]
    for k in list(contacts):
        if k != 'a':
            contacts[k] = tuple({x for x in contacts[k] for x in features_mapping[k][x]})
        else:
            contacts[k] = tuple(contacts[k])

    rdkit_contacts = {k: tuple(tuple(reverse_mapping[x] for x in x) for x in v) for k, v in contacts.items()}

    p = Pharmacophore()
    p.load_from_atom_ids(rdkit_lig, rdkit_contacts, confId=1)
    return p


def read_pdb_models(fname):
    n_returns = 0
    with open(fname) as f:
        lines = []
        for line in f:
            if not line.startswith("MODEL") and not line.startswith("ENDMDL"):
                lines.append(line)
            if line.startswith("ENDMDL"):
                lines.append(line)
                yield ''.join(lines)
                n_returns += 1
                lines = []
        if n_returns == 0:
            yield ''.join(lines)


def entry_point():

    parser = create_parser()
    args = parser.parse_args()
    args.lig_id = args.lig_id.upper()

    if args.output:
        os.makedirs(os.path.dirname(args.output), exist_ok=True)

    if args.ncpu > 1:
        pool = Pool(max(min(args.ncpu, cpu_count()), 1))
    else:
        pool = None

    non_water = _SOLVENT_TYPES - _WATER_RESIDUES

    if args.pdb_input is None:
        if args.verbose:
            sys.stderr.write('input trajectory is converting to pdb format\n')
        traj = md.load(args.input, top=args.topology, standard_names=True)
        traj.remove_solvent(inplace=True, exclude=non_water)
        traj[args.first:args.last:args.stride].save_pdb(args.output)
        pdb_input = args.output
    else:
        pdb_input = args.pdb_input

    if pool:
        for i, p in enumerate(pool.imap(partial(get_pharmacophore, ligand=args.lig_id),
                                                   read_pdb_models(pdb_input))):
            p.save_to_xyz(os.path.splitext(pdb_input)[0] + f'{i:05d}.xyz')
            if args.verbose:
                sys.stderr.write(f'\r{i + 1} pharmacophores were retrieved')
    else:
        for i, pdb_string in enumerate(read_pdb_models(pdb_input)):
            p = get_pharmacophore(pdb_string, args.lig_id)
            p.save_to_xyz(os.path.splitext(pdb_input)[0] + f'{i:05d}.xyz')
            if args.verbose:
                sys.stderr.write(f'\r{i + 1} pharmacophores were retrieved')
    if args.verbose:
        sys.stderr.write('\n')


if __name__ == '__main__':
    entry_point()
