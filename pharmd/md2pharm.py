#!/usr/bin/env python3

import argparse
import os
import sys
import mdtraj as md
from multiprocessing import Pool, cpu_count
from mdtraj.core.residue_names import _SOLVENT_TYPES, _WATER_RESIDUES

from pmapper.customize import load_smarts
from pmapper.pharmacophore import Pharmacophore
from io import TextIOWrapper
from collections import defaultdict
from CGRtools import PDBRead, to_rdkit_molecule
from PharmaContacts import (PiStack, PiCation, MetalComplex, Hydrophobic, CationPi, Salts,
                            HydrogenAcceptor, HydrogenDonor, HalogenDonor, HalogenAcceptor, PharmaContacts)
from CGRtools.files import XYZrw
from CGRtools import SDFRead
from CGRtools.files.XYZrw import get_possible_bonds
from itertools import combinations, chain
from math import sqrt
from pkg_resources import resource_stream


PDBRead.MoleculeContainer = PharmaContacts  # set molecule class to extended inherited


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
smarts = load_smarts()


with SDFRead(TextIOWrapper(resource_stream(__package__, 'queries.sdf'))) as f:
    queries = f.read()


def fix_disulphide(dist=4):
    def w(atoms, conformer, multiplier):
        possible_bonds = get_possible_bonds(atoms, conformer, multiplier)
        sulf = []
        for n, a in atoms.items():
            if a.atomic_number == 16:
                sulf.append(n)

        for n, m in combinations(sulf, 2):
            nx, ny, nz = conformer[n]
            mx, my, mz = conformer[m]
            d = sqrt((nx - mx) ** 2 + (ny - my) ** 2 + (nz - mz) ** 2)
            if d <= dist:
                possible_bonds[n][m] = possible_bonds[m][n] = d
        return possible_bonds
    XYZrw.get_possible_bonds = w


def get_pharmacophores(pdb, ligand, *, radius_multiplier=1.4):
    with PDBRead(pdb, radius_multiplier=radius_multiplier, ignore=True, element_name_priority=True,
                 parse_as_single=True) as f:
        cmol = next(f)  # load graph

        lig = [n for n, r in cmol.meta['RESIDUE'].items() if r == ligand]
        prt = [n for n, r in cmol.meta['RESIDUE'].items() if r != ligand]

        # check for ligand-protein bonding
        sprt = set(prt)
        if not all(sprt.isdisjoint(cmol._bonds[n]) for n in lig):
            raise Exception('ligand-protein bond found')

        # check charges
        valid = set()
        for q in queries:
            for m in q.get_mapping(cmol):
                valid.update(m.values())
        for n, c in cmol._charges.items():
            if n in sprt and c and n not in valid:
                raise Exception(f'invalid charges found on atom: {n} [{cmol.augmented_substructure([n], deep=5)}]')

        cmol.thiele()  # aromatize
        cgr_lig = cmol.substructure(lig)
        cgr_prt = cmol.substructure(prt)
        rdkit_lig = to_rdkit_molecule(cgr_lig)
        features = Pharmacophore._get_features_atom_ids(rdkit_lig, smarts)

        # apply pmapper pp to ligand
        mapping = dict(enumerate(cgr_lig))
        reverse_mapping = {v: k for k, v in mapping.items()}
        for t, ids in features.items():
            if t == 'D':
                cgr_lig.__dict__[legend[t]] = tuple((mapping[x], n) for x in ids for x in x
                                                    for n in cgr_lig._bonds[mapping[x]]
                                                    if cgr_lig._atoms[n].atomic_number == 1)
            elif t == 'a':
                cgr_lig.__dict__[legend[t]] = tuple(tuple(mapping[x] for x in x) for x in ids)
            else:
                cgr_lig.__dict__[legend[t]] = tuple(mapping[x] for x in ids for x in x)
        cgr_lig.__dict__['halogen_acceptor_centers'] = cgr_lig.__dict__['halogen_donor_centers'] = ()

        # mapping of atoms into groups of pp
        features_mapping = {}
        for k, v in features.items():
            if k != 'a':
                features_mapping[k] = tmp = defaultdict(list)
                for x in v:
                    x = tuple(mapping[x] for x in x)
                    for y in x:
                        tmp[y].append(x)

        # find contacts
        for conf in chain(cmol._conformers, f):
            cgr_lig._conformers[0] = {n: conf[n] for n in lig}
            cgr_prt._conformers[0] = conf

            found_contacts = cgr_lig.find_contacts(cgr_prt)

            # extract ligand atoms and groups
            contacts = {'a': set(), 'D': set(), 'A': set(), 'P': set(), 'N': set(), 'H': set(), 'M': set(), None: set()}
            for c in found_contacts:
                if isinstance(c, Salts):
                    if c.source in cgr_lig.positive_charged_centers:
                        contacts['P'].add(c.source)
                    else:
                        contacts['N'].add(c.source)
                else:
                    contacts[reverse_legend[type(c)]].add(c.source)

            del contacts[None]
            for k in list(contacts):
                # remove all`M` feature
                if k == 'M':
                    continue
                if k != 'a':
                    contacts[k] = tuple({x for x in contacts[k] for x in features_mapping[k][x]})
                else:
                    contacts[k] = tuple(contacts[k])

            rdkit_contacts = {k: tuple(tuple(reverse_mapping[x] for x in x) for x in v) for k, v in contacts.items()}

            rdkit_lig = to_rdkit_molecule(cgr_lig)
            p = Pharmacophore()
            p.load_from_atom_ids(rdkit_lig, rdkit_contacts, confId=1)
            yield p


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

    # call this function ONCE!!!
    fix_disulphide(dist=5)  # before forks!

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

    for i, p in enumerate(get_pharmacophores(pdb_input, args.lig_id)):
        p.save_to_xyz(os.path.splitext(pdb_input)[0] + f'{i:05d}.xyz')
        sys.stderr.write(f'\r{i + 1} pharmacophores were retrieved')

    # if pool:
    #     for i, p in enumerate(pool.imap(partial(get_pharmacophore, ligand=args.lig_id),
    #                                                read_pdb_models(pdb_input))):
    #         p.save_to_xyz(os.path.splitext(pdb_input)[0] + f'{i:05d}.xyz')
    #         if args.verbose:
    #             sys.stderr.write(f'\r{i + 1} pharmacophores were retrieved')
    # else:
    #     for i, pdb_string in enumerate(read_pdb_models(pdb_input)):
    #         p = get_pharmacophore(pdb_string, args.lig_id)
    #         p.save_to_xyz(os.path.splitext(pdb_input)[0] + f'{i:05d}.xyz')
    #         if args.verbose:
    #             sys.stderr.write(f'\r{i + 1} pharmacophores were retrieved')
    # if args.verbose:
    #     sys.stderr.write('\n')


if __name__ == '__main__':
    entry_point()
