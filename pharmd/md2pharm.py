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
from CGRtools.files import XYZrw
from CGRtools.files.XYZrw import get_possible_bonds
from itertools import combinations, chain
from math import sqrt
from pickle import loads


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
# smarts = load_smarts(load_metal_chelators=True)
smarts = load_smarts()


queries = b'\x80\x03]q\x00(cCGRtools.containers.query\nQueryContainer\nq\x01)\x81q\x02}q\x03(X\x0c\x00\x00\x00atoms_stereoq\x04}q\x05X\t\x00\x00\x00neighborsq\x06}q\x07(KLK\x03\x85q\x08KMK\x01\x85q\tKNK\x01\x85q\nuX\x0e\x00\x00\x00hybridizationsq\x0b}q\x0c(KLK\x02\x85q\rKMK\x02\x85q\x0eKNK\x01\x85q\x0fuX\x05\x00\x00\x00atomsq\x10}q\x11(KLcCGRtools.periodictable.groupXIV\nQueryC\nq\x12)\x81q\x13}q\x14X\x07\x00\x00\x00isotopeq\x15NsbKMcCGRtools.periodictable.groupXVI\nQueryO\nq\x16)\x81q\x17}q\x18h\x15NsbKNh\x16)\x81q\x19}q\x1ah\x15NsbuX\x05\x00\x00\x00bondsq\x1b}q\x1c(KL}q\x1d(KMcCGRtools.containers.bonds\nBond\nq\x1e)\x81q\x1fN}q X\x0c\x00\x00\x00_Bond__orderq!K\x02s\x86q"bKNh\x1e)\x81q#N}q$h!K\x01s\x86q%buKM}q&KLh\x1fsKN}q\'KLh#suX\x04\x00\x00\x00metaq(}q)X\x05\x00\x00\x00planeq*}q+(KLcnumpy.core.multiarray\nscalar\nq,cnumpy\ndtype\nq-X\x02\x00\x00\x00f8q.K\x00K\x01\x87q/Rq0(K\x03X\x01\x00\x00\x00<q1NNNJ\xff\xff\xff\xffJ\xff\xff\xff\xffK\x00tq2bC\x08P3\x90\r\xd7\xeb\xea?q3\x86q4Rq5h,h0C\x08\xec3\x06\xa9C\xe6\xe2\xbfq6\x86q7Rq8\x86q9KMh,h0C\x08*\xc4\xab\x8f\xe6y\xf3?q:\x86q;Rq<h,h0C\x08\'\xf43\xcd\xb4*\xf5\xbfq=\x86q>Rq?\x86q@KNh,h0C\x08\x00\x00\x00\x00\x00\x00\x00\x00qA\x86qBRqCh,h0C\x08\xcad\x04z\x1a-\xe2\xbfqD\x86qERqF\x86qGuX\x0e\x00\x00\x00parsed_mappingqH}qIX\x07\x00\x00\x00chargesqJ}qK(KLK\x00KMK\x00KNJ\xff\xff\xff\xffuX\x08\x00\x00\x00radicalsqL}qM(KL\x89KM\x89KN\x89uX\x04\x00\x00\x00nameqNX\x00\x00\x00\x00qOubh\x01)\x81qP}qQ(h\x04}qRh\x06}qS(K\x01K\x01\x85qTK\x02K\x01\x85qUK\x03K\x01\x85qVK\x04K\x01\x85qWuh\x0b}qX(K\x01K\x01\x85qYK\x02K\x01\x85qZK\x03K\x01\x85q[K\x04K\x01\x85q\\uh\x10}q](K\x01cCGRtools.periodictable.groupXV\nQueryN\nq^)\x81q_}q`h\x15NsbK\x02cCGRtools.periodictable.groupI\nQueryH\nqa)\x81qb}qch\x15NsbK\x03ha)\x81qd}qeh\x15NsbK\x04ha)\x81qf}qgh\x15Nsbuh\x1b}qh(K\x01}qi(K\x02h\x1e)\x81qjN}qkh!K\x01s\x86qlbK\x03h\x1e)\x81qmN}qnh!K\x01s\x86qobK\x04h\x1e)\x81qpN}qqh!K\x01s\x86qrbuK\x02}qsK\x01hjsK\x03}qtK\x01hmsK\x04}quK\x01hpsuh(}qvh*}qw(K\x01h,h0C\x08f\xaaa\xf5+h\t@qx\x86qyRqzh,h0C\x08Z\x7f\xd9 \x18(\xeb\xbfq{\x86q|Rq}\x86q~K\x02h,h0C\x08\x9d\xc0\xca\'\xd7f\t@q\x7f\x86q\x80Rq\x81h,h0C\x08\xce\xdd&\xdb-\xa3\xfb\xbfq\x82\x86q\x83Rq\x84\x86q\x85K\x03h,h0C\x08\xf2\xdc\xcc\xb7d\xc2\x0c@q\x86\x86q\x87Rq\x88h,h0C\x08\xa0dcA\xae\xa1\xbd\xbfq\x89\x86q\x8aRq\x8b\x86q\x8cK\x04h,h0C\x08g\xd3\xbe\xfc&\x1c\x10@q\x8d\x86q\x8eRq\x8fh,h0C\x08\xe3\xf7\x0f\xadKD\xf1\xbfq\x90\x86q\x91Rq\x92\x86q\x93uhH}q\x94hJ}q\x95(K\x01K\x01K\x02K\x00K\x03K\x00K\x04K\x00uhL}q\x96(K\x01\x89K\x02\x89K\x03\x89K\x04\x89uhNhOubh\x01)\x81q\x97}q\x98(h\x04}q\x99h\x06}q\x9aMO\x12K\x00\x85q\x9bsh\x0b}q\x9cMO\x12K\x01\x85q\x9dsh\x10}q\x9eMO\x12cCGRtools.periodictable.groupXVII\nQueryCl\nq\x9f)\x81q\xa0}q\xa1h\x15Nsbsh\x1b}q\xa2MO\x12}q\xa3sh(}q\xa4h*}q\xa5MO\x12G\x00\x00\x00\x00\x00\x00\x00\x00G\x00\x00\x00\x00\x00\x00\x00\x00\x86q\xa6shH}q\xa7hJ}q\xa8MO\x12J\xff\xff\xff\xffshL}q\xa9MO\x12\x89shNhOubh\x01)\x81q\xaa}q\xab(h\x04}q\xach\x06}q\xad(M\xd5\x02K\x02\x85q\xaeM\xd6\x02K\x01\x85q\xafM\xd7\x02K\x03\x85q\xb0M\xd8\x02K\x01\x85q\xb1M\xd9\x02K\x01\x85q\xb2M\xda\x02K\x01\x85q\xb3M\xdb\x02K\x01\x85q\xb4M\xdc\x02K\x01\x85q\xb5M\xdd\x02K\x01\x85q\xb6uh\x0b}q\xb7(M\xd5\x02K\x01\x85q\xb8M\xd6\x02K\x01\x85q\xb9M\xd7\x02K\x02\x85q\xbaM\xd8\x02K\x01\x85q\xbbM\xd9\x02K\x01\x85q\xbcM\xda\x02K\x01\x85q\xbdM\xdb\x02K\x02\x85q\xbeM\xdc\x02K\x01\x85q\xbfM\xdd\x02K\x01\x85q\xc0uh\x10}q\xc1(M\xd5\x02h^)\x81q\xc2}q\xc3h\x15NsbM\xd6\x02ha)\x81q\xc4}q\xc5h\x15NsbM\xd7\x02h\x12)\x81q\xc6}q\xc7h\x15NsbM\xd8\x02h^)\x81q\xc8}q\xc9h\x15NsbM\xd9\x02ha)\x81q\xca}q\xcbh\x15NsbM\xda\x02ha)\x81q\xcc}q\xcdh\x15NsbM\xdb\x02h^)\x81q\xce}q\xcfh\x15NsbM\xdc\x02ha)\x81q\xd0}q\xd1h\x15NsbM\xdd\x02ha)\x81q\xd2}q\xd3h\x15Nsbuh\x1b}q\xd4(M\xd5\x02}q\xd5(M\xd6\x02h\x1e)\x81q\xd6N}q\xd7h!K\x01s\x86q\xd8bM\xd7\x02h\x1e)\x81q\xd9N}q\xdah!K\x01s\x86q\xdbbuM\xd6\x02}q\xdcM\xd5\x02h\xd6sM\xd7\x02}q\xdd(M\xd5\x02h\xd9M\xd8\x02h\x1e)\x81q\xdeN}q\xdfh!K\x01s\x86q\xe0bM\xdb\x02h\x1e)\x81q\xe1N}q\xe2h!K\x02s\x86q\xe3buM\xd8\x02}q\xe4(M\xd7\x02h\xdeM\xd9\x02h\x1e)\x81q\xe5N}q\xe6h!K\x01s\x86q\xe7bM\xda\x02h\x1e)\x81q\xe8N}q\xe9h!K\x01s\x86q\xeabuM\xd9\x02}q\xebM\xd8\x02h\xe5sM\xda\x02}q\xecM\xd8\x02h\xe8sM\xdb\x02}q\xed(M\xd7\x02h\xe1M\xdc\x02h\x1e)\x81q\xeeN}q\xefh!K\x01s\x86q\xf0bM\xdd\x02h\x1e)\x81q\xf1N}q\xf2h!K\x01s\x86q\xf3buM\xdc\x02}q\xf4M\xdb\x02h\xeesM\xdd\x02}q\xf5M\xdb\x02h\xf1suh(}q\xf6h*}q\xf7(M\xd5\x02h,h0C\x08\xb9\x87\xde9\x12U\xfd?q\xf8\x86q\xf9Rq\xfah,h0C\x08\xb8y\xf0\xc8%\x84\xcf?q\xfb\x86q\xfcRq\xfd\x86q\xfeM\xd6\x02h,h0C\x08\x18a\xf6\xd5\x87I\x05@q\xff\x86r\x00\x01\x00\x00Rr\x01\x01\x00\x00h,h0C\x08\x80O[\xf9)\xe1\x82\xbfr\x02\x01\x00\x00\x86r\x03\x01\x00\x00Rr\x04\x01\x00\x00\x86r\x05\x01\x00\x00M\xd7\x02h,h0C\x08\x18\x94\x8d\xaae\x95\xf4?r\x06\x01\x00\x00\x86r\x07\x01\x00\x00Rr\x08\x01\x00\x00h,h0C\x08ha\x88Wa\x83\xd6\xbfr\t\x01\x00\x00\x86r\n\x01\x00\x00Rr\x0b\x01\x00\x00\x86r\x0c\x01\x00\x00M\xd8\x02h,h0C\x08\xbc\x15v\xf5\xb2K\xde?r\r\x01\x00\x00\x86r\x0e\x01\x00\x00Rr\x0f\x01\x00\x00h,h0C\x08\xe0\x1df4y\x98\xd4\xbfr\x10\x01\x00\x00\x86r\x11\x01\x00\x00Rr\x12\x01\x00\x00\x86r\x13\x01\x00\x00M\xd9\x02h,h0C\x08`\xb55\x0c8\x94\xac?r\x14\x01\x00\x00\x86r\x15\x01\x00\x00Rr\x16\x01\x00\x00h,h0C\x08\xcc\xca\x0b\xd5\xb1-\xf1\xbfr\x17\x01\x00\x00\x86r\x18\x01\x00\x00Rr\x19\x01\x00\x00\x86r\x1a\x01\x00\x00M\xda\x02h,h0C\x08\x00\x00\x00\x00\x00\x00\x00\x00r\x1b\x01\x00\x00\x86r\x1c\x01\x00\x00Rr\x1d\x01\x00\x00h,h0C\x08\xc8\x14\x0cA\xb9\x82\xd6?r\x1e\x01\x00\x00\x86r\x1f\x01\x00\x00Rr \x01\x00\x00\x86r!\x01\x00\x00M\xdb\x02h,h0C\x08\xaa\xda4\xc7u\xe3\xf9?r"\x01\x00\x00\x86r#\x01\x00\x00Rr$\x01\x00\x00h,h0C\x08\xf0\xc3\xf1\xaa\xe4N\xf2\xbfr%\x01\x00\x00\x86r&\x01\x00\x00Rr\'\x01\x00\x00\x86r(\x01\x00\x00M\xdc\x02h,h0C\x08(\xd2.\x03M\x81\x03@r)\x01\x00\x00\x86r*\x01\x00\x00Rr+\x01\x00\x00h,h0C\x08\xf66\x8a\xa2\xd7\xa5\xf3\xbfr,\x01\x00\x00\x86r-\x01\x00\x00Rr.\x01\x00\x00\x86r/\x01\x00\x00M\xdd\x02h,h0C\x08\xa07\x87\xde5\xa1\xf1?r0\x01\x00\x00\x86r1\x01\x00\x00Rr2\x01\x00\x00h,h0C\x08\xe6/{C\xbb\x98\xfc\xbfr3\x01\x00\x00\x86r4\x01\x00\x00Rr5\x01\x00\x00\x86r6\x01\x00\x00uhH}r7\x01\x00\x00hJ}r8\x01\x00\x00(M\xd5\x02K\x00M\xd6\x02K\x00M\xd7\x02K\x00M\xd8\x02K\x00M\xd9\x02K\x00M\xda\x02K\x00M\xdb\x02K\x01M\xdc\x02K\x00M\xdd\x02K\x00uhL}r9\x01\x00\x00(M\xd5\x02\x89M\xd6\x02\x89M\xd7\x02\x89M\xd8\x02\x89M\xd9\x02\x89M\xda\x02\x89M\xdb\x02\x89M\xdc\x02\x89M\xdd\x02\x89uhNhOube.'
queries = loads(queries)


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
