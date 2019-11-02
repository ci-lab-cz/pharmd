#!/usr/bin/env python3

import argparse
import itertools
import os
import numpy
import pybel
import mdtraj as md
from collections import defaultdict
from plip.modules.preparation import PDBComplex
#from rdkit import Chem


def create_parser():
    parser = argparse.ArgumentParser(description='Extract pharmacophore models from an MD trajectory of a '
                                                 'protein-ligand complex.')
    parser.add_argument('-i', '--input', metavar='input.xtc', required=True, type=str,
                        help='Input file with MD trajectory. Formats are the same as MDTraj supports.')
    parser.add_argument('-t', '--topology', metavar='input.pdb', required=True, type=str,
                        help='PDB file with topology')
    parser.add_argument('-f', '--first', metavar='INTEGER', required=False, type=int,
                        help='Staring frame number.')
    parser.add_argument('-l', '--last', metavar='INTEGER', required=False, type=int,
                        help='Last frame number.')
    parser.add_argument('-s', '--stride', metavar='INTEGER', required=False, type=int,
                        help='Step using to extract frames.')
    parser.add_argument('-g', '--lig_id', metavar='STRING', required=True, type=str,
                        help='three-letter ligand ID')
    parser.add_argument('-o', '--output', metavar='output.pdb', required=True,
                        help='output PDB file with all extracted frames. Solvent will be omitted.')
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='print progress to STDERR.')
    return parser


def read_smarts(smarts_file):
    with open(smarts_file, 'r') as f:
        dict_smarts = defaultdict(list)  # contains all smarts
        sma = []
        new_smarts_features = list(f)
        for line in new_smarts_features:
            if line[0] != '#' and line != '\n':
                sma.append(line.strip().split())
        for feature, label in sma:
            if label == 'D':
                dict_smarts['hbd'].append(feature)
            elif label == 'A':
                dict_smarts['hba'].append(feature)
            elif label == 'a':
                dict_smarts['pi'].append(feature)
            elif label == 'H':
                dict_smarts['h'].append(feature)
            elif label == 'N':
                dict_smarts['n'].append(feature)
            elif label == 'P':
                dict_smarts['p'].append(feature)
    return dict_smarts


def load_complex_from_string(pdb_string, lig_id):
    complex = PDBComplex()  # plip molecule
    complex.load_pdb(pdb_string, as_string=True)
    lig_pdb_string = []
    for line in pdb_string.split('\n'):
        if line[17:20] == lig_id:
            lig_pdb_string.append(line)
    lig_pdb_string = '\n'.join(lig_pdb_string)
    lig = pybel.readstring('pdb', lig_pdb_string)  # openbabel molecule
    # lig = Chem.MolFromPDBBlock(lig_pdb_string, removeHs=False)  # rdkit molecule
    return complex, lig


def plip_analysis(my_mol, lig_index):
    true_lig = my_mol.ligands[lig_index]
    hetid = true_lig.hetid
    chain = true_lig.chain
    position = true_lig.position
    my_bsid = f'{hetid}:{chain}:{position}'
    my_mol.analyze()
    my_interactions = my_mol.interaction_sets[my_bsid]  # find interactions only for one ligand
    return my_interactions


def find_smarts_hb(dict_smarts, smarts_hb, lig_mol):
    # find atoms which can form a hydrogen bond and matches smarts
    hb_ph = set()
    for smart in dict_smarts[smarts_hb]:
        sma = pybel.Smarts(smart)
        all_coords = sma.findall(lig_mol)
        for coords in all_coords:
            for coor in coords:
                hb_mol = lig_mol.atoms[coor - 1]
                hb_ph.add(hb_mol.coords)
    return hb_ph


def find_smarts(dict_smarts, smarts, lig_mol):
    # find groups of atoms which matches smarts
    coords_interaction = []
    for smart in dict_smarts[smarts]:
        sma = pybel.Smarts(smart)
        all_coords = sma.findall(lig_mol)
        for coords in all_coords:
            c_ph = []
            for coor in coords:
                c_mol = lig_mol.atoms[coor - 1]
                c_ph.append(c_mol.coords)
            coords_interaction.append(c_ph)
    return coords_interaction


def filter_features(atom_ids):
    # function exclude tuples if they are a complete subset of a bigger tuples to avoid superfluous features
    # for example the hydrophobic features has ids (1,2) and (1,2,3,4,5), only the second will be kept
    # another example the hydrophobic features has ids (1,2) and (2,3,4,5), both will be kept
    # atom ids is a list of tuples with atom ids [(1,), (1,2), (2,3,4)]
    # OUTPUT: tuple of tuples ((1,2), (2,3,4))
    tmp = sorted(tuple(set(atom_ids)), key=len, reverse=True)
    res = [tmp[0]]
    for item in tmp[1:]:
        item = set(item)
        if any(item.issubset(bigger) for bigger in res):
            continue
        else:
            res.append(tuple(item))
    return tuple(res)


def center_interactions(interaction):
    # find geometry center of the group
    x = []
    y = []
    z = []
    for coords in interaction:
        x.append(coords[0])
        y.append(coords[1])
        z.append(coords[2])
    l = len(interaction)
    center_coords = (round(sum(x)/l, 3), round(sum(y)/l, 3), round(sum(z)/l, 3))
    return center_coords


def center_interactions_el(all_coords):
    # find the center of group participating in electrostatic interactions
    centers = []
    for interact in all_coords:
        center_coords = center_interactions(interact)
        centers.append(center_coords)
    return centers


def dist_euclidean(check):
    # returns the list of the distances between ligand and protein
    # and the dict of the ligand coordinates with the distance to the protein
    distance = {}
    d = []
    for center, prot in check.items():
        xl = center[0]
        yl = center[1]
        zl = center[2]
        xp = prot[0]
        yp = prot[1]
        zp = prot[2]
        dist = numpy.sqrt((xl-xp)**2+(yl-yp)**2+(zl-zp)**2)
        d.append(dist)
        distance[center] = dist
    return d, distance


def idx_coords_ligand(my_mol, coords_interactions):
    # Assignment the indexes for atoms coordinates
    # returns the list of the dict where key - index of ligand atom in pdb, value - coords of ligand atom
    dict_atoms_index = {}  # keys - atoms coords, values - atoms index in pdb
    for atom in my_mol.atoms.values():
        dict_atoms_index[atom.coords] = atom.idx
    all_index_ligand = []
    for group in coords_interactions:
        index_ligand = {}
        for coords in group:
            index_ligand[dict_atoms_index[coords]] = coords
        all_index_ligand.append(index_ligand)
    return all_index_ligand


def find_right_distance(all_one_type, all_index_ligand):
    # returns the list of dicts with idex and coords of ligand atom
    # with distance between ligand atom and protein below the threshold
    dist = []
    for Type_atoms in all_one_type:
        for group_idxes in all_index_ligand:
            for Type_atom, idx in itertools.product(Type_atoms, group_idxes):
                if 0.5 < Type_atom.OBAtom.GetDistance(idx) < 4:
                    dist.append(group_idxes)
    return dist


def true_coords(dist):
    # find coordinates of ligands which satisfy right distance
    all_coords = []
    for interaction in dist:
        coords = []
        for crds in interaction.values():
            coords.append(crds)
        if coords not in all_coords:
            all_coords.append(coords)
    return all_coords


def find_positive(mol):
    # returns the list of list of nitrogen atoms in residue which can interact with ligand
    positive_reses = []
    all_mol_residues = mol.residues
    for res in all_mol_residues:
        if res.name == 'LYS' or res.name == 'ARG' or res.name == 'HIS':
            positive_reses.append(res)
    groups_atoms = []
    for res in positive_reses:
        groups_atoms.append(res.atoms)
    all_type_N = []
    for group in groups_atoms:
        type_N_atom = []
        for atom in group:
            if 'N' in atom.type and not 'Nam' in atom.type:
                type_N_atom.append(atom)
        all_type_N.append(type_N_atom)
    return all_type_N


def find_negative(mol):
    # returns the list of list of oxygen  atoms in residue which can interact with ligand
    name_res_neg = []
    negative_reses = []
    all_mol_residues = mol.residues
    for res in all_mol_residues:
        if res.name == 'GLU' or res.name == 'ASP':  # find GLU and ASP in protein
            negative_reses.append(res)
            name_res_neg.append(res.name)
    groups_atoms_neg = []  # add pybel.atoms in each residue GLU and ASP
    for res in negative_reses:
        groups_atoms_neg.append(res.atoms)
    all_type_O = []  # add atoms of GLU and ASP which can participate in electrostatic interations as positive center
    for group in groups_atoms_neg:
        type_O_atom = []
        for atom in group:
            if 'O.co2' in atom.type:
                type_O_atom.append(atom)
        all_type_O.append(type_O_atom)
    return all_type_O


def find_hydrogen_bonds(interactions, lig_mol):
    hb_ldon = set(hb_ldon.d.coords for hb_ldon in interactions.hbonds_ldon)
    hb_pdon = set(hb_pdon.a.coords for hb_pdon in interactions.hbonds_pdon)
    d_ph = find_smarts_hb(dict_smarts, 'hbd', lig_mol)
    all_hba_ph = find_smarts_hb(dict_smarts, 'hba', lig_mol)
    donor = hb_ldon.intersection(d_ph)
    acceptor = hb_pdon.intersection(all_hba_ph)
    return donor, acceptor


def find_aromatic(interactions, lig_mol):
    pistacking = [tuple(pistacking.ligandring.center) for pistacking in interactions.pistacking]
    true_pistacking = set()
    for interaction in pistacking:
        rounded = []
        for coords in interaction:
            rounded.append(round(coords, 3))
        true_pistacking.add(tuple(rounded))
    pi_ph_group = find_smarts(dict_smarts, 'pi', lig_mol)
    centers_aromatic = set()
    for interaction in pi_ph_group:
        center_coords = center_interactions(interaction)
        centers_aromatic.add(center_coords)
    aromatic = true_pistacking.intersection(centers_aromatic)
    return aromatic


def find_hydrophobic(interactions, lig_mol):
    hc = [hydrophobic_contacts for hydrophobic_contacts in interactions.hydrophobic_contacts]
    hydrophobic = {}  # dict ligand coords, which participates in hydrophobic interactions - protein coords
    for contact in hc:
        hydrophobic[contact.ligatom.coords] = contact.bsatom.coords
    hydrophobic_group = find_smarts(dict_smarts, 'h', lig_mol)
    h_g = []
    for i in hydrophobic_group:
        h_g.append(tuple(i))
    f_hydrophobic_group = filter_features(h_g)
    true_hydrophobic_group = {}  # key hydrophobic atom's coords in group - value protein coords
    for atom, protein in hydrophobic.items():  # matching with smarts
        for group in f_hydrophobic_group:
            if atom in group:
                true_hydrophobic_group[group] = protein
    hyd_centers = {}  # key coordinates of center of group - value coordinates of protein
    for group, prot in true_hydrophobic_group.items():
        hyd_centers[center_interactions(group)] = prot
    need_check = {}
    tmp = []
    for protein in hyd_centers.values():
        tmp.append(protein)
    for group, protein in hyd_centers.items():
        if tmp.count(protein) > 1:
            need_check[group] = protein
    d, distance = dist_euclidean(need_check)
    min_distance = []
    for center, dist in distance.items():
        if min(d) == dist:
            min_distance.append(center)
    hydrophobic = set(hyd_centers).difference(set(need_check))
    hydrophobic.update(set(min_distance))
    return hydrophobic


def find_electrostatic(dict_smarts, complex, lig_mol, lig):
    all_n_ph = find_smarts(dict_smarts, 'n', lig_mol)
    all_p_ph = find_smarts(dict_smarts, 'p', lig_mol)
    idx_p = idx_coords_ligand(complex, all_p_ph)
    idx_n = idx_coords_ligand(complex, all_n_ph)
    all_type_N = find_positive(lig)
    all_type_O = find_negative(lig)
    dist_negative = find_right_distance(all_type_N, idx_n)
    all_coords_negative = true_coords(dist_negative)
    dist_positive = find_right_distance(all_type_O, idx_p)
    all_coords_positive = true_coords(dist_positive)
    negative = center_interactions_el(all_coords_negative)
    positive = center_interactions_el(all_coords_positive)
    return positive, negative


def check_interactions(interactions, complex, lig, lig_index):
    checked_interactions = defaultdict(list)
    aromatic = find_aromatic(interactions, complex.ligands[lig_index].mol)
    if aromatic:
        checked_interactions['a'].append(aromatic)
    donor, acceptor = find_hydrogen_bonds(interactions, complex.ligands[lig_index].mol)
    if acceptor:
        checked_interactions['A'].append(acceptor)
    if donor:
        checked_interactions['D'].append(donor)
    hydrophobic = find_hydrophobic(interactions, complex.ligands[lig_index].mol)
    if hydrophobic:
        checked_interactions['H'].append(hydrophobic)
    positive, negative = find_electrostatic(dict_smarts, complex, complex.ligands[lig_index].mol, lig)
    if positive:
        checked_interactions['P'].append(positive)
    if negative:
        checked_interactions['N'].append(negative)
    return checked_interactions


def writeInteractions(fname, interactions):
    # writes xyz file with all interactions for one frame
    with open(fname, 'wt') as xyz:
        xyz.write('\n\n')
        for label, features in interactions.items():
            for feature in features:
                for coords in feature:
                    xyz.write('{0:s}{1:12.5f}{2:12.5f}{3:12.5f}\n'.format(label, coords[0], coords[1], coords[2]))


def read_pdb_models(fname):
    with open(fname) as f:
        lines = []
        for line in f:
            if not line.startswith("MODEL") and not line.startswith("ENDMDL"):
                lines.append(line)
            if line.startswith("ENDMDL"):
                yield ''.join(lines)
                lines = []


def get_lig_index(complex, lig_id):
    index = None
    for ligand in complex.ligands:
        if lig_id in ligand:
            index = complex.ligands.index(ligand)
    return index


def entry_point():

    global dict_smarts

    dict_smarts = read_smarts(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'smarts_features.txt'))

    parser = create_parser()
    args = parser.parse_args()
    args.lig_id = args.lig_id.upper()

    traj = md.load(args.input, top=args.topology)
    traj.remove_solvent(inplace=True)
    traj[args.first:args.last:args.stride].save_pdb(args.output)

    for i, pdb_string in enumerate(read_pdb_models(args.output)):
        complex, lig = load_complex_from_string(pdb_string, args.lig_id)
        lig_index = get_lig_index(complex, args.lig_id)
        interactions = plip_analysis(complex, lig_index)
        interactions = check_interactions(interactions, complex, lig, lig_index)
        writeInteractions(os.path.splitext(args.output)[0] + f'{i:05d}.xyz', interactions)


if __name__ == '__main__':
    entry_point()
