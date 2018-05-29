import sys
import os
import argparse
import re
from collections import OrderedDict

from rdkit import Chem
from rdkit.Chem import rdMMPA
import sqlite3

from mutate import __frag_replace
from mutate import smiles_to_smarts, get_canon_context_core

from find_frags_auto_rdkit import replace_no2

cycle_pattern = re.compile("(?<!:)[1-9]+")


def read_worst_and_ids(input_worst, input_ids):
    """
    Prepare list of fragments with their ids
    :param input_worst: file name of selected fragments
    :param input_ids: file name with fragment ids
    :return: list of fragments with ids, e.g. [[compound_id, fragment_id, core, env, (fragment_ids)], [...], ...]
    """

    list_of_fragments = []
    d = OrderedDict()

    # prepare list of the worst fragments
    with open(input_worst, 'r') as f_worst:
        f_worst.readline()
        for line in f_worst:
            line = line.split('\t')
            frag_core = line[2].split('|')[0]
            frag_env = line[2].split('|')[1]
            d[int(line[1])] = line[:2] + [frag_core, frag_env]

    # prepare list of ids and connect them with fragments
    with open(input_ids, 'r') as f_ids:
        for i, line in enumerate(f_ids):
            line = line.strip().split('\t')
            if i in d:
                d[i].append(tuple(j-1 for j in map(int, line[2:])))
    return list(d.values())


def get_Hs(mol, radius=3, keep_stereo=False):

    def get_atom_prop(molecule, prop="Index"):
        res = []
        for a in molecule.GetAtoms():
            try:
                res.append(a.GetIntProp(prop))
            except KeyError:
                continue
        return tuple(sorted(res))

    output = {}
    for atom in mol.GetAtoms():
        atom.SetIntProp("Index", atom.GetIdx())
    frags = rdMMPA.FragmentMol(mol, pattern="[#1]!@!=!#[!#1]", maxCuts=1, resultsAsMols=True, maxCutBonds=100)
    for _, chains in frags:
        components = list(Chem.GetMolFrags(chains, asMols=True))
        ids_0 = get_atom_prop(components[0])
        ids_1 = get_atom_prop(components[1])
        if Chem.MolToSmiles(components[0]) != '[H][*:1]':  # context cannot be H
            env, frag = get_canon_context_core(components[0], components[1], radius, keep_stereo)
            output[ids_1[0]] = env
        if Chem.MolToSmiles(components[1]) != '[H][*:1]':  # context cannot be H
            env, frag = get_canon_context_core(components[1], components[0], radius, keep_stereo)
            output[ids_0[0]] = env
    return output


def replace(in_mol, frag_core, frag_env, frag_ids, db_cur, min_atoms=0, max_atoms=0, radius=3):
    new_mols = {}

    frag_sma_core = smiles_to_smarts(frag_core)

    if radius == 3:
        db_cur.execute("""SELECT core_smi, core_sma
                          FROM radius3
                          WHERE env IN (SELECT env FROM radius3 WHERE env = ?)
                                AND
                                core_num_atoms BETWEEN ? AND ?""", (frag_env, min_atoms, max_atoms))
    elif radius == 2:
        db_cur.execute("""SELECT core_smi, core_sma
                          FROM radius2
                          WHERE env IN (SELECT env FROM radius2 WHERE env = ?)
                                AND
                                core_num_atoms BETWEEN ? AND ?""", (frag_env, min_atoms, max_atoms))
    rep = db_cur.fetchall()
    for core_smi, core_sma in rep:
        if core_smi != frag_core:
            frag_replace_output = __frag_replace(in_mol, frag_sma_core, core_sma, frag_ids)

            for new_mol in frag_replace_output:
                smi = Chem.MolToSmiles(new_mol, isomericSmiles=True)
                if smi not in new_mols:
                    new_mols[smi] = new_mol
    return new_mols



def make_replacements(input_sdf, input_worst, input_ids, db_cur, n_compounds, radius=3,
                      min_size=0, max_size=7, min_rel_size=0, max_rel_size=0.5,
                      min_inc=-2, max_inc=+2, replace_cycles=False):

    new_products = []
    id_mol = n_compounds

    compounds = Chem.SDMolSupplier(input_sdf, removeHs=False, sanitize=False)
    list_of_fragments = read_worst_and_ids(input_worst, input_ids)

    for mol in compounds:
        mol = replace_no2(mol)
        # mol.UpdatePropertyCache()
        Chem.SanitizeMol(mol)
        mol_hac = mol.GetNumHeavyAtoms()
        mol_id = str(mol.GetProp('ID'))
        mol_name = str(mol.GetProp('_Name'))

        for fragment in list_of_fragments:
            if fragment[0] == mol_id:       # if we have same fragment from coresponging mol

                d = {}

                if min_size == 0:
                    h = get_Hs(mol)
                    adj_ids = set()
                    for atom_id in fragment[-1]:
                        for nei in mol.GetAtomWithIdx(atom_id).GetNeighbors():
                            adj_ids.add(nei.GetIdx())
                    inter = adj_ids.intersection(h)
                    for i in inter:
                        d.update(replace(mol, '[H][*:1]', h[i], (i,), db_cur, min_atoms=1, max_atoms=3))

                num_heavy_atoms = Chem.MolFromSmiles(fragment[2]).GetNumHeavyAtoms()
                hac_ratio = num_heavy_atoms / mol_hac
                if ((min_size <= num_heavy_atoms <= max_size) and (min_rel_size <= hac_ratio <= max_rel_size)) or (replace_cycles and cycle_pattern.search(core)):
                    min_atoms = num_heavy_atoms + min_inc
                    max_atoms = num_heavy_atoms + max_inc

                    d.update(replace(mol, fragment[2], fragment[3], fragment[4], db_cur, min_atoms, max_atoms))

                for new_mol in d.values():
                    new_mol.SetProp('_Name', 'ID{}'.format(id_mol))
                    new_mol.SetProp('ID', 'ID{}'.format(id_mol))
                    new_mol.SetProp('parent_name', mol_name)
                    new_products.append(new_mol)
                    id_mol += 1

    return new_products, id_mol


def main(input_sdf, input_worst, input_ids, input_connection_db, output_product_file, n_compounds):

    print('Replacing fragments ...')

    # prepare database
    conn = sqlite3.connect(input_connection_db)
    db_cur = conn.cursor()

    products, last_id_mol = make_replacements(input_sdf, input_worst, input_ids, db_cur, n_compounds)

    w = Chem.SDWriter(output_product_file)
    for m in products: w.write(m)
    w.close()
    conn.close()
    return last_id_mol


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=
                                      'Create new file with replaced fragments.')
    parser.add_argument('-is', '--in_sdf', metavar='output_process_predictions.sdf', required=True,
                         help='path to file which contains selected compounds from pareto/desirability/...')
    parser.add_argument('-iw', '--in_worst', metavar='worst_fragments.txt', required=True,
                         help='path to the file where you store the worst fragments')
    parser.add_argument('-id', '--in_ids', metavar='fragments_ids.txt', required=True,
                        help='path to the file where you store ids of fragments')
    parser.add_argument('-ic', '--in_con', metavar='database.db', required=True,
                        help='path to the database with fragment replacements')
    parser.add_argument('-oc', '--out_compounds', metavar='new_compounds.sdf', required=True,
                        help='file name where you want to store new compounds')
    parser.add_argument('-n', '--n_compounds', action='store', type=int,
                        help='specifies number of used compounds, for new IDs')

    args = vars(parser.parse_args())

    main(args['input_sdf'], args['input_worst'], args['input_ids'],
         args['input_connection_db'], args['output_product_file'],
         args['n_compounds'])
