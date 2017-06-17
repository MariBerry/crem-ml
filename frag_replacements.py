import sys
import os
import argparse

from rdkit import Chem
import sqlite3

from mutate import __frag_replace
from mutate import smiles_to_smarts



def read_worst_and_ids(input_worst, input_ids):
    """
    Prepare list of fragments with their ids
    :param input_worst: file name of selected fragments
    :param input_ids: file name with fragment ids
    :return: list of fragments with ids, e.g. [[compound_id, fragment_id, core, env, (fragment_ids)], [...], ...]
    """

    list_of_fragments = []

    # prepare list of the worst fragments
    with open(input_worst, 'r') as f_worst:
        for line in f_worst:
            line = line.split('\t')
            frag_core = line[2].split('|')[0]
            frag_env = line[2].split('|')[1]
            list_of_fragments.append(line[:-2])
            list_of_fragments[-1].append(frag_core)
            list_of_fragments[-1].append(frag_env)

    # prepare list of ids and connect them with fragments
    with open(input_ids, 'r') as f_ids:
        for line in f_ids:
            line = line.strip().split('\t')
            # join ids with fragments
            for i, fragment in enumerate(list_of_fragments):
                if line[1] == fragment[2] + '|' + fragment[3] and line[0] == fragment[0]:                  # compare fragment_name|fragment_context
                    list_of_fragments[i].append(tuple(map(int, line[2:])))       # append list of fragments ids
                    break
    return list_of_fragments


def make_replacements(input_sdf, input_worst, input_ids, db_cur, radius=3, min_size=1, max_size=10, min_rel_size=0,
                      max_rel_size=1, min_inc=-2, max_inc=+2):

    products = {}

    compounds = Chem.SDMolSupplier(input_sdf, removeHs=False, sanitize=False)
    list_of_fragments = read_worst_and_ids(input_worst, input_ids)

    for mol in compounds:
        mol.UpdatePropertyCache()
        # mol = Chem.AddHs(mol)
        mol_hac = mol.GetNumHeavyAtoms()
        mol_id = str(mol.GetProp('ID'))

        for fragment in list_of_fragments:
            if fragment[0] == mol_id:       # if we have same fragment from coresponging mol
                num_heavy_atoms = Chem.MolFromSmiles(fragment[2]).GetNumHeavyAtoms()
                hac_ratio = num_heavy_atoms / mol_hac

                if (min_size <= num_heavy_atoms <= max_size) or (min_rel_size <= hac_ratio <= max_rel_size):

                    frag_sma = smiles_to_smarts(fragment[2])

                    min_atoms = num_heavy_atoms + min_inc
                    max_atoms = num_heavy_atoms + max_inc

                    if radius == 3:
                        db_cur.execute("""SELECT core_smi, core_sma
                                          FROM radius3
                                          WHERE env IN (SELECT env FROM radius3 WHERE env = ?)
                                                AND
                                                core_num_atoms BETWEEN ? AND ?""", (fragment[3], min_atoms, max_atoms))
                    elif radius == 2:
                        db_cur.execute("""SELECT core_smi, core_sma
                                          FROM radius2
                                          WHERE env IN (SELECT env FROM radius2 WHERE env = ?)
                                                AND
                                                core_num_atoms BETWEEN ? AND ?""", (fragment[3], min_atoms, max_atoms))
                    rep = db_cur.fetchall()
                    id_mol = 0
                    for core_smi, core_sma in rep:
                        if core_smi != fragment[2]:
                            frag_replace_output, id_mol = __frag_replace(mol, frag_sma, core_sma, id_mol, fragment[-1])
                            for new_mol in frag_replace_output:
                                smi = Chem.MolToSmiles(new_mol, isomericSmiles=True)
                                if smi not in products:
                                    products[smi] = new_mol
    return list(products.values())




# input_sdf = 'output_process_predictions.sdf'
# input_worst = 'worst_fragments.txt'
# input_ids = 'fragment_ids.txt'
#
# output_product = 'new_compounds.sdf'
#
# conn = sqlite3.connect('../../replacement_chembl_cuts4_H.db')


def main_params(input_sdf, input_worst, input_ids, input_connection_db, output_product_file):
    # prepare database
    conn = sqlite3.connect(input_connection_db)
    db_cur = conn.cursor()

    products = make_replacements(input_sdf, input_worst, input_ids, db_cur)

    w = Chem.SDWriter(output_product_file)
    for m in products: w.write(m)

def main():
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

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "in_sdf": input_sdf = v
        if o == "in_worst": input_worst = v
        if o == "in_ids": input_ids = v
        if o == "in_con": input_connection_db = v
        if o == "out_compounds": output_product_file = v

    main_params(input_sdf, input_worst, input_ids, input_connection_db, output_product_file)


if __name__ == '__main__':
    main()
