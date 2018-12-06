import sys
import argparse
from collections import OrderedDict

import numpy as np
from rdkit import Chem


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


def make_replacements(input_sdf, input_worst, input_ids, path_to_db):
    new_products = []
    id_mol = 0

    compounds = Chem.SDMolSupplier(input_sdf, removeHs=False, sanitize=False)
    list_of_fragments = read_worst_and_ids(input_worst, input_ids)

    for a in list_of_fragments:
        print(a)

    for mol in compounds:
        mol_id = str(mol.GetProp('ID'))
        for frag in list_of_fragments:
            bad_mol_name, bad_frag_id = frag[0], list(frag[4])
            if mol_id == bad_mol_name:
                protected_ids = np.delete(np.arange(mol.GetNumAtoms()), bad_frag_id)
                from crem.crem import mutate_mol

                out = mutate_mol(
                    mol,
                    path_to_db,
                    radius=3,
                    min_size=0,
                    max_size=10,
                    min_inc=-2,
                    max_inc=2,
                    min_freq=0,
                    return_rxn=False,
                    ncores=1,
                    protected_ids=list(protected_ids)
                )
                print(out)
                for o in out:
                    print(o)
                sys.exit()
    #     # mol.UpdatePropertyCache()
    #     Chem.SanitizeMol(mol)
    #     mol_hac = mol.GetNumHeavyAtoms()
    #     mol_id = str(mol.GetProp('ID'))
    #     mol_name = str(mol.GetProp('_Name'))
    #
    #     for fragment in list_of_fragments:
    #         if fragment[0] == mol_id:       # if we have same fragment from coresponging mol
    #
    #             d = {}
    #
    #             if min_size == 0:
    #                 h = get_Hs(mol)
    #                 adj_ids = set()
    #                 for atom_id in fragment[-1]:
    #                     for nei in mol.GetAtomWithIdx(atom_id).GetNeighbors():
    #                         adj_ids.add(nei.GetIdx())
    #                 inter = adj_ids.intersection(h)
    #                 for i in inter:
    #                     d.update(replace(mol, '[H][*:1]', h[i], (i,), db_cur, min_atoms=1, max_atoms=3))
    #
    #             num_heavy_atoms = Chem.MolFromSmiles(fragment[2]).GetNumHeavyAtoms()
    #             hac_ratio = num_heavy_atoms / mol_hac
    #             if ((min_size <= num_heavy_atoms <= max_size) and (min_rel_size <= hac_ratio <= max_rel_size)) or (replace_cycles and cycle_pattern.search(core)):
    #                 min_atoms = num_heavy_atoms + min_inc
    #                 max_atoms = num_heavy_atoms + max_inc
    #
    #                 d.update(replace(mol, fragment[2], fragment[3], fragment[4], db_cur, min_atoms, max_atoms))
    #
    #             for new_mol in d.values():
    #                 new_mol.SetProp('_Name', 'ID_{}'.format(id_mol))
    #                 new_mol.SetProp('ID', 'ID_{}'.format(id_mol))
    #                 new_mol.SetProp('parent_name', mol_name)
    #                 new_products.append(new_mol)
    #                 id_mol += 1

    sys.exit()

    return new_products


def main(input_sdf, input_worst, input_ids, path_to_db, output_product_file):

    print('Replacing fragments ...')

    products = make_replacements(input_sdf, input_worst, input_ids, path_to_db)

    w = Chem.SDWriter(output_product_file)
    for m in products: w.write(m)
    w.close()

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

    args = vars(parser.parse_args())

    main(args['input_sdf'], args['input_worst'], args['input_ids'],
         args['input_connection_db'], args['output_product_file'])
