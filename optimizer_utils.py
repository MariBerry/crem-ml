import os
import sys
import shutil

from subprocess import call
import sqlite3 as lite
from rdkit import Chem

from typing import List

sys.path.insert(1, os.path.join(sys.path[0], 'spci'))
import calc_atomic_properties_chemaxon
import filter_descriptors
import predict

sys.path.insert(1, os.path.join(sys.path[0], 'spci/sirms'))
import sirms


def create_database(working_dir: str, parameter_to_optimize: List) -> str:
    """
    Define new database and save it to working_dir.

    :param working_dir: path to working directory, new db will be stored there
    :param parameter_to_optimize: list of all parameters
    :return: path_to_database
    """
    parameter_to_optimize = ["predicted_{}".format(parameter) for parameter in parameter_to_optimize]
    path_to_database = os.path.join(working_dir, 'output.db')

    table_str = "CREATE TABLE optimizer_table"\
                    "(id TEXT NOT NULL,"\
                    "smi TEXT NOT NULL,"\
                    "generation INTEGER NOT NULL,"\
                    "parent TEXT,"\
                    "transformation TEXT,"\
                    "fit INTEGER,"
    table_str += " REAL,".join(parameter_to_optimize) + " REAL, overall_prediction REAL)"

    if os.path.isfile(path_to_database):
        os.remove(path_to_database)

    con = lite.connect(path_to_database)
    with con:
        cursor = con.cursor()
        cursor.execute(table_str)
        cursor.execute("CREATE INDEX idx ON optimizer_table (id)")
        cursor.execute("CREATE INDEX smi_idx ON optimizer_table (smi)")
        cursor.execute("DELETE FROM optimizer_table")

    return path_to_database

def add_mols_into_db(num_of_compounds: int, input_sdf: str, database: str, gen: int) -> int:
    """
    Read input sdf file, convert all mols into smiles, check if they are in DB,
    and if not add them with all possible options, such as transformation rules,
    parents, number of generation and so on.

    :param num_of_compounds: actual number of compounds in output database
    :param input_sdf: path to input sdf file
    :param database: path to output database
    :param gen: actual generation of optimization
    :return: number of compounds in database
    """
    # get generator of mols in sdf file
    supplier = Chem.SDMolSupplier(input_sdf)

    if gen == 0:
        new_sdf_path = os.path.join(os.path.dirname(input_sdf), 'tmp.sdf')
        new_sdf = Chem.SDWriter(new_sdf_path)

    con = lite.connect(database)
    with con:
        cursor = con.cursor()

        cursor.execute("SELECT smi FROM optimizer_table")
        mols_in_db = [mol[0] for mol in cursor.fetchall()]

        insert = []

        for mol in supplier:
            smile = Chem.MolToSmiles(mol)

            # mol doesn't have parent and transformation prop if it is in zero gen
            if gen == 0:
                mol.SetProp("_Name", "ID{}".format(str(num_of_compounds)))
                mol.SetProp("ID", "ID{}".format(str(num_of_compounds)))
                mol.SetProp("parent", "None")
                mol.SetProp("transformation", "None")
                new_sdf.write(mol)

            if smile not in mols_in_db:
                insert.append(("ID" + str(num_of_compounds),
                               smile,
                               gen,
                               mol.GetProp('parent'),
                               mol.GetProp('transformation')))
                mols_in_db.append(smile)
                num_of_compounds += 1

        cursor.executemany("INSERT INTO optimizer_table (id, smi, generation, parent, transformation) VALUES (?, ?, ?, ?, ?)", insert)
        con.commit()

        if gen == 0:
            new_sdf.close()
            os.remove(input_sdf)
            os.rename(new_sdf_path, input_sdf)

    return num_of_compounds

def quote_str(s: str) -> str:
    """
    Quote string

    :param s: input string
    :return: quoted string
    """

    return "'%s'" % s

def standardize_sdf(input_sdf_file: str, std_rules_path: str,
                    chemaxon_path: str, copy_rules: bool=False) -> str:
    """
    Create file with standardized compounds

    :param input_sdf_file: path to sdf file with compounds
    :param std_rules_path: path to file with rules for standardization
    :param chemaxon_path: path to chemaxon bin
    :param copy_rules: if specified, copy rules to output directory
    :return: path to new sdf file
    """

    print('Standardization is in progress...')

    # copy xml-rules if specified
    if copy_rules:
        shutil.copyfile(
            std_rules_path, os.path.join(os.path.dirname(input_sdf_file), std_rules_path.split("/")[-1]))

    # run standardize
    std_sdf = os.path.join(os.path.dirname(input_sdf_file), 'input_dataset_std.sdf')
    run_params = [os.path.join(chemaxon_path, 'standardize'),
                  '-c',
                  quote_str(std_rules_path),  # path to rules
                  quote_str(input_sdf_file),  # path to input file
                  '-f',
                  'sdf',  # type of output file
                  '-o',
                  quote_str(std_sdf)]  # name of output file
    call(' '.join(run_params), shell=True)

    return std_sdf

def calculate_atomic_prop(input_sdf_file: str, chemaxon_path: str, properties: List) -> str:
    """
    Calculate atomic properties with Chemaxon, it creates file with labeled compounds

    :param input_sdf_file: path to sdf file with compounds
    :param chemaxon_path: path to chemaxon bin
    :param properties: list of properties, e.g. ['charge', 'refractivity', 'logp', ...]
    :return: path to new sdf file
    """
    print('Atomic properties calculation is in progress...')
    lbl_sdf = os.path.join(os.path.dirname(input_sdf_file), 'input_dataset_std_lbl.sdf')
    calc_atomic_properties_chemaxon.main_params(input_sdf_file,
                                                lbl_sdf,
                                                properties,
                                                None,
                                                os.path.join(chemaxon_path, 'cxcalc'))
    return lbl_sdf

def calculate_sirms_descriptors(input_sdf_file: str, setup_file: str,
                                properties: List, output_format: str,
                                n_cores: int, copy_setup: bool=True,
                                fragments_fname=None):
    """
    Create files with descriptors

    :param input_sdf_file: path to standardized and labeled sdf file
    :param setup_file: path to file with setup for calculation of sirms descriptors
    :param properties: list of properties, e.g. ['CHARGE', 'REFRACTIVITY', 'LOGP', ...]
    :param output_format: svm
    :param n_cores: number of cores for computing
    :param copy_setup: if specified, copy setup file to output directory
    :param fragments_fname: if specified, use fragments ids
    """

    print("Descriptors calculation started. Please wait it can take some time")

    # copy setup file for sirms into generation dir
    if copy_setup:
        shutil.copyfile(setup_file,
                        os.path.join(os.path.dirname(input_sdf_file), os.path.basename(setup_file))
                        )

    # define output files
    if fragments_fname is not None:
        x_fname = os.path.join(os.path.dirname(input_sdf_file), 'new_x.txt')
    else:
        x_fname = os.path.join(os.path.dirname(input_sdf_file), 'x.txt')

    sirms.main_params(in_fname=input_sdf_file,    # input
                      out_fname=x_fname,        # output
                      opt_diff=properties,
                      min_num_atoms=2,
                      max_num_atoms=4,
                      min_num_components=1,
                      max_num_components=2,
                      min_num_mix_components=2,
                      max_num_mix_components=2,
                      mix_fname=None,
                      descriptors_transformation='num',
                      mix_type='abs',
                      opt_mix_ordered=False,
                      opt_verbose=False,
                      opt_noH=True,
                      frag_fname=fragments_fname,
                      per_atom_fragments=False,
                      self_association_mix=False,
                      reaction_diff=False,
                      quasimix=False,
                      id_field_name=None,
                      output_format=output_format,
                      ncores=n_cores)

    # filter sirms descriptors
    filter_descriptors.main_params(in_fname=x_fname,
                                   out_fname=x_fname,
                                   file_format=output_format)

def predict_properties(parameters: List, fragments_fname: str, output_format: str):
    """
    Creates summarized file with predictions

    :parama parameters: list of dicts with parameters
    :parama fragmens_fname: path to file with calculated descriptors
    :parama output_format: svm/txt/...
    """

    for parameter in parameters:
        print("Prediction for {} started".format(parameter['name']))
        output_file_name = os.path.join(os.path.dirname(fragments_fname),
                                        'predictions_{}.txt'.format(parameter['name']))
        predict.main_params(x_fname=fragments_fname,
                            input_format=output_format,
                            out_fname=output_file_name,
                            model_names=parameter['types_of_alg'],
                            model_dir=parameter['path'],
                            model_type=parameter['type_of_model'],
                            ad=['bound_box'],
                            verbose=False,
                            title=parameter['name'])
