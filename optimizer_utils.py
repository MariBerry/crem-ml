import os

from subprocess import call
import sqlite3 as lite
from rdkit import Chem

from typing import List


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
    table_str += " REAL,".join(parameter_to_optimize) + " REAL, prediction REAL)"

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
                mol.SetProp("parent", "None")
                mol.SetProp("transformation", "None")

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

    return num_of_compounds

def quote_str(s: str) -> str:
    """
    Quote string

    :param s: input string
    :return: quoted string
    """

    return "'%s'" % s

def standardize_sdf(input_sdf_file: str, std_rules_path: str, chemaxon_path: str, copy_rules: bool=False) -> str:
    """
    Create file with standardized compounds

    :param input_sdf_file: path to sdf file with compounds
    :param std_rules_path: path to file with rules for standardization
    :param chemaxon_path: path to chemaxon bin
    :param copy_rules: if specified, copy rules to output directory
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

    print('Standardization finished!')
    return std_sdf
