import os

from subprocess import call
import sqlite3 as lite

from typing import List


def create_database(working_dir: str, parameter_to_optimize: List) -> str:
    """
    Define new database and save it to working_dir.

    :param working_dir: path to working directory, new db will be stored there
    :param parameter_to_optimize: list of all parameters
    :return: path_to_database
    """
    parameter_to_optimize = ["_{}".format(parameter) for parameter in parameter_to_optimize]
    path_to_database = os.path.join(working_dir, 'output.db')

    table_str = "CREATE TABLE optimizer_table"\
                    "(id TEXT NOT NULL,"\
                    "smi TEXT NOT NULL,"\
                    "generation INTEGER NOT NULL,"\
                    "parent TEXT,"\
                    "transformation TEXT,"\
                    "fit INTEGER,"
    table_str += " REAL,".join(parameter_to_optimize) + " REAL)"

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
