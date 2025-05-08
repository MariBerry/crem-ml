#!/usr/bin/env python

import sys
import os
from os.path import exists, isfile
import yaml

from typing import Dict
from typing import TextIO
from typing import List


CONFIG_STRUCTURE = [['working_dir', 'path_to_output_dir'],
                    ['seed_structure', 'path_to_seed_structure'],
                    ['num_of_generations', 'number of generations'],
                    ['num_output_compounds', 'number of  compounds to generate'],
                    ['n_cores', '1'],
                    ['number_of_selected_compounds', 'number of compounds selected for optimization in one generation '],
                    ['random_compounds_selection', '0'],  # from 0 - 1, 0.2 means 20% of selected compounds are chosen randomly (floored)
                    ['optimization_method', 'pareto or desirability'],  # ' one of: pareto desirability'
                    ['smarts_string', "'[#6+0;!$(*=,#[!#6])]!@!=!#[*]'"],
                    ['max_cuts', '1'],
                    ['protected_ids', 'None'], # field in seed sdf, containing atom ids that should not be touched by replacements (default name, or specify as ar
                    ['replacement_database', 'path_to_database_with_replacement'],
                    ['number_of_worst_fragments', ''],
                    ['random_fragments_selection', '0'], # from 0 - 1, 0.2 means 20% of selected fragments are chosen randomly (floored)
                    ['max_frag_size', '7'],
                    ['min_inc', '-2'],
                    ['max_inc', '2'],
                    ['radius', '2'],
                    ['descriptors_type','type of descriptors to use'],
                    ['bounding_box', 'True'],
                    ['multitask', 'True or False'],  # False
                    ['variance_threshold', 'None'],  # threshold for variance when calculating applicability domain
                    ]

PARAMETER_STRUCTURE = [['name', 'name_of_parameter'],
                       ['path', 'path_to_parameter_model_folder'],
                       ['types_of_alg', 'space_separated_list_of_all_algs_used_for_predictions'],
                       ['type_of_model', 'class or reg'],
                       ['threshold', ''],
                       ['range', ''],
                       ['desirability', 'None']
                       ]

def create_config(output_dir: str, num_of_parameters: int, file_name: str='config.yaml') -> None:
    """
    Create config file for input settings for optimizer; set useful defaults

    :param output_dir: path to output dir
    :param num_of_parameters: number of parameters to optimize
    :param file_name: name of config file
    """

    assert exists(output_dir), "Output directory doesn't exist"

    # set indexes
    working_dir_index = 0
    num_parameters_index = 1

    # working dir
    CONFIG_STRUCTURE[working_dir_index][1] = os.path.abspath(output_dir)

    conf_file = os.path.join(output_dir, file_name)
    with open(conf_file, 'w') as config:
        for item in CONFIG_STRUCTURE:
            config.write("{}: {}\n".format(item[0], item[1]))
        for i in range(num_of_parameters):
            config.write("param_{}: \n".format(str(i)))
            for item in PARAMETER_STRUCTURE:
                config.write("  {}: {}\n".format(item[0], item[1]))
    print("Structure of configuration file was created in {} directory. "
          "Please, fill in non-default values. Feel free to rewrite default values.".format(CONFIG_STRUCTURE[working_dir_index][1]))

def str_to_bool(val):
    return str(val).strip().lower() in ("1", "true", "yes", "on")
def parse_none(val):
    return None if str(val).strip().lower() in ("none", "", "null", "nil") else val

def test_config(input_config: str) -> Dict:
    """
    Test and process all settings in config file

    :param input_config: path to config file
    :return: dictionary with checked and processed settings
    """

    # check if config file exists
    assert isfile(input_config), "{} is not a file".format(input_config)
    with open(input_config, 'r') as stream:
        try:
            config = yaml.load(stream, Loader=yaml.FullLoader)
        except yaml.YAMLError as exc:
            print(exc)

    print(set([i for i in config.keys() if 'param' not in i]) - set([item[0] for item in CONFIG_STRUCTURE]))
    assert len( set([i for i in config.keys() if 'param' not in i]) - set([item[0] for item in CONFIG_STRUCTURE])) <=0 # todo finish this with print

    # set defaults for  settings absent in config
    if 'n_cores' not in config:
        config['n_cores'] = 1
    if 'random_compounds_selection' not in config:
        config['random_compounds_selection'] = 0
    if 'smarts_string' not in config:
        config[
            'smarts_string'] = "'[#6+0;!$(*=,#[!#6])]!@!=!#[*]'"  # needs to be in quotes because of special characters
    if 'max_cuts' not in config:
        config['max_cuts'] = 1
    if 'protected_ids' not in config:
        config['protected_ids'] = None
    if 'random_fragments_selection' not in config:
        config['random_fragments_selection'] = 0
    if 'max_frag_size' not in config:
        config['max_frag_size'] = 7  # enables replacing 6-membered ring plus one atom in it
    if 'min_inc' not in config:
        config['min_inc'] = -2
    if 'max_inc' not in config:
        config['max_inc'] = 2
    if 'radius' not in config:
        config['radius'] = 2
    if 'boundnig_box' not in config:
        config['boundnig_box'] = True
    if 'multitask' not in config:
        config['multitask'] = False
    if 'variance_threshold' not in config:
        config['variance_threshold'] = None

    # check general settings
    for key, value in config.items():
        # check directories
        if (key == 'working_dir') :
            assert exists(config[key]), "{} doesn't exists".format(key)
            config[key] = value
        # check files
        elif  (key == 'seed_structure') \
            or (key == 'replacement_database'):
                assert isfile(config[key]), "{} doesn't exists".format(key)
                config[key] = value
        # check parameters
        elif "param_" in key:
            # check if models dir exists
            assert exists(config[key]['path']), "{} doesn't exists".format(config[key]['path'])
            # check if types of models exists, e.g. model + type + ".pkl"
            config[key]['types_of_alg'] = config[key]['types_of_alg'].split(" ")
            for alg in config[key]['types_of_alg']:
                if alg != "MPNN":
                    assert isfile(os.path.join(config[key]['path'], alg + ".pkl")), \
                        "{} doesn't exists".format(os.path.join(config[key]['path'], alg + ".pkl"))
            if config['optimization_method'] == 'desirability':
                   assert  ('desirability'  in config[key]) and (parse_none(config[key]['desirability']) is not None)
        # check settings with  int numbers
        elif (key == 'n_cores') or (key == 'max_cuts') \
            or (key == 'radius') or (key == 'number_of_worst_fragments') \
            or (key == 'min_inc') or (key == 'max_inc')\
            or (key == 'max_frag_size') or (key == 'num_of_generations') \
            or (key == 'num_output_compounds'):
                try:
                    num = int(value)
                except:
                    num = value
                assert isinstance(num, int), "{} is not defined properly".format(key)
                config[key] = num
        # check settings with float numbers
        elif (key == 'random_compounds_selection') or (key == 'random_fragments_selection'):
            try:
                num = float(value)
            except:
                num = value
            assert isinstance(num, float), "{} is not defined properly".format(key)
            config[key] = num
        # check settings with  bool values
        elif (key == 'bounding_box') or (key == 'multitask'):
            num = str_to_bool(value)
            assert isinstance(num, bool), "{} is not defined properly".format(key)
            config[key] = num
        # check str settings and settings with possible None
        elif (key == 'protected_ids') or (key == 'variance_threshold'):
            num = parse_none(value)
            assert (isinstance(num, str) or num  is None), "{} is not defined properly".format(key)
            config[key] = num
        # assure we use only compatible descriptors
        elif (key == 'descriptors_type'):
            acceptable_descr = ['sirms', 'MG2', 'AP', 'RDK', 'TT', 'bMG2', 'bAP', 'bRDK', 'MPNN_fingerprint']
            assert (value in acceptable_descr), "{} is not defined properly".format(key)



    # warnings
    if 'bounding_box' in config  and config ['descriptors_type'] == "MPNN_fingerprint":
        print( "Note, bounding box is ignored when MPNN models are used.")
    if 'multitask' in config and  config ['multitask']  and config ['descriptors_type'] != "MPNN_fingerprint":
        print( "Note, parameter 'multitask' is ignored  when models other than MPNN are used.")


    return config
