#!/usr/bin/env python

import sys
import os
from os.path import exists, isfile
import yaml

from typing import Dict
from typing import TextIO
from typing import List


CONFIG_STRUCTURE = [['script_dir', 'path_to_folder_of_optimizer'],
                    ['working_dir', 'path_to_output_dir'],
                    ['spci_dir', 'path_to_folder_with_spci'],
                    ['sirms_dir', 'path_to_folder_with_sirms'],
                    ['num_parameters', 'number_of_parameters_to optimize'],
                    ['setup_file', 'path_to_setup_file'],
                    ['std_rules', 'path_to_file_with_std_rules'],
                    ['chemaxon', 'path_to_chemaxon_bin_folder'],
                    ['seed_structure', 'path_to_seed_structure'],
                    ['number_of_selected_compounds', 'fill'],
                    ['bounded_box', 'True or False'],              # True
                    ['properties_chemaxon', 'fill'],      # 'charge logp acc don refractivity'
                    ['properties_sirms', 'fill'],         # 'CHARGE LOGP HB REFRACTIVITY'
                    ['properties_calc_contrib', 'fill'],  # 'overall'
                    ['smart_string', "'[#6+0;!$(*=,#[!#6])]!@!=!#[*]'"],
                    ['max_cuts', 'fill'],
                    ['radius', 'fill'],
                    ['replacement_database', 'path_to_database_with_replacement'],
                    ['number_of_worst_fragments', 'fill'],
                    ['max_frag_size', 'fill'],
                    ['output_format', 'svm'],
                    ['num_of_generation', 'fill'],
                    ['optimization_methods', 'fill'],     # 'pareto desirability'
                    ['store_all_files', 'True or False']          # If false, it deletes all temp files, only db will be stored
                    ]

PARAMETER_STRUCTURE = [['name', 'name_of_parameter'],
                       ['path', 'path_to_parameter_model_folder'],
                       ['types_of_alg', 'space_separated_list_of_all_algs_used_for_predictions'],
                       ['type_of_model', 'class or reg'],
                       ['threshold', 'fill'],
                       ['range', 'fill'],
                       ['desirability', 'fill']
                       ]

def create_config(output_dir: str, num_of_parametrs: int, file_name: str='config.yaml') -> None:
    """
    Create empty file for input settings for optimizer

    :param output_dir: path to output dir
    :param num_of_parametrs: number of parameters to optimize
    :param file_name: name of config file
    """

    assert exists(output_dir), "Output directory doesn't exist"

    # set indexes
    script_dir_index = 0
    working_dir_index = 1
    spci_dir_index = 2
    sirms_dir_index = 3
    num_parameters_index = 4

    CONFIG_STRUCTURE[script_dir_index][1] = os.path.abspath(os.path.dirname(sys.argv[0]))
    # working dir
    CONFIG_STRUCTURE[working_dir_index][1] = os.path.abspath(output_dir)
    # spci dir
    CONFIG_STRUCTURE[spci_dir_index][1] = os.path.join(CONFIG_STRUCTURE[script_dir_index][1], 'spci')
    # sirms dir
    CONFIG_STRUCTURE[sirms_dir_index][1] = os.path.join(CONFIG_STRUCTURE[spci_dir_index][1], 'sirms')
    # num param
    CONFIG_STRUCTURE[num_parameters_index][1] = num_of_parametrs

    conf_file = os.path.join(output_dir, file_name)
    with open(conf_file, 'w') as config:
        for item in CONFIG_STRUCTURE:
            config.write("{}: {}\n".format(item[0], item[1]))
        for i in range(num_of_parametrs):
            config.write("param_{}: \n".format(str(i)))
            for item in PARAMETER_STRUCTURE:
                config.write("  {}: {}\n".format(item[0], item[1]))
    print("Structure of configuration file was created in {} directory. "
          "Fell free to rewrite default values.".format(CONFIG_STRUCTURE[working_dir_index][1]))

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
            config = yaml.load(stream)
        except yaml.YAMLError as exc:
            print(exc)

    # check non parameters settings
    for key, value in config.items():
        # check directories
        if (key == 'script_dir') or (key == 'working_dir') \
            or (key == 'spci_dir') or (key == 'sirms_dir') \
            or (key == 'chemaxon'):
            assert exists(config[key]), "{} doesn't exists".format(key)
            config[key] = value
        # check files
        elif (key == 'setup_file') \
            or (key == 'seed_structure') \
            or (key == 'std_rules') \
            or (key == 'replacement_database'):
                assert isfile(config[key]), "{} doesn't exists".format(key)
                config[key] = value
        # check param
        elif "param_" in key:
            # check if models dir exists
            assert exists(config[key]['path']), "{} doesn't exists".format(config[key]['path'])
            # check if types of models exists, e.g. model + type + ".pkl"
            config[key]['types_of_alg'] = config[key]['types_of_alg'].split(" ")
            for alg in config[key]['types_of_alg']:
                assert isfile(os.path.join(config[key]['path'], alg + ".pkl")), \
                        "{} doesn't exists".format(os.path.join(config[key]['path'], alg + ".pkl"))
        # check settings with numbers
        elif (key == 'number_of_selected_compounds') or (key == 'max_cuts') \
            or (key == 'radius') or (key == 'number_of_worst_fragments') \
            or (key == 'max_frag_size') or (key == 'num_of_generation') \
            or (key == 'num_parameters'):
                try:
                    num = int(value)
                except:
                    num = value
                assert type(num) == type(1), "{} is not defined properly".format(key)
                config[key] = num
        # strip properties for chemaxon and sirms
        elif (key == 'properties_chemaxon') or (key == 'properties_sirms'):
            config[key] = value.split(" ")

    config['bounded_box'] = True if config['bounded_box'] == 'True' or config['bounded_box'] == 'true' else False
    config['store_all_files'] = True if config['store_all_files'] == 'True' or config['store_all_files'] == 'true' else False

    # print(138 * "_")
    # for key, value in config.items():
    #     if type(value) == type({}):
    #         print("{}:".format(key))
    #         for key_nested, value_nested in value.items():
    #             print("  {}: {}".format(key_nested, value_nested))
    #     else:
    #         print("{}: {}".format(key, value))
    # print(138 * "_")

    return config
