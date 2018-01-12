#!/usr/bin/env python

import sys
import os
from os.path import exists, isfile
from typing import Dict
from typing import TextIO
from typing import List


CONFIG_STRUCTURE = {'script_dir': 'fill',
                    'working_dir': 'fill',
                    'spci_dir': 'fill',
                    'sirms_dir': 'fill',
                    'path_to_setup_file': 'fill',
                    'path_to_std_rules_file': 'fill',
                    'path_to_chemaxon_bin': 'fill',
                    'path_to_seed_structure': 'fill',
                    'parameter_to_optimize': 'fill',    # ['5HT1A', 'HLM']
                    'path_to_parameter_model': 'fill',  # ['path_to_5HT1A_model_dir', 'path_to_HLM_model_dir']
                    'types_of_algorithms': 'fill',      # [['rf', 'svm'], [gbm]]
                    'type_of_model': 'fill',            # ['reg', 'class']
                    'thresholds': 'fill',               # ['more7', 'less0.5']
                    'number_of_selected_compounds': 'fill',
                    'bounded_box': 'fill',              # True
                    'properties_chemaxon': 'fill',      # ['charge', 'logp', 'acc', 'don', 'refractivity']
                    'properties_sirms': 'fill',         # ['CHARGE', 'LOGP', 'HB', 'REFRACTIVITY']
                    'properties_calc_contrib': 'fill',  # 'overall'
                    'smart_string': '[#6+0;!$(*=,#[!#6])]!@!=!#[*]',
                    'max_cuts': 'fill',
                    'radius': 'fill',
                    'number_of_worst_fragments': 'fill',
                    'max_frag_size': 'fill',
                    'output_format': 'svm',
                    'num_of_generation': 'fill',
                    'optimization_methods': 'fill',     # ['desirability', 'pareto']
                    'store_all_files': 'fill'           # If false, it deletes all temp files, only db will be stored
                    }

def create_config(output_dir: str, file_name: str='config') -> None:
    """
    Create empty file for input settings for optimizer

    :param output_dir: path to output dir
    :param file_name: name of config file
    """

    assert exists(output_dir), "Output directory doesn't exist"

    CONFIG_STRUCTURE['script_dir'] = os.path.abspath(os.path.dirname(sys.argv[0]))
    CONFIG_STRUCTURE['spci_dir'] = os.path.join(CONFIG_STRUCTURE['script_dir'], 'spci')
    CONFIG_STRUCTURE['sirms_dir'] = os.path.join(CONFIG_STRUCTURE['spci_dir'], 'sirms')
    CONFIG_STRUCTURE['working_dir'] = os.path.abspath(output_dir)

    conf_file = os.path.join(output_dir, file_name)
    with open(conf_file, 'w') as config:
        for key, value in CONFIG_STRUCTURE.items():
            config.write(key + "\n")
            config.write(value + "\n\n")
    print("Structure of configuration file was created in {} directory. "
          "Fell free to rewrite default values.".format(CONFIG_STRUCTURE['working_dir']))

def test_config(config: Dict) -> Dict:
    """
    Test and process all settings in config file

    :param config: dictionary with settings
    :return: dictionary with checked and processed settings
    """
    # check if files and directories exists and parse numbers
    for key, value in config.items():
        # check directories
        if (key == 'script_dir') or (key == 'working_dir') \
            or (key == 'spci_dir') or (key == 'sirms_dir') \
            or (key == 'path_to_chemaxon_bin'):
            assert exists(config[key][0]), "{} doesn't exists".format(key)
            config[key] = value[0]
        # check files
        elif (key == 'path_to_setup_file') \
            or (key == 'path_to_seed_structure') \
            or (key == 'path_to_std_rules_file'):
                assert isfile(config[key][0]), "{} doesn't exists".format(key)
                config[key] = value[0]
        # check if models dir exists
        elif key == 'path_to_parameter_model':
            for model in config[key]:
                assert exists(model), "{} doesn't exists".format(key)
        # check if types of models exists, e.g. model + type + ".pkl"
        elif key == 'types_of_algorithms':
            for i, model in enumerate(config['path_to_parameter_model']):
                model_alg = []
                for alg in config[key][i]:
                    alg = alg.strip()
                    assert isfile(os.path.join(model, alg + ".pkl")), \
                        "{} doesn't exists".format(os.path.join(model, alg + ".pkl"))
                    model_alg.append(alg)
                config[key][i] = model_alg
        # check settings with numbers
        elif (key == 'number_of_selected_compounds') or (key == 'max_cuts') \
            or (key == 'radius') or (key == 'number_of_worst_fragments') \
            or (key == 'max_frag_size') or (key == 'num_of_generation'):
                try:
                    num = int(value[0])
                except:
                    num = value[0]
                assert type(num) == type(1), "{} is not defined properly".format(key)
                config[key] = num
        # strip properties for chemaxon and sirms
        elif (key == 'properties_chemaxon') or (key == 'properties_sirms'):
            properties = []
            for prop in config[key][0]:
                properties.append(prop.strip())
            config[key] = properties

    # check if we have all settings for all models
    n_model = set()
    n_model.add(len(config['parameter_to_optimize']))
    n_model.add(len(config['path_to_parameter_model']))
    n_model.add(len(config['types_of_algorithms']))
    n_model.add(len(config['type_of_model']))
    n_model.add(len(config['thresholds']))
    assert len(n_model) == 1, "You have to have same number of parameters for model"

    config['bounded_box'] = True if config['bounded_box'][0] == 'True' or config['bounded_box'][0] == 'true' else False
    config['store_all_files'] = True if config['store_all_files'][0] == 'True' or config['store_all_files'][0] == 'true' else False
    config['smart_string'] = config['smart_string'][0]
    config['output_format'] = config['output_format'][0]
    config['properties_calc_contrib'] = config['properties_calc_contrib'][0]

    return config

def check_config(input_config: str) -> Dict:
    """
    Check config file. If all files exist and if all settings are specified correctly.

    :param input_config: path to config file
    :return: filled CONFIG_STRUCTURE
    """

    def parse_blocks(config: TextIO) -> List:
        """
        Parse blocks of settings

        :param config: config file object
        :return: dict of settings
        """

        lines, key = {}, config.readline().strip()
        lines[key] = []

        for line in config.readlines():
            line = line.strip()
            if line == '':
                key = None
                pass
            else:
                if key:
                    if ',' in line and 'desirability' not in line and key != 'smart_string':
                        lines[key].append([item for item in line.split(',')])
                    else:
                        lines[key].append(line)
                else:
                    key = line
                    lines[key] = []

        return lines

    # check if config exists
    assert isfile(input_config), "{} is not a file".format(input_config)

    with open(input_config, 'r') as config:
        settings = parse_blocks(config)

    settings = test_config(settings)

    return settings
