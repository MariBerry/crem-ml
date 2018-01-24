#!/usr/bin/env python

import argparse
import os
import shutil
import sys

from typing import Dict

import datetime

import process_config
import optimizer_utils


def optimize(settings: Dict) -> None:
    """
    Run all optimization tasks

    :param settings: dictionary with input settings
    """
    # for new IDs
    num_of_compounds = 0

    # create database
    settings['output_database'] = optimizer_utils.create_database(
        settings['working_dir'], settings['parameter_to_optimize'])

    # create output folder
    settings['working_dir'] = os.path.join(settings['working_dir'], 'out')
    if os.path.exists(settings['working_dir']):
        shutil.rmtree(settings['working_dir'])

    for gen in range(settings['num_of_generation']):

        # create generation dir
        generation_dir = os.path.join(settings['working_dir'], 'generation_{}'.format(gen))
        if os.path.exists(generation_dir):
           shutil.rmtree(generation_dir)
        os.makedirs(generation_dir)

        # add new unique compounds into database
        tmp_num_comp = num_of_compounds
        num_of_compounds = optimizer_utils.add_mols_into_db(num_of_compounds,
                                         settings['path_to_seed_structure'],
                                         settings['output_database'],
                                         gen)
        # check if we have new compounds in new generation
        if tmp_num_comp == num_of_compounds:
            print("\nIn generation {} aren't new compounds".format(gen))
            sys.exit()

        # start generation
        start = datetime.datetime.now()
        print(50 * '_', '\nGeneration {}: {}'.format(gen, start))

        # #standardization
        # settings['path_to_seed_structure'] = optimizer_utils.standardize_sdf(
        #     input_sdf_file=settings['path_to_seed_structure'],
        #     std_rules_path=settings['path_to_std_rules_file'],
        #     chemaxon_path=settings['path_to_chemaxon_bin']
        # )






def main():
    parser = argparse.ArgumentParser(description='System for designing new drugs')
    parser.add_argument('-i', '--input_config',
                        help='path to config with input settings')
    parser.add_argument('-o', '--output_location', required=True,
                        help='path to output location, place where all outputs are saved')
    parser.add_argument('-n', '--n_params', action='store', type=int,
                        help='specifies number of parameters to optimize')
    parser.add_argument('-d', '--define_config_structure', action='store_true', default=False,
                        help='define config structure and save it to output location')
    args = vars(parser.parse_args())

    if args['define_config_structure']:
        assert args['n_params'], "You have to specify number of parameters with -n or --n_params option"
        process_config.create_config(args['output_location'], args['n_params'])
    else:
        settings = process_config.test_config(args['input_config'])
        # optimize(settings)

if __name__ == '__main__':
    main()
