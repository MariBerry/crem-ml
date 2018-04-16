#!/usr/bin/env python

import argparse
import os
import shutil
import sys

from typing import Dict

import datetime

import process_config
import optimizer_utils
import process_predictions


def optimize(settings: Dict, input_config: str) -> None:
    """
    Run all optimization tasks

    :param settings: dictionary with input settings
    :param input_config: path to input config file
    """
    # for new IDs
    num_of_compounds = 0

    # for fitted compounds
    num_of_fitted_compounds = 0

    parameters_list_dicts = [settings[parameter] for parameter in settings if "param_" in parameter]

    # create database
    settings['output_database'] = optimizer_utils.create_database(
        settings['working_dir'], [parameter['name'] for parameter in parameters_list_dicts])

    shutil.copyfile(input_config, os.path.join(settings['working_dir'], 'config.yaml'))

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
                                         settings['seed_structure'],
                                         settings['output_database'],
                                         gen)
        # check if we have new compounds in new generation
        if tmp_num_comp == num_of_compounds:
            print("\nIn generation {} aren't new compounds".format(gen))
            sys.exit()

        # copy input_sdf_file into gen directory
        new_sdf = os.path.join(generation_dir, 'input_dataset.sdf')
        shutil.copyfile(settings['seed_structure'], new_sdf)
        settings['seed_structure'] = new_sdf

        # start generation
        start = datetime.datetime.now()
        print(50 * '_', '\nGeneration {}: {}'.format(gen, start))

        # standardization
        settings['seed_structure'] = optimizer_utils.standardize_sdf(
            input_sdf_file=settings['seed_structure'],
            std_rules_path=settings['std_rules'],
            chemaxon_path=settings['chemaxon']
        )

        # calc atomic properties
        settings['seed_structure'] = optimizer_utils.calculate_atomic_prop(
            input_sdf_file=settings['seed_structure'],
            chemaxon_path=settings['chemaxon'],
            properties=settings['properties_chemaxon']
        )

        # calculation of sirms descriptors
        optimizer_utils.calculate_sirms_descriptors(settings['seed_structure'],
                                    settings['setup_file'],
                                    settings['properties_sirms'],
                                    settings['output_format'],
                                    settings['n_cores']
                                    )
        fragments_fname = os.path.join(generation_dir, 'x.txt')

        # predict properties of std_lbl_sdf file
        optimizer_utils.predict_properties(parameters_list_dicts,
                                           fragments_fname,
                                           settings['output_format']
                                           )

        # process predictions
        list_of_prediction_files = []   # prepare list of file paths with predictions
        for parameter in parameters_list_dicts:
            list_of_prediction_files.append(
                os.path.join(generation_dir, 'predictions_{}.txt'.format(parameter['name'])))

        settings['processed_predictions_file'] = os.path.join(generation_dir, 'processed_predictions.txt')
        num_of_fitted_compounds += process_predictions.main(settings['seed_structure'],
                                     list_of_prediction_files,
                                     settings['processed_predictions_file'],
                                     [parameter['name'] for parameter in parameters_list_dicts],
                                     settings['optimization_methods'],
                                     [parameter['threshold'] for parameter in parameters_list_dicts],
                                     settings['bounded_box'],
                                     [parameter['desirability'] for parameter in parameters_list_dicts],
                                     settings['number_of_selected_compounds']
                                     )

        if num_of_fitted_compounds >= settings['num_output_compounds']:
            print("Optimizer reached number of fitted compounds specified in config.")

        # find fragments
        settings['fragments_ids_file'] = os.path.join(generation_dir, 'fragments_ids.txt')
        error_fname_frag = os.path.join(generation_dir, 'fragments_log.log')
        optimizer_utils.find_frags_rdkit(settings['processed_predictions_file'],
                                         settings['fragments_ids_file'],
                                         settings['smart_string'],
                                         settings['max_cuts'],
                                         settings['radius'],
                                         settings['keep_stereo'],
                                         error_fname_frag)



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
        optimize(settings, args['input_config'])

if __name__ == '__main__':
    main()
