#!/usr/bin/env python

import argparse
import os
import shutil

from typing import Dict

import datetime

import process_config
import optimizer_utils


def optimize(settings: Dict) -> None:
    """
    Run all optimization tasks

    :param settings: dictionary with input settings
    """

    settings['working_dir'] = os.path.join(settings['working_dir'], 'out')

    for gen in range(settings['num_of_generation']):

        # create generation dir
        generation_dir = os.path.join(settings['working_dir'], 'generation_{}'.format(gen))
        if os.path.exists(generation_dir):
           shutil.rmtree(generation_dir)
        os.makedirs(generation_dir)

        # copy input sdf file into generation dir
        shutil.copy2(settings['path_to_seed_structure'], os.path.join(generation_dir, 'input_dataset.sdf'))
        settings['path_to_seed_structure'] = os.path.join(generation_dir, 'input_dataset.sdf')

        # change dir to generation dir
        os.chdir(generation_dir)

        # start generation
        start = datetime.datetime.now()
        print(50 * '_', '\nGeneration {}: {}'.format(gen, start))

        #standardization
        settings['path_to_seed_structure'] = optimizer_utils.standardize_sdf(
            input_sdf_file=settings['path_to_seed_structure'],
            std_rules_path=settings['path_to_std_rules_file'],
            chemaxon_path=settings['path_to_chemaxon_bin']
        )






def main():
    parser = argparse.ArgumentParser(description='System for designing new drugs')
    parser.add_argument('-i', '--input_config',
                        help='path to config with input settings')
    parser.add_argument('-o', '--output_location', required=True,
                        help='path to output location, place where all outputs are saved')
    parser.add_argument('-d', '--define_config_structure', action='store_true', default=False,
                        help='define config structure and save it to output location')
    args = vars(parser.parse_args())

    if args['define_config_structure']:
        process_config.create_config(args['output_location'])
    else:
        settings = process_config.check_config(args['input_config'])
        optimize(settings)


if __name__ == '__main__':
    main()
