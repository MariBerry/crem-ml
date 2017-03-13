#!/usr/bin/env python
import os
import sys
import shutil

from subprocess import call

import calc_atomic_properties_chemaxon
import filter_descriptors
import predict
import process_predictions

sys.path.insert(1, os.path.join(sys.path[0], 'sirms'))
import sirms


def quote_str(s):
    return "'%s'" % s


# standardization, it creates file with standardized compounds
# input format: input_sdf_file: path to sdf file with compounds
#               std_rules_path: path to file with rules for standardization
#               copy_rules: boolean value, if true, it copies std_rules_file to directory with input sdf file
def standardize_sdf(input_sdf_file, std_rules_path, copy_rules=True):
    print('Standardization is in progress...')

    # copy xml-rules if the file is absent in the sdf folder
    if copy_rules:
        shutil.copyfile(
            std_rules_path, os.path.join(os.path.dirname(input_sdf_file), std_rules_path.split("/")[-1]))

    # run standardize
    std_sdf = input_sdf_file.split("/")[-1].split(".")[0] + '_std.sdf'
    run_params = ['standardize',
                  '-c',
                  quote_str(std_rules_path),  # path to rules
                  quote_str(input_sdf_file),  # path to input file
                  '-f',
                  'sdf',  # type of output file
                  '-o',
                  quote_str(std_sdf)]  # name of output file
    call(' '.join(run_params), shell=True)
    print('Standardization finished')


# calculate atomic properties with Chemaxon, it creates file with labeled compounds
# input format: input_std_sdf: path to sdf file with standardized compounds
#               chemaxon_path: path to bin directory in chemaxon
#               properties: list of properties, e.g. ['charge', 'refractivity', 'logp', ...]
def calculate_atomic_prop(input_std_sdf, chemaxon_path, properties):
    print('Atomic properties calculation is in progress...')
    output_std_lbl_sdf = input_std_sdf.split("/")[-1].split(".")[0] + "_lbl.sdf"
    calc_atomic_properties_chemaxon.main_params(input_std_sdf,
                                                output_std_lbl_sdf,
                                                properties,
                                                None,
                                                os.path.join(chemaxon_path, 'cxcalc'))
    print('Atomic properties calculation finished')


# calculate sirms descriptors
# input format: working_file: path to standardized and labeled sdf file
#               setup_path: path to file with setup for calculation of sirms descriptors
#               properties: list of properties, e.g. ['charge', 'refractivity', 'logp', ...]
#               output_format: string, e.g. 'txt'/'svm'
#               copy_setup: boolean value, if true, it copies file with setup to calculation of sirms descriptor
#                           to directory with input sdf file

def calculate_sirms_descriptors(working_file, setup_path, properties, output_format, copy_setup=True):
    print("Descriptors calculation started. Please wait it can take some time")

    # copy setup.txt to folder with sdf file
    if copy_setup:
        shutil.copyfile(
            setup_path, os.path.join(os.path.dirname(working_file), setup_path.split("/")[-1]))

    # calc sirms descriptors
    x_fname = os.path.join(os.path.dirname(working_file), 'x.txt')
    sirms.main_params(in_fname=working_file,    # input
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
                      frag_fname=None,
                      per_atom_fragments=False,
                      self_association_mix=False,
                      reaction_diff=False,
                      quasimix=False,
                      id_field_name=None,
                      output_format=output_format)

    # filter sirms descriptors
    filter_descriptors.main_params(in_fname=x_fname,
                                   out_fname=x_fname,
                                   file_format=output_format)
    print("Descriptors calculation finished")


# it creates file with predicted values
# input format: parameters: list of parameters which are going to be predicted, e.g. ['logBB', 'solubility', ...]
#               x_fname: path to file with calculated descriptors
#               input_format: string, e.g. 'txt'/'svm'
#               models: type of models in list of lists, but in the same order as parameters,
#                       e.g. [['rf', 'gbm', ...], ['gbm', 'svm', ...], [...], ...]
#               models_dir: list of paths to directories where models are stored,
#                           e.g. ['some_path/models/logBB', 'some_path/models/solubility', ...]
#               models_type: list of types for models, but in the same order as parameters,
#                            e.g. ['reg', 'class', 'class', ...]
# output format: it creates files with predictions for specified parameters
def predict_properties(parameters, x_fname, input_format, models, models_dir, models_type):
    try:
        os.remove('predictions.txt')
    except OSError:
        pass
    for parameter, model, model_dir, model_type in zip(parameters, models, models_dir, models_type):
        print("Prediction for {} started".format(parameter))
        predict.main_params(x_fname=x_fname,                                # input file with descriptors
                            input_format=input_format,
                            out_fname='predictions_' + parameter + '.txt',   # output files with predictid values
                            model_names=model,                              # model for prediction of property
                            model_dir=model_dir,
                            model_type=model_type,
                            ad=['bound_box'],
                            verbose=False)
        print("Prediction for {} finished".format(parameter))
        f = open('predictions.txt', 'a')
        f.write(open('predictions_' + parameter + '.txt', 'r').read())
        f.close()
        os.remove('predictions_' + parameter + '.txt')


def process_prediction(input_sdf, input_pred, parameters, output_file, methods, thresholds, use_bounded_box=True):
    print("Processing prediction has started")
    process_predictions.main_params(input_sdf,
                                    input_pred,
                                    parameters,
                                    output_file,
                                    methods,
                                    thresholds,
                                    use_bounded_box)
    print("Processing prediction is finished")


pathname = os.path.dirname(sys.argv[0])
# paths
script_dir = os.path.abspath(pathname)
working_dir = os.getcwd()
chemaxon_path = "/home/david/ChemAxon/JChem/bin"                                                    # have to be specified
input_sdf_file = working_dir + "/dataset_FW.sdf"                                                    # have to be specified
setup_file = script_dir + "/setup.txt"                                                              # derived
std_rules_file = script_dir + "/std_rules.xml"                                                      # derived
std_sdf_file = working_dir + "/" + input_sdf_file.split("/")[-1].split(".")[0] + '_std.sdf'         # derived
std_lbl_sdf_file = working_dir + "/" + input_sdf_file.split("/")[-1].split(".")[0] + '_std_lbl.sdf' # derived
x_fname = working_dir + "/x.txt"                                                                    # derived
predictions = working_dir + "/predictions.txt"

properties = ['charge', 'logp', 'refractivity']
output_format = 'svm'

# standardization
standardize_sdf(input_sdf_file, std_rules_file)

# calculation of atomic properties
calculate_atomic_prop(std_sdf_file, chemaxon_path, [property_name.lower() for property_name in properties])


# calculation of sirms descriptors
calculate_sirms_descriptors(std_lbl_sdf_file, setup_file, [property_name.upper() for property_name in properties], output_format)

# preparing settings predictions for one property - example
"""
prediction_properties = ['LOGBB']
models_dir = [script_dir + '/models/logBB/']
models = [['gbm', 'svm']]
models_type = ['class']
"""

# with more properties to predict
paramaters_to_predict = ['LOGBB', 'solubility']
models_dir = [script_dir + '/models/logBB/', script_dir + '/models/solubility/']
models = [['gbm', 'svm'], ['rf', 'svm']]
models_type = ['class', 'reg']

# predict properties of std_lbl_sdf file
predict_properties(paramaters_to_predict, x_fname, output_format, models, models_dir, models_type)

# process predictions
methods = ['filtering', 'pareto']
thresholds =['more0.5', 'more-4']
# for testing bounded_box = False
bounded_box = False
process_prediction(std_lbl_sdf_file, predictions, paramaters_to_predict, 'output.sdf', methods, thresholds, bounded_box)

