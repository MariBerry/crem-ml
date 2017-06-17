#!/usr/bin/env python
import os
import sys
import shutil

from subprocess import call
sys.path.insert(1, os.path.join(sys.path[0], 'spci'))

import calc_atomic_properties_chemaxon
import filter_descriptors
import predict
import find_frags_auto_rdkit as find_frags
import calc_frag_contrib as calc_contrib
import process_predictions
import process_contributions_rdkit
import frag_replacements

sys.path.insert(1, os.path.join(sys.path[0], 'spci/sirms'))
import sirms

import datetime

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


# calculate sirms descriptors
# input format: working_file: path to standardized and labeled sdf file
#               setup_path: path to file with setup for calculation of sirms descriptors
#               properties: list of properties, e.g. ['charge', 'refractivity', 'logp', ...]
#               output_format: string, e.g. 'txt'/'svm'
#               copy_setup: boolean value, if true, it copies file with setup to calculation of sirms descriptor
#                           to directory with input sdf file

def calculate_sirms_descriptors(working_file, setup_path, properties, output_format, ncores, fragments_fname=None,
                                copy_setup=True):
    print("Descriptors calculation started. Please wait it can take some time")

    # copy setup.txt to folder with sdf file
    if copy_setup:
        shutil.copyfile(
            setup_path, os.path.join(os.path.dirname(working_file), setup_path.split("/")[-1]))

    # calc sirms descriptors
    if fragments_fname is not None:
        x_fname = os.path.join(os.path.dirname(working_file), 'new_x.txt')
    else:
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
                      frag_fname=fragments_fname,
                      per_atom_fragments=False,
                      self_association_mix=False,
                      reaction_diff=False,
                      quasimix=False,
                      id_field_name=None,
                      output_format=output_format,
                      ncores=ncores)

    # filter sirms descriptors
    filter_descriptors.main_params(in_fname=x_fname,
                                   out_fname=x_fname,
                                   file_format=output_format)


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


def find_frags_rdkit(input_sdf_file, fragment_ids_file, smarts_string, max_cuts, radius, keep_stereo, verbose, error_fname):
    print("Finding fragments has started")
    find_frags.main_params(in_sdf=input_sdf_file,
                                out_txt=fragment_ids_file,
                                query=smarts_string,
                                max_cuts=max_cuts,
                                radius = radius,
                                keep_stereo = keep_stereo,
                                verbose=verbose,
                                error_fname=error_fname)


def calc_frag_contrib(x_fname, parameters, models, models_dir, properties, models_type, in_format):

    for parameter, model, model_dir, model_type in zip(parameters, models, models_dir, models_type):
        print("Fragment contribution for {} started".format(parameter))
        calc_contrib.main_params(x_fname=x_fname,
                                      out_fname='contrib_' + parameter + '.txt',
                                      model_names=model,
                                      model_dir=model_dir,
                                      prop_names=properties,
                                      model_type=model_type,
                                      verbose=False,
                                      save_pred=True,
                                      input_format=in_format,
                                      long_format=True,
                                      save_frag_ids=True)
        """
        f = open('predictions.txt', 'a')
        f.write(open('predictions_' + parameter + '.txt', 'r').read())
        f.close()
        os.remove('predictions_' + parameter + '.txt')
        """

def process_contributions(input_sdf, frag_norm_output_file, worst_output_file, parameters_to_predict, models, thresholds,
                number_of_worst_fragments):
    print("Processing contributions has started")
    process_contributions_rdkit.main_params(input_sdf,
                                        frag_norm_output_file,
                                        worst_output_file,
                                        parameters_to_predict,
                                        models,
                                        thresholds,
                                        number_of_worst_fragments)

def fragment_replacements(input_sdf, input_worst, input_ids, input_connection_db, output_product_file):
    print("Fragments replacements has started")
    frag_replacements.main_params(input_sdf,
                                    input_worst,
                                    input_ids,
                                    input_conn_db,
                                    output_product_file)


# number of generation
num_of_gen = 3

pathname = os.path.dirname(sys.argv[0])
script_dir = os.path.abspath(pathname)
working_dir = os.getcwd()
chemaxon_path = "/home/david/ChemAxon/JChem/bin"
input_sdf_file = working_dir + "/dataset_FW.sdf"
home_dir = working_dir
setup_file = script_dir + "/setup.txt"
std_rules_file = script_dir + "/std_rules.xml"

for gen in range(num_of_gen):
    generation_dir = working_dir+'/generation_'+str(gen)
    if os.path.exists(generation_dir):
       shutil.rmtree(generation_dir)
    os.makedirs(generation_dir)
    shutil.copy2(input_sdf_file, generation_dir+'/')


    if gen == 0:
        os.rename(generation_dir + '/' + input_sdf_file.split('/')[-1], generation_dir + '/input_dataset.sdf')
    else:
        os.rename(input_sdf_file, generation_dir + '/input_dataset.sdf')

    input_sdf_file = generation_dir + '/input_dataset.sdf'
    working_dir = generation_dir
    os.chdir(working_dir)

    std_sdf_file = working_dir + "/" + input_sdf_file.split("/")[-1].split(".")[0] + '_std.sdf'
    std_lbl_sdf_file = working_dir + "/" + input_sdf_file.split("/")[-1].split(".")[0] + '_std_lbl.sdf'
    x_fname = working_dir + "/x.txt"
    predictions = working_dir + "/predictions.txt"
    fragment_ids_file = working_dir + '/fragment_ids.txt'
    fragment_x_fname = working_dir + '/new_x.txt'

    properties_chemaxon = ['charge', 'logp', 'acc', 'don', 'refractivity']
    properties_sirms = ['CHARGE', 'LOGP', 'HB', 'REFRACTIVITY']
    properties_calc_contrib = ['overall']
    output_format = 'svm'
    print(50*'_')
    print('Generation: ', gen)

    start = datetime.datetime.now()
    print(start)

    # standardization
    standardize_sdf(input_sdf_file, std_rules_file)

    # calculation of atomic properties
    calculate_atomic_prop(std_sdf_file, chemaxon_path, properties_chemaxon)

    # calculation of sirms descriptors
    ncores = 8
    calculate_sirms_descriptors(std_lbl_sdf_file, setup_file, properties_sirms, output_format, ncores)

    # with more properties to predict
    paramaters_to_predict = ['LOGBB', 'solubility']
    models_dir = [script_dir + '/models/logBB/', script_dir + '/models/solubility/']
    models = [['gbm', 'svm'], ['rf', 'svm']]
    models_type = ['class', 'reg']

    # predict properties of std_lbl_sdf file
    predict_properties(paramaters_to_predict, x_fname, output_format, models, models_dir, models_type)

    # process predictions
    # methods = ['filtering', 'pareto', 'desirability']
    methods = ['desirability']
    thresholds =['more0.5', 'more-2', 'desirability_0.45:0,0.55:10*x-4.5,1000:1_-2.1:0,-1.9:5*x+10.5,1000:1#5']
    # thresholds =['more0.5', 'more-2']
    # for testing bounded_box = False
    bounded_box = True
    process_prediction(std_lbl_sdf_file, predictions, paramaters_to_predict, 'output.sdf', methods, thresholds, bounded_box)

    # input to modification part
    output_sdf_file = os.path.dirname(input_sdf_file) + '/output_process_predictions.sdf'

    # calculate fragment ids
    smarts_string = "[#6+0;!$(*=,#[!#6])]!@!=!#[*]"
    max_cuts = 3
    find_frags_verbose = False
    error_fname = os.path.dirname(input_sdf_file) + '/fragments_log.log'
    radius = [3]
    keep_stereo = False
    find_frags_rdkit(output_sdf_file, fragment_ids_file, smarts_string, max_cuts, radius, keep_stereo, find_frags_verbose, error_fname)

    # calculate sirms descriptors using fragments
    calculate_sirms_descriptors(output_sdf_file, setup_file, properties_sirms,
                                output_format, ncores, fragments_fname=fragment_ids_file)

    # calculate fragments contributions
    calc_frag_contrib(fragment_x_fname, paramaters_to_predict, models, models_dir, properties_calc_contrib, models_type,
                      output_format)

    # find worst fragments
    models_contrib = ['gbm_svm', 'rf_svm']  # different format of models, because of calling from terminal
    number_of_worst_fragments = 1
    process_contributions(output_sdf_file, 'fragment_contrib_norm.txt', 'worst_fragments.txt', paramaters_to_predict,
                                 models_contrib, thresholds, number_of_worst_fragments)

    # fragment replacements
    input_sdf = 'output_process_predictions.sdf'
    input_worst = 'worst_fragments.txt'
    input_ids = 'fragment_ids.txt'
    output_product = 'new_compounds.sdf'
    input_conn_db = '/home/david/Documents/projects/dp/optimizer/replacement_chembl_cuts4_H.db'

    fragment_replacements(input_sdf, input_worst, input_ids, input_conn_db, output_product)

    end = datetime.datetime.now()
    print(end)
    print(end - start)

    input_sdf_file = os.path.abspath('new_compounds.sdf')
    # end - for creating new folder with new generation we have to prepare some path variables
    os.chdir(home_dir)
    working_dir = home_dir
