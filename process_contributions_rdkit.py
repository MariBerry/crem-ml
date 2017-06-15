#!/usr/bin/env python

from rdkit import Chem
from process_predictions import parse_threshold_pareto_filtering
import math
import numpy as np
from collections import defaultdict
import argparse


def get_ranges():
    """
    Now it only returns hardcoded ranges of experimental values for models. Later need to be changed to read it from
    some external filed stored in coresponding model directory.
    :return: dictionary with parameter: e.g. {'LOGBB': 1, 'solubility': 13}
    """
    return {'LOGBB': 1, 'solubility': 13}

def create_file_with_processed_data(file_norm_contrib, parameters):
    """
    Create file with header
    :param file_norm_contrib: name of file which you want to create
    :param parameters: optimalized parameters in list e.g. ['LOGBB', 'solubility', ...]
    """
    with open(file_norm_contrib, 'w+') as f:
        first_line = 'compound_id\tfragment_id\tfragment\t'
        first_line += ''.join(['norm_contrib_'+parameter+'\t' for parameter in parameters])
        first_line += 'average\n'
        f.write(first_line)
    return f

def get_predicted_value_for_whole_compound(in_file, parameters):
    """
    Get predicted values for whole compounds
    :param in_file: input file name, e.g. output_process_prediction.sdf
    :param parameters: optimalized parameters in list e.g. ['LOGBB', 'solubility', ...]
    :return: dict with predicted values, e.g. {21: [0.38, 1.15], 22: [0.75, -3.4]}, list of all ids
    """

    number_of_compounds = []
    predicted_values = {}
    predicted_parameters = {}
    compounds = Chem.SDMolSupplier(in_file, removeHs=False, sanitize=False)
    for mol in compounds:
        id = str(mol.GetProp('ID'))
        number_of_compounds.append(int(id))
        for parameter in parameters:
            predicted_parameters[parameter] = float(mol.GetProp(parameter))
        predicted_values[id] = predicted_parameters
        predicted_parameters = {}
    return predicted_values, number_of_compounds


def compute_normalized_value(in_value, predicted_value, threshold, range):
    """
    Get normalized value for fragment contributions
    :param in_value: fragment contributions
    :param predicted_value: predicted value for whole compound
    :param threshold: threshold values in list, e.g. ['more0.5', 'more-2']
    :param range: range of model
    :return: normalized value
    """

    if threshold[0] == 'more':
        if threshold[1] <= predicted_value:
            return abs(in_value)
        else:
            return 2 * ((1 / (1 + math.exp((7 / range) * -in_value))) - 1) + 1
    elif threshold[0] == 'less':
        if threshold[1] >= predicted_value:
            return abs(in_value)
        else:
            return -(2 * ((1 / (1 + math.exp((7 / range) * -in_value))) - 1) + 1)
    # between
    else:
        if predicted_value <= threshold[2] and predicted_value >= threshold[1]:
            return abs(in_value)
        elif predicted_value > threshold[2]:
            return -(2 * ((1 / (1 + math.exp((7 / range) * -in_value))) - 1) + 1)
        else:
            return 2 * ((1 / (1 + math.exp((7 / range) * -in_value))) - 1) + 1


def compute_average_from_models(frags, models):
    """
    Compute average from values from different models for one parameter
    :param frags: dictionary with compound id, fragment id and list with fragment name and values for all models
    :param models: list of models, e.g. [svm, rf]
    :return: frags without values for models, but with average from them
    """

    for compound_id, frag_id in frags.items():
        for frag_id, values in frag_id.items():
            property_avg = sum(values[-len(models):])/len(models)
            for i in range(len(models)):
                frags[compound_id][frag_id].pop()
            frags[compound_id][frag_id].append(property_avg)
    return frags


def get_string_of_dict(frags):
    """
    Generate formated string for output
    :param frags: dictionary with compound id, fragment id and list with fragment name and all averages
    :return: generated string formated for txt file
    """

    number_of_fragments = 0
    output_string = ""
    for compound_id, frag_id in frags.items():
        for frag_id, values in frag_id.items():
            output_string += compound_id + '\t' + frag_id + '\t'
            output_string += \
                ''.join([str(norm_value) + '\t' for norm_value in frags[compound_id][frag_id]])
            output_string += '\n'
            number_of_fragments += 1
    return output_string, number_of_fragments


def normalize_contributions(file_norm_contrib, parameters, predicted_values, thresholds, ranges, models):
    """
    Create file with normalized values
    :param file_norm_contrib: name of file which you want to create
    :param parameters: optimalized parameters in list e.g. ['LOGBB', 'solubility', ...]
    :param predicted_values: dict with predicted values, e.g. {21: [0.38, 1.15], 22: [0.75, -3.4]}
    :param thresholds: threshold values in list, e.g. ['more0.5', 'more-2']
    :param ranges: dictionary with parameter: e.g. {'LOGBB': 1, 'solubility': 13}
    :param models: list of models, e.g. [[svm, rf], [gbm]]
    :return: number of all fragments
    """

    with open(file_norm_contrib, 'a') as in_file:

        frags = defaultdict(dict)

        with open('contrib_' + parameters[0] + '.txt', 'r') as contrib_f:
            contrib_f.readline()
            for line in contrib_f:
                line = [l.strip() for l in line.split('\t')]

                norm_value = compute_normalized_value(
                    float(line[5]),                                 # fragment contributions
                    predicted_values[line[0]][parameters[0]],       # predicted parameter value of compound
                    thresholds[0],                                  # threshold for parameter
                    ranges[parameters[0]])                          # range of parameter

                try:
                    frags[line[0]][line[2]].append(norm_value)
                except:                                             # first record of fragment
                    frags[line[0]][line[2]] = []
                    frags[line[0]][line[2]].append(line[1])         # fragment name|context
                    frags[line[0]][line[2]].append(norm_value)      # normalized value for fragment contrib

            frags = compute_average_from_models(frags, models[0])

        for i, parameter in enumerate(parameters[1:]):
            file_name = 'contrib_' + parameter + '.txt'
            with open(file_name, 'r') as contrib_f:
                contrib_f.readline()
                for line in contrib_f:
                    line = [l.strip() for l in line.split('\t')]

                    norm_value = compute_normalized_value(
                        float(line[5]),  # fragment contributions
                        predicted_values[line[0]][parameter],  # predicted parameter value of compound
                        thresholds[i+1],  # threshold for parameter
                        ranges[parameter])  # range of parameter

                    frags[line[0]][line[2]].append(norm_value)
                frags = compute_average_from_models(frags, models[i + 1])

        for compound_id, frag_id in frags.items():
            for frag_id, values in frag_id.items():
                property_avg = sum(values[-len(parameters):]) / len(parameters)
                frags[compound_id][frag_id].append(property_avg)
        output_string, number_of_fragments = get_string_of_dict(frags)
        in_file.write(output_string)
    return number_of_fragments



def pick_worst_ones(file_norm_frag, file_worst_frag, number_of_fragments, number_of_compounds):
    """
    Find the worst fragments and save them into file
    :param file_worst_frag: name of file which you want to create with the worst fragments
    :param number_of_fragments: number of the worst fragments for one compounds
    :param number_of_compounds: list of all compounds (just ids)
    :return: return array of the worst fragments and save it to the file
    """

    normalized_array = np.genfromtxt(
        file_norm_frag,
        skip_header=1,
        usecols=(0, 1, 2, -2), # compound_id, fragment_id, fragment(name|context), average
        dtype=['<i8', '<i8', '<U150', '<f16'],
        names=['compound_id', 'fragment_id', 'fragment', 'average_contrib'],
        delimiter='\t')

    normalized_array.sort(order=['compound_id', 'average_contrib'])

    worst_list = np.asarray(normalized_array[:number_of_fragments])
    number_of_compounds.remove(worst_list[0][0])

    for i in range(len(number_of_compounds)):
        for index_fragment, fragment in enumerate(normalized_array[number_of_fragments:]):
            if fragment[0] == number_of_compounds[i]:
                from_ = number_of_fragments + index_fragment
                to = 2 * number_of_fragments + index_fragment

                # test if to is not to big
                if to >= normalized_array.shape[0]:
                    to = normalized_array.shape[0] - 1

                # test if to points to the same compound
                to_still_same_compound = False
                while not to_still_same_compound:
                    if fragment[0] == normalized_array[to][0]:
                        to_still_same_compound = True
                    else:
                        to -= 1

                worst_list = np.append(worst_list, normalized_array[from_:to], axis=0)
                break

    # np.savetxt(file_worst_frag, worst_list, delimiter='\t', fmt="%d %d %s %.8f")
    with open(file_worst_frag, 'w') as worst_file:
        for row in worst_list:
            worst_file.write(str(row[0]) + '\t')
            worst_file.write('\t'.join([str(item) for item in list(row)[1:]]) + '\n')

    return worst_list


def main_params(input_sdf, frag_norm_output_file, worst_output_file, parameters_to_predict, models, thresholds,
                number_of_worst_fragments):
    # prepare models to format [[svm, rf], [gbm]]
    for i, model in enumerate(models):
        models[i] = model.split('_')

    # prepare output file
    create_file_with_processed_data(frag_norm_output_file, parameters_to_predict)

    # prepare predicted values for whole compounds
    predicted_values, number_of_compounds = get_predicted_value_for_whole_compound(input_sdf, parameters_to_predict)

    # parse thresholds
    thresholds = parse_threshold_pareto_filtering(thresholds)

    # create frag_norm_contrib file and normalize all fragments
    normalize_contributions(frag_norm_output_file, parameters_to_predict, predicted_values, thresholds, get_ranges(), models)

    # create file with the worst fragments
    worst_list = pick_worst_ones(frag_norm_output_file, worst_output_file, number_of_worst_fragments, number_of_compounds)


def main():
    parser = argparse.ArgumentParser(description=
                                     'Normalize fragments contributions and pick the worst ones.')
    parser.add_argument('-is', '--in_sdf', metavar='path_to_file_with_selected_compounds', required=True,
                        help='path to file which contains selected compounds from pareto/desirability/...')
    parser.add_argument('-of', '--out_frag', metavar='fragment_contrib_norm.txt', required=True,
                        help='name of the while where you want to store normalized contributions for fragments')
    parser.add_argument('-ow', '--out_worst', metavar='worst_fragments.txt', required=True,
                        help='name of the while where you want to store the worst fragments')
    parser.add_argument('-p', '--parameters', metavar='[LOGBB solubility]', required=True, nargs='*',
                        help='predicted parameters')
    parser.add_argument('-m', '--models', metavar='[rf svm_knn]', required=True, nargs='*',
                        help='models which you used for prediction in same order as parameters connected with _')
    parser.add_argument('-t', '--thresholds', metavar='[more1.78 less12 between-2.05to2.97]', required=True, nargs='*',
                        help='thresholds to be match, written in the same order as properties')
    parser.add_argument('-n', '--n_worst', action='store', type=int,
                        help='specifies number of worst fragments')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "in_sdf": input_sdf = v
        if o == "out_frag": frag_norm_output_file = v
        if o == "out_worst": worst_output_file = v
        if o == "parameters": parameters_to_predict = v
        if o == "models": models = v
        if o == "thresholds": thresholds = v
        if o == "n_worst": number_of_worst_fragments = v

    main_params(input_sdf, frag_norm_output_file, worst_output_file, parameters_to_predict, models, thresholds,
                number_of_worst_fragments)


if __name__ == '__main__':
    main()

