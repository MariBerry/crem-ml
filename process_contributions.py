#!/usr/bin/env python

from process_predictions import parse_threshold_pareto_filtering
import math
import numpy as np
from collections import defaultdict


def get_ranges():
    """
    Now it only returns hardcoded ranges of experimental values for models. Later need to be changed to read it from
    some external filed stored in coresponding model directory.
    :return: dictionary with parameter: range
    """
    return {'LOGBB': 1, 'solubility': 13}


def create_file_with_processed_data(in_fname, parameters):
    with open(in_fname, 'w+') as f:
        first_line = 'compound_id\tfragment_id\t'
        first_line += ''.join(['norm_contrib_' + parameter + '\t' for parameter in parameters])
        first_line += 'average\n'
        f.write(first_line)
    return f


def get_predicted_value_for_whole_compound(in_file, parameters):
    number_of_compounds = []

    predicted_values = {}
    predicted_parameters = {}

    with open(in_file, 'r') as in_f:
        iter_file = iter(in_f)
        for line in iter_file:
            if 'ID' in line.rstrip():
                id = next(iter_file).rstrip()
                number_of_compounds.append(int(id))
                for parameter in parameters:
                    while not parameter in next(iter_file).rstrip():
                        pass
                    predicted_parameters[parameter] = float(next(iter_file).rstrip())
                predicted_values[id] = predicted_parameters
        return predicted_values, number_of_compounds


def compute_normalized_value(in_value, predicted_value, threshold, range):
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
    for compound_id, frag_id in frags.items():
        for frag_id, values in frag_id.items():
            property_avg = sum(values[-len(models):])/len(models)
            for i in range(len(models)):
                frags[compound_id][frag_id].pop()
            frags[compound_id][frag_id].append(property_avg)
    return frags


def get_string_of_dict(frags):
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


def normalize_contributions(in_fname, parameters, predicted_values, thresholds, ranges, models):
    with open(in_fname, 'a') as in_file:

        frags = defaultdict(dict)

        with open('contrib_' + parameters[0] + '.txt', 'r') as contrib_f:
            contrib_f.readline()
            for line in contrib_f:
                line = [l.strip() for l in line.split('\t')]

                norm_value = compute_normalized_value(
                    float(line[5]),  # fragment contributions
                    predicted_values[line[0]][parameters[0]],  # predicted parameter value of compound
                    thresholds[0],  # threshold for parameter
                    ranges[parameters[0]])  # range of parameter

                try:
                    frags[line[0]][line[2]].append(norm_value)
                except:
                    frags[line[0]][line[2]] = []
                    frags[line[0]][line[2]].append(norm_value)

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



def pick_worst_ones(in_fname, number_of_fragments, number_of_compounds):

    normalized_array = np.genfromtxt(
        in_fname,
        skip_header=1,
        usecols=(0, 1, -2),
        dtype=None,
        names=['compound_id', 'fragment_id', 'average_contrib'],
        delimiter='\t')

    normalized_array.sort(order=['compound_id', 'average_contrib'])

    worst_list = np.asarray(normalized_array[:number_of_fragments])
    number_of_compounds.remove(worst_list[0][0])

    for i in range(len(number_of_compounds)):
        for index_fragment, fragment in enumerate(normalized_array[number_of_fragments:]):
            if fragment[0] == number_of_compounds[i]:
                from_ = number_of_fragments + index_fragment
                to = 2 * number_of_fragments + index_fragment
                worst_list = np.append(worst_list, normalized_array[from_:to], axis=0)
                break

    np.savetxt('worst_fragments.txt',worst_list, delimiter='\t')

    return worst_list


in_fname = 'fragment_contrib_norm.txt'
models = [['svm'], ['rf', 'svm']]
parameters = ['LOGBB', 'solubility']
pred_file = 'output_pareto.sdf'
thresholds =['more0.5', 'more-2']

processed_file = create_file_with_processed_data(in_fname, parameters)
predicted_values, number_of_compounds = get_predicted_value_for_whole_compound(pred_file, parameters)
thresholds = parse_threshold_pareto_filtering(thresholds)
number_of_fragments = normalize_contributions(in_fname, parameters, predicted_values, thresholds, get_ranges(), models)

# we can't specify more worst fragments then we have
number_of_worst_fragments = 10

worst_list = pick_worst_ones(in_fname, number_of_worst_fragments, number_of_compounds)

print(worst_list)