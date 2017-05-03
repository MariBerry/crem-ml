#!/usr/bin/env python

from process_predictions import parse_threshold
import math
import numpy as np


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
        first_line += ''.join(['norm_contrib_'+parameter+'\t' for parameter in parameters])
        first_line += 'average\n'
        f.write(first_line)
    return f


def get_predicted_value_for_whole_compound(in_file, parameters):
    predicted_values = {}
    predicted_parameters = {}

    with open(in_file, 'r') as in_f:
        iter_file = iter(in_f)
        for line in iter_file:
            if 'ID' in line.rstrip():
                id = next(iter_file).rstrip()
                for parameter in parameters:
                    while not parameter in next(iter_file).rstrip():
                        pass
                    predicted_parameters[parameter] = float(next(iter_file).rstrip())
                predicted_values[id] = predicted_parameters
        return predicted_values


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


def normalize_contributions(in_fname, parameters, predicted_values, thresholds, ranges):
    with open(in_fname, 'a') as in_file:
        list_of_contributions = list()
        for number_of_parameter, parameter in enumerate(parameters):
            with open('contrib_' + parameter + '.txt', 'r') as contrib_f:
                string_to_write_to_file = ''
                contrib_f.readline()

                for line_number, line in enumerate(contrib_f.readlines()):
                    line = line.split('\t')

                    norm_value = compute_normalized_value(
                        float(line[5]),                         # fragment contributions
                        predicted_values[line[0]][parameter],   # predicted parameter value of compound
                        thresholds[number_of_parameter],        # threshold for parameter
                        ranges[parameter])                      # range of parameter

                    if number_of_parameter == 0:
                        list_of_contributions.append([norm_value])
                    else:
                        list_of_contributions[line_number].append(norm_value)

                    # if it is last parameter, we need to compute average of normalized values and construct string
                    # which will be added to file
                    if number_of_parameter == len(parameters) - 1:
                        string_to_write_to_file = line[0] + '\t' + line[2] + '\t'
                        string_to_write_to_file += \
                            ''.join([str(norm_value) + '\t' for norm_value in list_of_contributions[line_number]])
                        average = sum(list_of_contributions[line_number])/len(list_of_contributions[line_number])
                        string_to_write_to_file += str(average)

                        in_file.write(string_to_write_to_file + '\n')
                        string_to_write_to_file = ''
                number_of_fragments = line_number + 1
    return line_number



def pick_worst_ones(in_fname, number_of_fragments):
    number_of_col = 3
    worst_list = np.zeros((number_of_fragments + 1, number_of_col))
    with open(in_fname, 'r') as in_f:
        in_f.readline()
        for i in range(number_of_fragments):
            line = in_f.readline().split('\t')
            worst_list[i] = np.array([line[0], line[1], float(line[-1])])
        worst_list = worst_list[worst_list[:, number_of_col-1].argsort()]

        for line in in_f.readlines():
            line = line.split('\t')
            worst_list[number_of_fragments] = np.array([line[0], line[1], float(line[-1])])
            worst_list = worst_list[worst_list[:, number_of_col-1].argsort()]

    return worst_list[:number_of_fragments]


in_fname = 'fragment_contrib_norm.txt'
parameters = ['LOGBB', 'solubility']
pred_file = 'output_pareto.sdf'
thresholds =['more0.5', 'more-4']

processed_file = create_file_with_processed_data(in_fname, parameters)
predicted_values = get_predicted_value_for_whole_compound(pred_file, parameters)
thresholds = parse_threshold(thresholds)
number_of_fragments = normalize_contributions(in_fname, parameters, predicted_values, thresholds, get_ranges())

# we can't specify more worst fragments then we have
number_of_worst_fragments = 50

worst_list = pick_worst_ones(in_fname, number_of_worst_fragments)


