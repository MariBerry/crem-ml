#!/usr/bin/env python

import argparse

import numpy as np
import pareto_simple_cull as pareto_alg

from sympy import symbols
from sympy.parsing.sympy_parser import parse_expr


# --------------------------Saving outputs--------------------------
# save output dictionary to file with predicted values
# input format: input_sdf: path to standardized and labeled sdf file which you used for prediction
#               output_file: name of output file, it will be processed according to predictions,
#                            e.g. output.sdf -> output_filtering.sdf/output_pareto.sdf
#               input_dict: output dictionary where are selected compounds stored,
#                            e.g. {'filtering': [[id1, predicted_value1, predicted_value2],[...]], 'pareto':[[...]]}
#               predictions: list of parameters name, e.g. ['logBB', 'solubility', ...]
def save_output(input_sdf, output_file, input_dict, parameters_to_predict):
    output_string = ''
    # prepare list of files
    list_of_files = []
    for parameter in input_dict.keys():
        list_of_files.append(open(output_file.split('.')[0] + '_' + parameter + '.' + output_file.split('.')[1], 'w'))
    for parameter, o_file in zip(input_dict.keys(), list_of_files):
        in_file = open(input_sdf, 'r')
        iter_file = iter(in_file)
        for item in input_dict[parameter]:
            for line in iter_file:
                line = line.rstrip()
                if str(int(item[0])) == line:
                    output_string += line + '\n'
                    for line in iter_file:
                        line = line.rstrip()
                        if '$$$$' in line:
                            for index_of_predictions, predict in enumerate(parameters_to_predict):
                                output_string += '>  <' + predict + '>\n' + str(item[1 + index_of_predictions]) + '\n\n'
                            output_string += line + '\n'
                            break
                        else:
                            output_string += line + '\n'
                    break
                else:
                   while not '$$$$' in next(iter_file).rstrip():
                           pass
        o_file.write(output_string)
        output_string = ''
        o_file.close()
        in_file.close()


# --------------------------Preparing variables--------------------------
# prepares numpy array of predicted values for compounds

# input format: input_pred: path to file with predictions,
#               bounded_box: True/False,
#               predictions: list of parameters which you want to predict, e.g. ['logBB', 'solubility']
# output format: return prepared array in format [id, predicted_value1, predicted_value2, ...]

# if ad is specified, return only compounds which are in ad
def prepare_array(input_pred, bounded_box, parameters_to_predict):
    # load and prepare data to numpy array
    working_list = []
    with open(input_pred, 'r') as in_f:
        for line in in_f:
            if '#' in line:
                pass
            else:
                line = [cell.strip() for cell in line.split('\t')]
                working_list.append([int(line[0])])
                working_list[-1].append(float(line[-2]))

                if bounded_box:
                    working_list[-1].append(line[-1] == 'True')
    working_arr = np.array(working_list)

    # split data after all predictions
    working_arr = np.split(working_arr, len(parameters_to_predict), axis=0)

    # add to final_array first column with ids
    final_ar = working_arr[0]

    # add columns to final_arr with predictions and bounded box if specified
    for ar in working_arr[1:]:
        if bounded_box:
            final_ar = np.hstack((final_ar, ar[:, [1, 2]]))
        else:
            final_ar = np.hstack((final_ar, ar[:, [1]]))

    # filter array if bound_box
    if bounded_box:
        for i in range(len(parameters_to_predict)):
            final_ar = final_ar[final_ar[:, 2 * i + 2] == 1]
        return final_ar[:, [0] + [i * 2 + 1 for i in range(len(parameters_to_predict))]]
    else:
        return final_ar


# convert input thresholds to parsed 2D list
# input format: thresholds: list of thresholds, e.g. ['more4', 'betwenn-0.5to1', less'-2']
# output format: return list of parsed threshold, e.g. [['more',4], ['between', -0.5, 1], ['less', -2]]
def parse_threshold_pareto_filtering(thresholds):
    threshold_match = []
    for threshold in thresholds:
        if 'less' in threshold:
            threshold_match.append(['less', float(threshold[4:])])
        elif 'more' in threshold:
            threshold_match.append(['more', float(threshold[4:])])
        elif 'between' in threshold:
            threshold_match.append(
                ['between', float(threshold[7:].split('to')[0]), float(threshold[7:].split('to')[1])])
    return threshold_match


def parse_threshold_desirability(thresholds):
    for threshold in thresholds:
        if 'desirability' in threshold:
            return threshold.split('#')[0].split('_')[1:], int(threshold.split('#')[1][0])


# --------------------------Filtering--------------------------
# filter rows according to threshold
# input format: working_ar: numpy array with format [id, predicted_value1, predicted_value2, ...],
#               thresholds: parsed thresholds with format e.g. [['more',4], ['between', -0.5, 1], ['less', -2]]
# output format: return filtered array in same format as working_ar
def filtering(working_ar, thresholds):
    for index_of_threshold, threshold in enumerate(thresholds):
        if threshold[0] == 'more':
            working_ar = working_ar[working_ar[:, index_of_threshold + 1] > threshold[1]]
        elif threshold[0] == 'less':
            working_ar = working_ar[working_ar[:, index_of_threshold + 1] < threshold[1]]
        else:
            working_ar = working_ar[(working_ar[:, index_of_threshold + 1] <= threshold[2])
                                    & (working_ar[:, index_of_threshold + 1] >= threshold[1])]
    return working_ar


# --------------------------Pareto--------------------------

# computes simple distance from predicted values, which are stored in 1D array
# input format: col: 1D numpy array with predicted values,
#                    threshold_type: string with type of threshold, e.g. 'more', 'less', 'between'
#                    threshold_value1: float value
#                    threshold_value2: if threshold_type == 'between': second value of threshold, else: None
# output format: return 1D numpy array with computed distances to coresponding threshold
#                if predicted value match the threshold, then return 0
def get_distance_from_threshold(col, threshold_type, threshold_value1, threshold_value2=None):
    if threshold_type == 'more':
        if col >= threshold_value1: return 0
        else: return threshold_value1 - col
    elif threshold_type == 'less':
        if col < threshold_value1: return 0
        return col - threshold_value1
    else:
        if col >= threshold_value1 and col <= threshold_value2:
            return 0
        else:
            return threshold_value1 - col if col < threshold_value1 else col - threshold_value2


# it finds compounds which lies on pareto frontier, also if some compounds match the threshold,
# it saves them into separate file called output_match_pareto.sdf
# input format: input_sdf: path to standardized and labeled sdf file which you used for prediction
#               working_ar: numpy array with format [id, predicted_value1, predicted_value2, ...],
#               thresholds: parsed thresholds with format e.g. [['more',4], ['between', -0.5, 1], ['less', -2]]
#               predictions: list of parameters which you want to predict, e.g. ['logBB', 'solubility']
# output format: return array of compounds which lies on pareto frontier in same format as working_ar
def pareto(input_sdf, working_ar, thresholds, parameters_to_predict):
    # vectorize calculation of distances
    get_distance = np.vectorize(get_distance_from_threshold, otypes=[np.float64])

    # output list of compounds which lies on pareto frontier
    output = []
    # output dictionary for compounds which match the thresholds
    match_threshold_output = {}

    # compute distances, then append columns with stored distances from thresholds
    # format [id, pred1, pred2, pred3, ... , dist1, dist2, dist3 ...]
    for index_of_threshold, threshold in enumerate(thresholds):
        # if the threshold is 'between'
        if len(threshold) == 3:
            working_ar = np.c_[working_ar, get_distance(working_ar[:, index_of_threshold + 1], threshold[0], threshold[1], threshold[2])]
        else:
            working_ar = np.c_[working_ar, get_distance(working_ar[:, index_of_threshold + 1], threshold[0], threshold[1])]

    # add to match_threshold_output rows which match the thresholds
    tmp = working_ar[np.logical_and.reduce([working_ar[:,i+1+len(thresholds)] == 0 for i in range(len(thresholds))])]
    for index in range(tmp.shape[0]):
        output.append([item for item in tmp[index][0:len(thresholds)+1]])
    # if some compounds match the threshold
    if len(output) != 0:
        match_threshold_output['pareto'] = output
        save_output(input_sdf, 'output_match.sdf', match_threshold_output, parameters_to_predict)
        # clear the output
        output = []
        # remove rows which are in output
        working_ar = working_ar[np.logical_or.reduce([working_ar[:,i+1+len(thresholds)] != 0 for i in range(len(thresholds))])]

    # prepare points to pareto function, format [[dist1, dist2, ...], [dist1, dist2, ...], ...]
    input_to_pareto_function = working_ar[:, 0 + len(thresholds) + 1]
    for index_of_threshold in range(len(thresholds) - 1):
        input_to_pareto_function = np.c_[input_to_pareto_function, working_ar[:, index_of_threshold + len(thresholds) + 2]]
    input_to_pareto_function = input_to_pareto_function.tolist()

    # get list of indexes from pareto frontier
    input_to_pareto_function = pareto_alg.simple_cull(input_to_pareto_function, pareto_alg.dominates_min)

    for index in input_to_pareto_function:
        output.append([item for item in working_ar[index][0:len(thresholds)+1]])
    return output


# --------------------------Desirability functions--------------------------
def process_function(in_function):
    functions = []
    for fnc in in_function:
        fnc = fnc.split(',')
        for index, fun in enumerate(fnc):
            fnc[index] = [fun.split(":")[0]]  + [parse_expr(fun.split(":")[1])]
        functions.append(fnc)
    return functions


def get_norm_value(function, x_input):
    x = symbols("x")
    for index, bound in enumerate(function):
        if round(x_input, 5) <= float(bound[0]): return float(function[index][1].subs(x, x_input))
    return 0


def desirability(input_sdf, working_ar, threshold_filtering, threshold_desire, number_of_compounds, parameters_to_predict):

    # filter compounds which match the threshold
    match_threshold_output = {}
    output = filtering(working_ar.copy(), threshold_filtering)
    match_threshold_output['match'] = output
    save_output(input_sdf, 'output.sdf', match_threshold_output, parameters_to_predict)

    # delete filtered compounds from working array
    for id in output[:, 1]:
        working_ar = working_ar[working_ar[:, 1] != id]

    working_ar = np.hstack((working_ar, np.zeros((working_ar.shape[0], 1))))
    functions = process_function(threshold_desire)

    if number_of_compounds > working_ar.shape[0]:
        number_of_compounds = working_ar.shape[0]

    number_of_parameters = len(functions)
    for index, row in enumerate(working_ar):
        row_average = 0
        for predictions, fnc in zip(row[1:-1], functions):
            row_average += get_norm_value(fnc, predictions)
        working_ar[index, -1] = row_average / number_of_parameters

    tmp_arr = working_ar[working_ar[:, -1].argsort()][::-1][:number_of_compounds, :-1]
    return tmp_arr[tmp_arr[:, 0].argsort()]

def main_params(input_sdf, input_pred, parameters_to_predict, output_file, methods, thresholds, use_bounded_box):

    # preparing threshold for filtering or pareto
    threshold_match = parse_threshold_pareto_filtering(thresholds.copy())

    #prepare array
    working_array = prepare_array(input_pred, use_bounded_box, parameters_to_predict)

    # output dictionary
    output = {}

    if 'filtering' in methods:
        output['process_predictions'] = filtering(working_array, threshold_match)
    elif 'pareto' in methods:
        output['process_predictions'] = pareto(input_sdf, working_array, threshold_match, parameters_to_predict)
    else: # desirability
        threshold_desire, number_of_compounds = parse_threshold_desirability(thresholds)
        output['process_predictions'] = desirability(input_sdf, working_array, threshold_match, threshold_desire, number_of_compounds, parameters_to_predict)

    # save outputs to coresponding files
    save_output(input_sdf, output_file, output, parameters_to_predict)


def main():
    parser = argparse.ArgumentParser(description=
                                     'Process predicted values and send compounds which match threshold to output pool.')
    parser.add_argument('-is', '--in_sdf', metavar='path_to_file_with_standardized_compounds', required=True,
                        help='path to file which contains standardized compounds')
    parser.add_argument('-ip', '--in_pred', metavar='path_to_file_with_predictions', required=True,
                        help='path to file which contains predicted values for properties')
    parser.add_argument('-o', '--out', metavar='output.txt', required=True,
                        help='processed predictions with compounds, which match the threshold')
    parser.add_argument('-p', '--properties', metavar='[LOGBB solubility]', required=True, nargs='*',
                        help='predicted properties')
    parser.add_argument('-m', '--methods', metavar='[filtering pareto]', required=True, nargs='*',
                        help='define which method will be used for getting output file')
    parser.add_argument('-t', '--thresholds', metavar='[>1.78 <12 -2.05to2.97]', required=True, nargs='*',
                        help='thresholds to be match, written in the same order as properties')
    parser.add_argument('-a', '--ad', action='store_true', default=False,
                        help='save to output file only if it is in the application domain')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == "in_sdf": input_sdf = v
        if o == "in_pred": input_pred = v
        if o == "out": output_file = v
        if o == "properties": predictions = v
        if o == "methods": methods = v
        if o == "thresholds": thresholds = v
        if o == "ad": use_bounded_box = v
        if o == "pareto": pareto = v

    main_params(input_sdf, input_pred, predictions, output_file, methods, thresholds, use_bounded_box)


if __name__ == '__main__':
    main()
