#!/usr/bin/env python

import argparse

import numpy as np
from rdkit import Chem

import pareto_simple_cull as pareto_alg


# --------------------------Saving outputs--------------------------
# save output dictionary to file with predicted values
# input format: input_sdf: path to standardized and labeled sdf file which you used for prediction
#               output_file: name of output file, it will be processed according to predictions,
#                            e.g. output.sdf -> output_filtering.sdf/output_pareto.sdf
#               input_dict: output dictionary where are selected compounds stored,
#                            e.g. {'filtering': [[id1, predicted_value1, predicted_value2],[...]], 'pareto':[[...]]}
#               predictions: list of parameters name, e.g. ['logBB', 'solubility', ...]
def save_output(input_sdf, output_file, input_dict, predictions):
    # store all compounds from sdf file
    compounds = Chem.SDMolSupplier(input_sdf, removeHs=False, sanitize=False)

    # prepare list of files
    list_of_files = []
    for parameter in input_dict.keys():
        list_of_files.append(open(output_file.split('.')[0] + '_' + parameter + '.' + output_file.split('.')[1], 'a'))
        list_of_files[-1] = Chem.SDWriter(list_of_files[-1])

    # iterate through all compounds
    for mol in compounds:
        # iterate through all keys in dictionary
        for index_of_file, output_list in enumerate(input_dict.values()):
            for item in output_list:
                if int(mol.GetProp('ID')) == int(item[0]):
                    for index_of_prediction, predict in enumerate(predictions):
                        mol.SetDoubleProp(predict, item[index_of_prediction + 1])
                    list_of_files[index_of_file].write(mol)
                    break
                else:
                    if int(mol.GetProp('ID')) < int(item[0]):
                        break

    for f in list_of_files: f.close()


# --------------------------Preparing variables--------------------------
# prepares numpy array of predicted values for compounds

# input format: input_pred: path to file with predictions,
#               bounded_box: True/False,
#               predictions: list of parameters which you want to predict, e.g. ['logBB', 'solubility']
# output format: return prepared array in format [id, predicted_value1, predicted_value2, ...]

# if ad is specified, return only compounds which are in ad
def prepare_array(input_pred, bounded_box, predictions):
    in_data = np.genfromtxt(input_pred, dtype=None, delimiter='\t')

    working_arr = np.asarray([[in_data[0][0], in_data[0][-2], in_data[0][-1]]])
    for item in in_data[1:]:
        working_arr = np.append(working_arr, [[item[0], item[-2], item[-1]]], axis=0)

    # -2 prediction column, -1 bounded box column
    if bounded_box:
        working_arr = working_arr[:, [0, -2, -1]]
    else:
        working_arr = working_arr[:, [0, -2]]

    # split data after all predictions
    working_arr = np.split(working_arr, len(predictions), axis=0)

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
        for i in range(len(predictions)):
            final_ar = final_ar[final_ar[:, 2 * i + 2] == 1]
        return final_ar[:, [0] + [i * 2 + 1 for i in range(len(predictions))]]
    else:
        return final_ar


# convert input thresholds to parsed 2D list
# input format: thresholds: list of thresholds, e.g. ['more4', 'betwenn-0.5to1', less'-2']
# output format: return list of parsed threshold, e.g. [['more',4], ['between', -0.5, 1], ['less', -2]]
def parse_threshold(thresholds):
    threshold_match = []
    for threshold in thresholds:
        if 'less' in threshold:
            threshold_match.append(['less', float(threshold[4:])])
        elif 'more' in threshold:
            threshold_match.append(['more', float(threshold[4:])])
        else:
            threshold_match.append(
                ['between', float(threshold[7:].split('to')[0]), float(threshold[7:].split('to')[1])])
    return threshold_match


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
def pareto(input_sdf, working_ar, thresholds, predictions):
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
        save_output(input_sdf, 'output_match.sdf', match_threshold_output, predictions)
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


def main_params(input_sdf, input_pred, predictions, output_file, methods, thresholds, use_bounded_box):
    # preparing variables

    # parse threshold
    threshold_match = parse_threshold(thresholds)

    # output dictionary
    output = {}
    if 'filtering' in methods:
        working_array = prepare_array(input_pred, use_bounded_box, predictions)
        output['filtering'] = filtering(working_array, threshold_match)
    if 'pareto' in methods:
        working_array = prepare_array(input_pred, use_bounded_box, predictions)
        output['pareto'] = pareto(input_sdf, working_array, threshold_match, predictions)

    # save outputs to coresponding files
    save_output(input_sdf, output_file, output, predictions)


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
