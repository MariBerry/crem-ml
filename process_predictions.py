#!/usr/bin/env python

import argparse
import os

import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None

import pareto_simple_cull as pareto_alg
from optimizer_utils import save_output_poll

from sympy import symbols
from sympy.parsing.sympy_parser import parse_expr

from typing import List
from typing import NewType
pandas_table = NewType('Processed pandas table with id of compound and predicted properties',
                      pd.DataFrame
                      )


def prepare_working_arr(in_pred: List, parameters: List, bounded_box: bool) -> pandas_table:
    """
    Reads file with predictions and process it into pandas table

    :param in_pred: list of paths to files with all predictions
    :param parameters: list of parameters to predict
    :param bounded_box: if True, then return only compounds within bounded box
    :return: pandas table with processed predictions
    """

    tables = [pd.read_table(file) for file in in_pred]

    for table, parameter in zip(tables, parameters):

        # if ad
        if bounded_box:
            if table[table.bound_box == 1].shape[0] == 0: # no compounds in ad
                return None
            else:
                table = table[table.bound_box == 1]

        table.drop('bound_box', axis=1, inplace=True)
        cols_to_drop = list(range(1,table.shape[1]-1))
        table.drop(table.columns[cols_to_drop], axis=1, inplace=True)
        table.rename(columns={table.columns[0]: 'id', 'consensus': parameter}, inplace=True)
        table.set_index('id', inplace=True)

    tables = pd.concat(tables, axis=1, join='inner')

    # if ad
    if bounded_box:
        if tables.shape[0] == 0:
            return None

    return tables

def parse_threshold(thresholds: List) -> List:
    """
    Convert input thresholds to parsed 2D list

    :param thresholds: list of thresholds, e.g. ['more4', 'betwenn-0.5to1', less'-2']
    :return: list of parsed threshold, e.g. [['more',4], ['between', -0.5, 1], ['less', -2]]
    """

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

def compute_distance_from_threshold(x: float, threshold: List) -> float:
    """
    Compute distance of predicted value from threshold value.
    This function is applied on pandas series.

    :param x: predicted value
    :param threshold: parsed threshold list, e.g. ['more', 7]
    :return: computed distance from one predicted value
    """

    if threshold[0] == 'more':
        if x > threshold[1]:
            return 0
        else:
            return threshold[1] - x

    elif threshold[0] == 'less':
        if x < threshold[1]:
            return 0
        else:
            return x - threshold[1]

    else:   # between
        if x >= threshold[1] and x <= threshold[2]:
            return 0
        else:
            return threshold[1] - x if x < threshold[1] else x - threshold[2]

def prepare_points_for_pareto(table: pandas_table) -> List:
    """
    Process pandas table into list of points which are distances from threshold.

    :param table: pandas table with distances from threshold of compounds
    :return: list of lists with distances, e.g. [[dist11, dist12], [dist21, dist22]]
    """

    points = []
    for index, row in table.iterrows():
        points.append(row.tolist())
    return points

def process_desirability_functions(desirabilities: List) -> List:
    """
    It parses list of desirabilities to special format

    :param desirabilities: list of desirabilities
    :return: special list of functions, e.g. [['0.45', 0], ['0.55', 10*x - 4.5], ['1000', 1]] - one desirability function
    """

    if len(desirabilities) == 1:   # if we have only one function
        function = desirabilities[0].split(",")
        for index, fun in enumerate(function):
            function[index] = [fun.split(":")[0]]  + [parse_expr(fun.split(":")[1])]
        return function
    else:
        functions = []
        for fnc in desirabilities:
            fnc = fnc.split(',')
            for index, fun in enumerate(fnc):
                fnc[index] = [fun.split(":")[0]]  + [parse_expr(fun.split(":")[1])]
            functions.append(fnc)
        return functions

def get_norm_value(x_input: float, function: List) -> float:
    """
    Compute norm value from predicted parameter. This function is applied on
    pandas series

    :param x: predicted value
    :param function: desirability function, e.g. [['0.45', 0], ['0.55', 10*x - 4.5], ['1000', 1]]
    :return: norm value
    """

    x = symbols("x")
    for index, bound in enumerate(function):
        if round(x_input, 5) <= float(bound[0]): return float(function[index][1].subs(x, x_input))
    return 0


def main(in_sdf, in_pred, out_fname, parameters, optimization_methods,
         thresholds, ad, desirabilities=[], n_compounds=0):

    selected_compounds_index = set()

    # process all predictions
    predictions = prepare_working_arr(in_pred, parameters, ad)
    if predictions is None:
        print('Compounds are not in ad. Calculating outside ad!')
        predictions = prepare_working_arr(in_pred, parameters, False)

    # filtering
    thresholds = parse_threshold(thresholds)

    # compute distances from thresholds and use it with pareto if specified
    distance_predictions = predictions.copy()
    for parameter, threshold in zip(parameters, thresholds):
        distance_predictions[parameter] = predictions[parameter].apply(compute_distance_from_threshold,
                                                              threshold=threshold)

    # find compounds which are in threshold
    output_filtering = distance_predictions[distance_predictions.sum(axis=1) == 0]
    output_filtering = predictions.loc[output_filtering.index].copy()

    for method in optimization_methods:
        if method == 'pareto':

            # use compounds which are not in threshold
            distance_predictions = distance_predictions[distance_predictions.sum(axis=1) > 0]

            input_to_pareto = prepare_points_for_pareto(distance_predictions)

            # get list of indexes from pareto frontier
            pareto = pareto_alg.simple_cull(input_to_pareto, pareto_alg.dominates_min)

            for index in predictions.loc[distance_predictions.iloc[pareto].index].index:
                selected_compounds_index.add(index)

        elif method == 'desirability':

            desirability_predictions = predictions.copy()

            # use compounds which are not in threshold
            desirability_predictions = desirability_predictions.loc[distance_predictions.sum(axis=1) > 0]

            # we have less or equal compounds in input sdf then we specified
            # that we need from this stage, so we use all of them
            if n_compounds >= desirability_predictions.shape[0]:
                for index in predictions.loc[desirability_predictions.index].index:
                    selected_compounds_index.add(index)
            else:
                functions = process_desirability_functions(desirabilities)

                for parameter, function in zip(parameters, functions):
                    desirability_predictions[parameter] = desirability_predictions[parameter].\
                        apply(get_norm_value, function=function)

                desirability_predictions['desirability'] = desirability_predictions.sum(axis=1)/(len(parameters))
                desirability_predictions = desirability_predictions.sort_values(by='desirability', ascending=False)

                for index in predictions.loc[desirability_predictions.head(n_compounds).index].index:
                    selected_compounds_index.add(index)

        else:
            print('Unspecified optimization method!')

    if output_filtering.shape[0] > 0:
        save_output_poll(in_sdf,
                         os.path.join(os.path.dirname(in_sdf), 'output_match.sdf'),
                         output_filtering)

    # save selected compounds
    save_output_poll(in_sdf,
                     out_fname,
                     predictions.loc[list(selected_compounds_index)])

    return output_filtering.shape[0]

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Process predicted values with specified optimization method')
    parser.add_argument('-is', '--in_sdf', required=True,
                        help='path to file which contains standardized compounds')
    parser.add_argument('-ip', '--in_pred', required=True, nargs='*',
                        help='path to files which contains predicted values for properties')
    parser.add_argument('-o', '--out', required=True,
                        help='processed predictions using specified opt. method')
    parser.add_argument('-p', '--parameters', required=True, nargs='*',
                        help='parameters for prediction')
    parser.add_argument('-m', '--methods', required=True, nargs='*',
                        help='method for processing predictions: pareto or desirability')
    parser.add_argument('-t', '--thresholds', required=True, nargs='*',
                        help='thresholds to be match, written in the same order as properties')
    parser.add_argument('-a', '--ad', action='store_true', default=False,
                        help='save to output file only if it is in the application domain')
    parser.add_argument('-d', '--desirabilities', nargs='*',
                        help='if desirability method specified, need to specify desirability string')
    parser.add_argument('-n', '--n_compounds', action='store', type=int,
                        help='if desirability method specified, need to specify number of selected compounds')
    args = vars(parser.parse_args())

    main(args['in_sdf'], args['in_pred'], args['out'], args['parameters'],
         args['methods'], args['thresholds'], args['ad'],
         args['desirabilities'], args['n_compounds'])
