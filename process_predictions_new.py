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
    This function is applied on pandas DataFrame.

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


def main(in_sdf, in_pred, out_fname, parameters, optimization_methods,
         desirabilities, thresholds, ad):

    selected_compounds_index = set()

    # process all predictions
    predictions = prepare_working_arr(in_pred, parameters, ad)
    if predictions is None:
        print('Compounds are not in ad. Calculating outside ad!')
        predictions = prepare_working_arr(in_pred, parameters, False)

    # filtering
    thresholds = parse_threshold(thresholds)

    # compute distances from thresholds and use it with pareto if specified
    pareto_predictions = predictions.copy()
    for parameter, threshold in zip(parameters, thresholds):
        pareto_predictions[parameter] = predictions[parameter].apply(compute_distance_from_threshold,
                                                              threshold=threshold)

    # find compounds which are in threshold
    output_filtering = pareto_predictions[pareto_predictions.sum(axis=1) == 0]
    output_filtering = predictions.loc[output_filtering.index].copy()

    for method in optimization_methods:
        if method == 'pareto':

            # use compounds which are not in threshold
            pareto_predictions = pareto_predictions[pareto_predictions.sum(axis=1) > 0]

            input_to_pareto = prepare_points_for_pareto(pareto_predictions)

            # get list of indexes from pareto frontier
            pareto = pareto_alg.simple_cull(input_to_pareto, pareto_alg.dominates_min)

            for index in predictions.loc[pareto_predictions.iloc[pareto].index].index:
                selected_compounds_index.add(index)

        elif method == 'desirability':

            print('desirability')
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
    parser.add_argument('-d', '--desirabilities', nargs='*',
                        help='if desirability method specified, need to specify desirability string')
    parser.add_argument('-t', '--thresholds', required=True, nargs='*',
                        help='thresholds to be match, written in the same order as properties')
    parser.add_argument('-a', '--ad', action='store_true', default=False,
                        help='save to output file only if it is in the application domain')
    args = vars(parser.parse_args())

    main(args['in_sdf'], args['in_pred'], args['out'], args['parameters'],
         args['methods'], args['desirabilities'], args['thresholds'], args['ad'])
