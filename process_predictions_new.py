#!/usr/bin/env python

import argparse

import numpy as np
import pandas as pd

import pareto_simple_cull as pareto_alg

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


def main(in_sdf, in_pred, out, parameters, method, desirability, threshold, ad):
    predictions = prepare_working_arr(in_pred, parameters, ad)
    if predictions is None:
        print('Compounds are not in ad. Calculating outside ad!')
        predictions = prepare_working_arr(in_pred, parameters, False)
    print(predictions)


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
    parser.add_argument('-m', '--method', required=True,
                        help='method for processing predictions: pareto or desirability')
    parser.add_argument('-d', '--desirability', nargs='*',
                        help='if desirability method specified, need to specify desirability string')
    parser.add_argument('-t', '--threshold', required=True, nargs='*',
                        help='thresholds to be match, written in the same order as properties')
    parser.add_argument('-a', '--ad', action='store_true', default=False,
                        help='save to output file only if it is in the application domain')
    args = vars(parser.parse_args())

    main(args['in_sdf'], args['in_pred'], args['out'], args['parameters'],
         args['method'], args['desirability'], args['threshold'], args['ad'])
