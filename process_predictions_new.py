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


def prepare_working_arr(in_pred: List, bounded_box: bool) -> pandas_table:
    """
    Reads file with predictions and process it into pandas table

    :param in_pred: list of paths to files with all predictions
    :param bounded_box: if True, then return only compounds within bounded box
    :return: pandas table with processed predictions
    """

    pass


def main(in_sdf, in_pred, out, parameter, method, desirability, threshold, ad):
    print(in_sdf, in_pred, out, parameter, method, desirability, threshold, ad)
    working_ar = prepare_working_arr(in_pred, ad)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Process predicted values with specified optimization method')
    parser.add_argument('-is', '--in_sdf', required=True,
                        help='path to file which contains standardized compounds')
    parser.add_argument('-ip', '--in_pred', required=True, nargs='*',
                        help='path to files which contains predicted values for properties')
    parser.add_argument('-o', '--out', required=True,
                        help='processed predictions using specified opt. method')
    parser.add_argument('-p', '--parameter', required=True, nargs='*',
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

    main(args['in_sdf'], args['in_pred'], args['out'], args['parameter'],
         args['method'], args['desirability'], args['threshold'], args['ad'])
