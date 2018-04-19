#!/usr/bin/env python

import argparse
import math

from rdkit import Chem
import pandas as pd

from optimizer_utils import parse_threshold


def main(in_sdf_f, in_contrib_f, out_frag_f, out_worst_f, parameters, ranges,
         types_of_alg, thresholds, n_worst):

    print(in_sdf_f, in_contrib_f, out_frag_f, out_worst_f, parameters, ranges,
          types_of_alg, thresholds, n_worst)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description=
                                     'Normalize fragments contributions and pick the worst ones.')
    parser.add_argument('-is', '--in_sdf', required=True,
                        help='path to file which contains selected compounds from pareto/desirability/...')
    parser.add_argument('-in', '--in_contrib', required=True, nargs='*',
                        help='list of paths to calculated contributions of fragments for specific parameter')
    parser.add_argument('-of', '--out_frag', required=True,
                        help='path to file where you want to store normalized contributions for fragments')
    parser.add_argument('-ow', '--out_worst', required=True,
                        help='path to file where you want to store the worst fragments')
    parser.add_argument('-p', '--parameters', required=True, nargs='*',
                        help='list of predicted parameters')
    parser.add_argument('-r', '--ranges', required=True, nargs='*',
                        help='list of ranges of all parameters')
    parser.add_argument('-m', '--models', required=True, nargs='*',
                        help='types of algorithms used for predictions, e.g. [svm_gbm rf rf_gbm_svm]')
    parser.add_argument('-t', '--thresholds', required=True, nargs='*',
                        help='thresholds to be matched')
    parser.add_argument('-n', '--n_worst', action='store', type=int,
                        help='specifies number of worst fragments')

    args = vars(parser.parse_args())

    main(args['in_sdf'], args['in_contrib'], args['out_frag'], args['out_worst'],
         args['parameters'], args['ranges'], args['models'], args['thresholds'],
         args['n_worst'])
