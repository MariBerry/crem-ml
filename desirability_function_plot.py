from sympy import *
from sympy.parsing.sympy_parser import parse_expr
import numpy as np
import matplotlib.pyplot as plt

import argparse

from process_predictions import get_norm_value, process_function

def main_params(in_function, in_range):

    fnc = process_function(in_function)
    print(fnc)
    in_values = np.arange(int(in_range[0]), int(in_range[1]), 0.01)
    out_values = in_values.copy()

    for index, value in enumerate(out_values):
        out_values[index] = get_norm_value(fnc, value)
    plt.plot(in_values, out_values)
    plt.grid()
    plt.show()


def main():
    parser = argparse.ArgumentParser(description=
                'Plots your function from string for desirability functions')
    parser.add_argument('-if', '--in_function', metavar='0.45:0,0.55:10*x-4.5,1000:1',
                        required=True, help='input string for plot')
    parser.add_argument('-ir', '--in_range', metavar='[-5 10]', nargs='*',
                        required=True, help='type two number which represents range your plot')

    args = vars(parser.parse_args())
    for o, v in args.items():
        if o == 'in_function': in_function = v
        if o == 'in_range': in_range = v

    main_params(in_function, in_range)

if __name__ == '__main__':
    main()
