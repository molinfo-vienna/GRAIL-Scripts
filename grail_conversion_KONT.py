"""
@author: Arthur Garon
department of pharmaceutical chemistry
university of vienna
"""

import sys
import os
import CDPL.Grid as Grid
import CDPL.Math as Math
import argparse

from common import *


def parseArguments():
    parser = argparse.ArgumentParser(description='>>> Generate kont files (GRID program format) from a grid.cdf file.')
    parser.add_argument('-gcdf',
                        dest='gcdf',
                        required=True,
                        help='[Required] The path of the grid.cdf file',
                        nargs=1)

    parser.add_argument('-o',
                        dest='output',
                        help='[Optional] The output folder where the kont files will be generated (Default: current directory)',
                        nargs=1,
                        default=None)
    parser.add_argument('-fi',
                        dest='frame_index',
                        help='[Optional] Index of the frame to consider (Default: 0)',
                        nargs=1,
                        default=None)

    parse_args = parser.parse_args()
    return parse_args

if __name__ == '__main__':
    args = parseArguments()

    gcdf = args.gcdf[0]

    if args.output is None:
        output = './'
    else:
        output = args.output[0]

    if output[-1] != '/':
        output += '/'

    if args.frame_index is None:
        frame_index = 0
    else:
        frame_index = int(args.frame_index[0])

    kont_conversion(gcdf, frame_index, output)

