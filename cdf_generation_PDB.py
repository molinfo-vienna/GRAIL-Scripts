"""
@author: Arthur Garon, Thomas Seidel
Department of Pharmaceutical Sciences
University of Vienna
"""

import argparse
import os

from common import *


def parseArguments():
    parser = argparse.ArgumentParser(description='>>> Generate a cdf file (CDPKIT format) and a text file listing all residues from a PDB file.')
    parser.add_argument('-pdb',
                        dest='pdb',
                        required=True,
                        help='[Required] The path of the PDB file',
                        nargs=1)

    parser.add_argument('-n',
                        dest='name',
                        help='[Optional] The name of the cdf file (Default: name of the pdb file)',
                        nargs=1,
                        default=None)
    parser.add_argument('-o',
                        dest='output',
                        help='[Optional] The output folder where the cdf file will be generated (Default: current directory)',
                        nargs=1,
                        default=None)

    parse_args = parser.parse_args()
    return parse_args

if __name__ == '__main__':
    args = parseArguments()

    pdb = args.pdb[0]

    if args.name is None:
        name = os.path.basename(pdb)[:-4]
    else:
        name = args.name[0]

    if args.output is None:
        output = './'
    else:
        output = args.output[0]

    if output[-1] != '/':
        output += '/'

    cdfMol_pdb(pdb, output, name)



