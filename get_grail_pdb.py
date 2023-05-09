"""
@author: Arthur Garon
Department of Pharmaceutical Chemistry
University of Vienna
"""

import argparse
import os

from common import *


def parseArguments():
    parser = argparse.ArgumentParser(description='>>> Generate the pharmacophore interaction grids and protein atom density grids from a PDB file.')
    parser.add_argument('-pdb',
                        dest='pdb',
                        required=True,
                        help='[Required] The path of the PDB file',
                        nargs=1)
    parser.add_argument('-o',
                        dest='output',
                        help='[Optional] The output folder where the files will be generated (Default: current directory)',
                        nargs=1,
                        default=None)
    parser.add_argument('-lig',
                        dest='ligand_code',
                        help='[Optional] The 3-letters code of the ligand to generate an atom density grid',
                        nargs=1,
                        default=None)
    parser.add_argument('-grid',
                        dest='grid_spec',
                        help='[Optional] The specification of the grid bounding box',
                        nargs=1,
                        default=None)
    parser.add_argument('-bur',
                        dest='bur',
                        help='[Optional] To generate a buriedness grid (Default: False)',
                        action='store_true',
                        default=False)

    parse_args = parser.parse_args()
    return parse_args

if __name__ == '__main__':
    args = parseArguments()

    pdb = args.pdb[0]
    bur = args.bur

    if args.ligand_code is None:
        lig = ''
    else:
        lig = args.ligand_code[0]

    if args.output is None:
        output = './'
    else:
        output = args.output[0]

    if output[-1] != '/':
        output += '/'

    name = os.path.basename(pdb)[:-4]

    cdfMol_pdb(pdb, output, name, False)

    if args.grid_spec is None:
        bbox_min, bbox_max = getGridInfo(output + name + ".cdf", None)
    else:
        bbox_min, bbox_max = convertToBBoxMinMax(args.grid_spec[0])
            
    grailGeneration(output + name + ".cdf", lig, bbox_min, bbox_max, output, bur)
    kont_conversion(output + name + '_grid.cdf', 0, output)
