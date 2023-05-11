"""
@author: Arthur Garon, Thomas Seidel
Department of Pharmaceutical Sciences
University of Vienna
"""

import math
import time

import CDPL.Biomol as Biomol
import CDPL.Chem as Chem
import CDPL.Math as Math
import argparse
import os

from common import *


def parseArguments():
    parser = argparse.ArgumentParser(description='>>> Generate a text file displaying the maximum fluctuation and the rmsf of all the residues during the simulation.')
    parser.add_argument('-cdf',
                        dest='cdf',
                        required=True,
                        help='[Required] The path of the cdf file',
                        nargs=1)

    parser.add_argument('-o',
                        dest='output',
                        help='[Optional] The output folder where the text file will be generated (Default: current directory)',
                        nargs=1,
                        default=None)

    parse_args = parser.parse_args()
    return parse_args

def calcMeanPos(positions):
    mean_pos = Math.Vector3D()

    for pos in positions:
        mean_pos += pos

    mean_pos /= len(positions)
    return mean_pos


def calcMaxFluct(positions, mean_pos):
    max_fluct = 0.0
    tmp = Math.Vector3D()

    for pos in positions:
        tmp.assign(pos)
        tmp -= mean_pos
        dev = Math.norm2(tmp)

        if dev > max_fluct:
            max_fluct = dev

    return max_fluct

def calcRMSF(positions, mean_pos):
    msf = 0.0
    tmp = Math.Vector3D()

    for pos in positions:
        tmp.assign(pos)
        tmp -= mean_pos
        msf += Math.innerProd(tmp, tmp)

    msf = msf / len(positions)
    return math.sqrt(msf)


def residueStability(cdf, output, name):
    initial_time = time.time()

    cdf_mol = loadCDFMolecule(cdf)
    residues = Biomol.ResidueList(cdf_mol)
    num_confs = Chem.getNumConformations(cdf_mol)

    res_dict = {}
    for res in residues:
        if Biomol.getResidueCode(res) != 'HOH':
            atoms = [atom for atom in res.atoms if Chem.getType(atom) != Chem.AtomType.H]
            positions = [calcAtomSetCentroid(atoms, i) for i in xrange(1,num_confs)]

            mean_pos = calcMeanPos(positions)
            fluct = calcMaxFluct(positions, mean_pos)
            rmsf = calcRMSF(positions, mean_pos)
            res_id = getResidueID(res)

            res_dict[res_id] = [fluct,rmsf]

    tmp_output = output + name + "_residue_info.txt"
    with open(tmp_output, 'w') as txt_writer:
        txt_writer.write('residue (name_resid_chain), max fluctuation, rmsf\n')
        for key, val in sorted(res_dict.items(), key=lambda i: (i[1],i[0])):
            txt_writer.write('{}, {:>15}, {:>15}\n'.format(key, val[0], val[1]))

    calc_time = time.time() - initial_time
    print('> Residue information file generated in {}s'.format(int(calc_time)))


if __name__ == '__main__':
    args = parseArguments()

    cdf = args.cdf[0]

    if args.output is None:
        output = './'
    else:
        output = args.output[0] + '/'

    residueStability(cdf, output, os.path.basename(cdf)[:-4])
