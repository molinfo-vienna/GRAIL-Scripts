"""
@author: Arthur Garon
department of pharmaceutical chemistry
university of vienna
"""

import CDPL.Biomol as Biomol
import CDPL.Chem as Chem
import CDPL.Math as Math
import argparse
import os

from common import *


def parseArguments():
    parser = argparse.ArgumentParser(description='>>> Generates an aligned cdf file (CDPKIT format) from both a cdf and a text file with the amino acids to use as anchor.')
    parser.add_argument('-cdf',
                        dest='cdf',
                        required=True,
                        help='[Required] The path of the cdf file',
                        nargs=1)
    parser.add_argument('-res',
                        dest='residue',
                        required=True,
                        help='[Required] The path of the text file with the residues information',
                        nargs=1)

    parser.add_argument('-o',
                        dest='output',
                        help='[Optional] The output folder where the aligned cdf file will be generated (Default: current directory)',
                        nargs=1,
                        default=None)

    parse_args = parser.parse_args()
    return parse_args


def alignResidues(cdf, res_file, output, name):
    mol = loadCDFMolecule(cdf)
    residue_subset = Biomol.ResidueList(mol)

    with open(res_file) as f:
        res_subset_ids = sorted([line.strip(' ').rstrip('\n') for line in f])

    filterResidues(residue_subset, res_subset_ids)
    num_confs = Chem.getNumConformations(mol)

    print('> Aligning frames...')

    res_positions = []
    for res in residue_subset:
        atoms = [atom for atom in res.atoms if Biomol.isPDBBackboneAtom(atom)]
        if atoms == []:
            print('> Skip amino acid, invalid name: {}'.format(getResidueID(res)))
            res_subset_ids.remove(getResidueID(res))
        else:
            positions = []
            i = 0

            while i < num_confs:
                positions.append(calcAtomSetCentroid(atoms, i))
                i += 1

            res_positions.append(positions)

    filterResidues(residue_subset, res_subset_ids)
    alignment = Math.DKabschAlgorithm()

    al_ref_positions = Math.DMatrix(3, residue_subset.getSize())
    al_positions = Math.DMatrix(3, residue_subset.getSize())
    i = 0

    while i < residue_subset.getSize():
        pos = res_positions[i][0]

        al_ref_positions.setElement(0, i, pos[0])
        al_ref_positions.setElement(1, i, pos[1])
        al_ref_positions.setElement(2, i, pos[2])
        i += 1

    i = 1
    xform = Math.Matrix4D()

    while i < num_confs:
        j = 0

        while j < residue_subset.getSize():
            pos = res_positions[j][i]

            al_positions.setElement(0, j, pos[0])
            al_positions.setElement(1, j, pos[1])
            al_positions.setElement(2, j, pos[2])
            j += 1

        if not alignment.align(al_positions, al_ref_positions):
            print('> Could not align frame {} ...'.format(i))
        else:
            print('> Aligning frame {}'.format(i))
            xform.assign(alignment.getTransform())
            Chem.transformConformation(mol, i, xform)


        i += 1

    tmp_output = output + name + "_align.cdf"
    try:
        Chem.FileCDFMolecularGraphWriter(tmp_output).write(mol)
    except:
        print('> Cdf_mol writing failure.')
        raise


if __name__ == '__main__':
    args = parseArguments()

    cdf = args.cdf[0]
    residue = args.residue[0]

    if args.output is None:
        output = './'
    else:
        output = args.output[0] + '/'

    name = os.path.basename(cdf)[:-4]
    alignResidues(cdf, residue, output, name)


