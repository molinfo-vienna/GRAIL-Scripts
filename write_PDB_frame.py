"""
@author: Arthur Garon
department of pharmaceutical chemistry
university of vienna
"""

import argparse
import os
import MDAnalysis
import warnings

warnings.filterwarnings("ignore") ### Avoid warning from reading gromacs files


def parseArguments():
    parser = argparse.ArgumentParser(description='>>> Generates a pdb file from a cdf file.')
    parser.add_argument('-trj',
                        dest='trajectory',
                        required=True,
                        help='[Required] The path of the trajectory file',
                        nargs=1)
    parser.add_argument('-top',
                        dest='topology',
                        required=True,
                        help='[Required] The path of the topology file',
                        nargs=1)

    parser.add_argument('-o',
                        dest='output',
                        help='[Optional] The output folder where the pdb file will be generated (Default: current directory)',
                        nargs=1,
                        default=None)
    parser.add_argument('-fi',
                        dest='frame_index',
                        help='[Optional] Frame index to consider (Default: 0)',
                        nargs=1,
                        default=None)

    parse_args = parser.parse_args()
    return parse_args

def writePDB(trajectory, topology, frame_index, output):

    u = MDAnalysis.Universe(topology, trajectory)
    structure = u.select_atoms("all")

    print('> Writing frame {} as PDB file ...'.format(frame_index))
    for step in u.trajectory[frame_index:frame_index + 1]:
        out_path = output + os.path.basename(trajectory)[:-4] + '_frame_{}.pdb'.format(frame_index)
        structure.write(out_path)

if __name__ == '__main__':
    args = parseArguments()

    trajectory = args.trajectory[0]
    topology = args.topology[0]

    if args.output is None:
        output = './'
    else:
        output = args.output[0]

    if args.frame_index is None:
        frame_index = 0
    else:
        frame_index = int(args.frame_index[0])

    writePDB(trajectory, topology, frame_index, output)
