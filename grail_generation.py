"""
@author: Arthur Garon, Thomas Seidel
Department of Pharmaceutical Sciences
University of Vienna
"""

import CDPL.Chem as Chem
import CDPL.Math as Math

import argparse
import os
from multiprocessing.dummy import Pool as ThreadPool

from common import *


def parseArguments():
    parser = argparse.ArgumentParser(description='>>> Generates a grid.cdf file (CDPKIT format) from both an aligned cdf and a text file listing the residues to cover by in the grids.')
    parser.add_argument('-cdf',
                        dest='cdf',
                        required=True,
                        help='[Required] The path of the cdf file',
                        nargs=1)
    parser.add_argument('-res',
                        dest='residue',
                        help='[Semi-optional] The path of the text file listing the residues of interest (if not present, -grid must be given!)',
                        nargs=1,
                        default=None)
    parser.add_argument('-grid',
                        dest='grid_spec',
                        help='[Semi-optional] The specification of the grid bounding box (if not present, -res must be given!)',
                        nargs=1,
                        default=None)
    parser.add_argument('-t',
                        dest='thread',
                        help='[Optional] The number of threads to launch (Default: no multithreading)',
                        nargs=1,
                        default=None)
    parser.add_argument('-o',
                        dest='output',
                        help='[Optional] The output folder where the grid.cdf file will be generated (Default: current directory)',
                        nargs=1,
                        default=None)
    parser.add_argument('-lig',
                        dest='ligand_code',
                        help='[Optional] The 3-letter code of the ligand to generate an atom density grid',
                        nargs=1,
                        default=None)
    parser.add_argument('-bur',
                        dest='bur',
                        help='[Optional] To generate a buriedness grid (Default: False)',
                        action='store_true',
                        default=False)
     
    parse_args = parser.parse_args()
    return parse_args

def cdfSplit(cdf_file ,chunk_number, output):
    main_cdf = loadCDFMolecule(cdf_file)
    length = len(Chem.get3DCoordinatesArray(main_cdf.atoms[1]))
    del main_cdf

    chunks_start = [0]
    chunks_end = []
    chunk_size = int(length / chunk_number) + (length % chunk_number > 0)

    for i in range(length):
        if i % chunk_size == 0 and i != 0 and i != (length - 1):
            chunks_end.append(i)
            chunks_start.append(i)
        elif i == (length - 1):
            chunks_end.append(i)

    for chunk_index in range(int(chunk_number)):
        print('> Cdf chunk {} generation ...'.format(chunk_index))
        tmp_cdf = loadCDFMolecule(cdf_file)
        for atom_index in range(tmp_cdf.numAtoms):
            tmp_coordinate_array = Math.Vector3DArray()
            for coord in [Chem.get3DCoordinatesArray(tmp_cdf.getAtom(atom_index))[i] for i in range(chunks_start[chunk_index], chunks_end[chunk_index])]:
                tmp_coordinate_array.addElement(coord)

            Chem.set3DCoordinatesArray(tmp_cdf.getAtom(atom_index), tmp_coordinate_array)

        try:
            Chem.FileCDFMolecularGraphWriter(output + os.path.basename(cdf_file)[:-4] + '_chunk_'+ str(chunk_index) +'.cdf').write(tmp_cdf)
        except:
            print('> Cdf_mol writing failure.')
            raise

def functMT(list_arg):
    grailGeneration(list_arg[0], list_arg[1], list_arg[2], list_arg[3], list_arg[4], list_arg[5])

def cat(infilenames, output):
    outfilename = output + '/' + os.path.basename(infilenames[0]).split('_chunk_')[0] + '_grid.cdf'
    with open(outfilename, 'w') as outfile:
        for infilename in infilenames:
            with open(infilename) as infile:
                outfile.write(infile.read())

            os.remove(infilename)

if __name__ == '__main__':
    args = parseArguments()

    if args.residue is None and args.grid_spec is None:
        print('Error: either -res or -grid option must be given')
        sys.exit(1)

    cdf = args.cdf[0]
    bur = args.bur

    if args.thread is None:
        thread = 0
    else:
        thread = int(args.thread[0])

    if args.output is None:
        output = './'
    else:
        output = args.output[0]

    if output[-1] != '/':
        output += '/'

    if args.ligand_code is None:
        lig = ' '
    else:
        lig = args.ligand_code[0]

    if args.grid_spec is None:
        bbox_min, bbox_max = getGridInfo(cdf, args.residue[0])
    else:
        bbox_min, bbox_max = convertToBBoxMinMax(args.grid_spec[0])

    if thread == 0:
        grailGeneration(cdf, lig, bbox_min, bbox_max, output, bur)
    else:
        cdfSplit(cdf, thread, output)

        arg_mt = [(output + os.path.basename(cdf)[:-4] + '_chunk_'+ str(i) +'.cdf', lig, bbox_min, bbox_max, output, bur) for i in range(thread)]

        pool = ThreadPool(thread)
        pool.map(functMT, arg_mt)
        pool.close()
        pool.join()

        tmp_grid = [output + '/' + os.path.basename(cdf)[:-4] + '_chunk_'+ str(i) +'_grid.cdf' for i in range(thread)]
        cat(tmp_grid, output)
        for x in [output + '/' + os.path.basename(cdf)[:-4] + '_chunk_'+ str(i) +'.cdf' for i in range(thread)]:
            os.remove(x)

    print('> Done')






