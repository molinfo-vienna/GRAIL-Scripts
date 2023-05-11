"""
@author: Arthur Garon, Thomas Seidel
Department of Pharmaceutical Sciences
University of Vienna
"""

import sys
import time

import CDPL.Chem as Chem
import CDPL.Grid as Grid
import CDPL.Math as Math
import numpy as np
import os
import argparse

def parseArguments():
    parser = argparse.ArgumentParser(description='>>> Generates an avg_grid.cdf file (CDPKIT format) from a grid.cdf file.')
    parser.add_argument('-gcdf',
                        dest='gcdf',
                        required=True,
                        help='[Required] The path of the grid.cdf file',
                        nargs=1)

    parser.add_argument('-o',
                        dest='output',
                        help='[Optional] The output folder where the average_grid.cdf file will be generated (Default: current directory)',
                        nargs=1,
                        default=None)
    parser.add_argument('-fr',
                        dest='frame_range',
                        help='[Optional] Frame range to consider, eg [0,55] (Default: all frames)',
                        nargs=1,
                        default=None)

    parse_args = parser.parse_args()
    return parse_args


def writeAverageGrids(gcdf, frame_range, output):
    initial_time = time.time()
    grid_set = Grid.DRegularGridSet()
    cdf_grid_reader = Grid.FileCDFDRegularGridSetReader(gcdf)
    if frame_range[1] == -1:
        frame_range[1] = cdf_grid_reader.getNumRecords()

    if not cdf_grid_reader.read(grid_set):
        print('> CDF reader failure.')

    new_grids = Grid.DRegularGridSet()
    for grid in grid_set:
        g1 = grid.getSize1()
        g2 = grid.getSize2()
        g3 = grid.getSize3()
        g_type = Grid.getName(grid)

        new_grids.addElement(Grid.DRegularGrid(g1, g2, g3).assign(grid))
        Grid.setName(new_grids.getLastElement(), "Avg_" + g_type)

        new_grids.addElement(Grid.DRegularGrid(g1, g2, g3).assign(grid))
        Grid.setName(new_grids.getLastElement(), "Std_" + g_type)

        for k in range(g3):
            for i in range(g1):
                for j in range(g2):
                    new_grids.getLastElement()[i, j, k] = 0

    print('> Calculating average grids...')
    index_frame = 0

    while True:
        grid_set.clear()
        if not cdf_grid_reader.read(grid_set):
            break

        if frame_range[0] <= index_frame <= frame_range[1]:
            n = 0
            for grid in grid_set:
                tmp_mean = Math.DGrid(grid.getSize1(), grid.getSize2(), grid.getSize3())
                g_mean = new_grids[2 * n]
                g_v = new_grids[2 * n + 1]

                tmp_mean.assign(g_mean + (grid - g_mean) / (index_frame + 1))

                g_v.assign(g_v + Math.elemProd((grid - g_mean), (grid - tmp_mean)))
                g_mean.assign(tmp_mean)

                n += 1
        index_frame += 1

    n = 0
    for grid in grid_set:
        g1 = grid.getSize1()
        g2 = grid.getSize2()
        g3 = grid.getSize3()
        g_v = new_grids[2 * n + 1]

        for k in range(g3):
            for i in range(g1):
                for j in range(g2):
                    variance = round(g_v(i, j, k) / (index_frame + 1),
                                     10)  #### Due to calculations, negative values closer to e-10 than 0 are rounded to avoid warnings
                    g_v[i, j, k] = np.sqrt(variance)

        n += 1

    out_fname = os.path.splitext(os.path.basename(gcdf))[0] + '_average_grid.cdf'
    out_path = os.path.join(output, out_fname)
    grid_set_writer = Grid.FileCDFDRegularGridSetWriter(out_path)

    if not grid_set_writer.write(new_grids):
        print('> Could not write to output file')

    calc_time = time.time() - initial_time
    print('> Average grids generated in {}s'.format(int(calc_time)))

    
if __name__ == '__main__':
    args = parseArguments()

    gcdf = args.gcdf[0]

    if args.output is None:
        output = './'
    else:
        output = args.output[0]

    if output[-1] != '/':
        output += '/'

    if args.frame_range is None:
        frame_range = [0,-1]
    else:
        frame_range = args.frame_range[0]

    writeAverageGrids(gcdf, frame_range, output)




