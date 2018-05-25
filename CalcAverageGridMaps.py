"""
@author: Arthur Garon
department of pharmaceutical chemistry
university of vienna
"""

import sys
import os.path as path
import CDPL.Grid as Grid
import CDPL.Math as Math
import numpy as np

def process():
    if len(sys.argv) < 3:
        print >> sys.stderr, 'Usage:', sys.argv[0], '[input.cdf] [output directory]'
        sys.exit(2)

    print >> sys.stderr, '> Processing grid CDF-file:', sys.argv[1], '...'

    grid_set = Grid.DRegularGridSet()
    cdf_grid_reader = Grid.FileCDFDRegularGridSetReader(sys.argv[1])

    if not cdf_grid_reader.read(grid_set):
        print '> CDF reader failure.'

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

    print >> sys.stderr, '> Calculating averaged grids...'

    index_frame = 0

    while True:
        grid_set.clear()
        if not cdf_grid_reader.read(grid_set):
            break
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
                    variance = round(g_v(i, j, k) / (index_frame+1), 10) #### Due to calculations, negative values closer than e-10 from 0 are rounded to avoid warnings
                    g_v[i, j, k] = np.sqrt(variance)

        n += 1

    print >> sys.stderr, '> Writing averaged grids to', sys.argv[2]

    out_fname = path.splitext(path.basename(sys.argv[1]))[0] + '_average_grid.cdf'
    out_path = path.join(sys.argv[2], out_fname)
    grid_set_writer = Grid.FileCDFDRegularGridSetWriter(out_path)

    if not grid_set_writer.write(new_grids):
        print >> sys.stderr, '> Could not write output file'

if __name__ == '__main__':
    process()
