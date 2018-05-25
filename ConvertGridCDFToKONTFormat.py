"""
@author: Arthur Garon
department of pharmaceutical chemistry
university of vienna
"""

import sys
import os.path as path
import CDPL.Grid as Grid
import CDPL.Math as Math


def process():
    if len(sys.argv) < 4:
        print >> sys.stderr, 'Usage:', sys.argv[0], '[input.cdf] [output directory] [frame index]'
        sys.exit(2)

    print >> sys.stderr, '> Processing grid CDF-file:', sys.argv[1], '...'

    index_frame = int(sys.argv[3])
    grid_set = Grid.DRegularGridSet()
    cdf_grid_reader = Grid.FileCDFDRegularGridSetReader(sys.argv[1])

    print >> sys.stderr, '> Writing grid data for frame', index_frame, '...\n'

    grid_set.clear()
    if not cdf_grid_reader.read(index_frame, grid_set):
        print '> No grid in file.'
        sys.exit(2)

    coord = Math.Vector3D()

    for grid in grid_set:
        int_descr = Grid.getName(grid)
        out_fname = path.splitext(path.basename(sys.argv[1]))[0] + '_' + int_descr + '_frame_no_' + sys.argv[3] + '.kont'
        out_path = path.join(sys.argv[2], out_fname)
        kont_file = open(out_path, 'w')

        index_mol = 0
        grid_values = []
        print '>', int_descr, 'grid'

        for k in range(grid.getSize3()):
            for i in range(grid.getSize1()):
                for j in range(grid.getSize2()):
                    index_mol += 1
                    grid_values.append(grid(i, j, k))
                    grid.getCoordinates(i, j, k, coord)
                    kont_file.write("{:>7}   {:>7.3f} {:>7.3f} {:>7.3f}\n".format(index_mol, coord[0], coord[1], coord[2]))

        kont_file.write("{:>8}\n".format(int_descr))

        for txt in grid_values:
            kont_file.write("{:>8.3f}\n".format(txt))

        kont_file.close()


if __name__ == '__main__':
    process()
