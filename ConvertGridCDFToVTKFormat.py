# -*- mode: python; tab-width: 4 -*- 

## 
# This file is part of the Chemical Data Processing Toolkit
#
# Copyright (C) 2003-2008 Thomas A. Seidel <thomas.seidel@univie.ac.at>
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program; see the file COPYING. If not, write to
# the Free Software Foundation, Inc., 59 Temple Place - Suite 330,
# Boston, MA 02111-1307, USA.
##


import sys
import os.path as path
import CDPL.Chem as Chem
import CDPL.Grid as Grid
import CDPL.Math as Math
import Util
import numpy
import pyevtk


def process():
    if len(sys.argv) < 3:
        print >> sys.stderr, 'Usage:', sys.argv[0], '[input.cdf] [output directory]'
        sys.exit(2)

    in_fname = path.splitext(path.basename(sys.argv[1]))[0]
    grid_set = Grid.DRegularGridSet()
    cdf_reader = Grid.FileCDFDRegularGridSetReader(sys.argv[1])
    pvd_file = open(path.join(sys.argv[2], in_fname + '.pvd'), 'w')

    Util.writePVDHeader(pvd_file)

    print >> sys.stderr, '- Processing grid CDF-file:', sys.argv[1], '...'

    frame_no = 0

    while True:
        grid_set.clear()

        if not cdf_reader.read(grid_set):
            break

        output_data = {}
        grid_spacing = None
        grid_origin = None

        for grid in grid_set:
            int_descr = Grid.getName(grid)
            grid_values = numpy.ndarray((grid.getSize1(), grid.getSize2(), grid.getSize3()), numpy.float32)
            output_data[int_descr] = grid_values

            for i in range(grid.getSize1()):
                for j in range(grid.getSize2()):
                    for k in range(grid.getSize3()):
                        grid_values[i, j, k] = grid(i, j, k)

            if not grid_spacing:
                grid_spacing = (grid.getXStepSize(), grid.getYStepSize(), grid.getZStepSize())
        
            if not grid_origin:
                grid_origin = Math.Vector3D()
                grid.getCoordinates(0, 0, 0, grid_origin)

        out_fname = in_fname + '_frame_no_' + str(frame_no)
        out_path = path.join(sys.argv[2], out_fname)

        print >> sys.stderr, '- Writing grid data for frame', frame_no, '...'

        if not pyevtk.hl.imageToVTK(out_path, origin = grid_origin, spacing = grid_spacing, pointData = output_data):
            print '!! Could not write grid output file'
            sys.exit(2)

        Util.writePVDEntry(pvd_file, frame_no, out_fname, 'vti')

        frame_no += 1

    Util.writePVDFooter(pvd_file)

if __name__ == '__main__':
    process()
