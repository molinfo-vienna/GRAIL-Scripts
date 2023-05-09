##
# Calculates the grid point vaule differences for two sequences of input grid sets
#
# Author: Thomas Seidel, Department of Pharmaceutical Sciences, University of Vienna
##

import sys
import argparse

import CDPL.Grid as Grid


def parseArguments():
    parser = argparse.ArgumentParser(description='Calculates happy/unhappy water regions for a series of input grids.')

    parser.add_argument('-i',
                        dest='input_file',
                        required=True,
                        help='Paths to the input grid set files (in CDF-format)',
                        nargs=2)

    parser.add_argument('-o',
                        dest='output_file',
                        required=True,
                        help='The output grid set file (written in CDF-format)',
                        nargs=1)

    parsed_args = parser.parse_args()

    return parsed_args
    
# called for each read input grid sets
def genOutputGridSet(in_grid_set1, in_grid_set2, out_grid_set):
    out_grid_set.clear()

    if len(in_grid_set1) != len(in_grid_set2):
        sys.exit('Error: differently sized input grid sets')

    for i in range(0, len(in_grid_set1)):
        in_grid1 = in_grid_set1[i]
        in_grid2 = in_grid_set2[i]

        if Grid.getName(in_grid1) != Grid.getName(in_grid2):
            sys.exit('Error: detected different order of grids in input grid sets')
        
        out_grid = Grid.DRegularGrid(in_grid1) # create copy of the first input grid
        out_grid -= in_grid2                   # subtract the second input grid from the copy of the first
        
        Grid.setName(out_grid, "Diff_" + Grid.getName(in_grid1)) # modify name

        out_grid_set.addElement(out_grid) # append to output grid set

# main procedure
def process(input_file1, input_file2, output_file):
    in_grid_set1 = Grid.DRegularGridSet()
    in_grid_set2 = Grid.DRegularGridSet()
    out_grid_set = Grid.DRegularGridSet()
     
    grid_set_reader1 = Grid.FileCDFDRegularGridSetReader(input_file1)
    grid_set_reader2 = Grid.FileCDFDRegularGridSetReader(input_file2)
    grid_set_writer = Grid.FileCDFDRegularGridSetWriter(output_file)

    print('Calculating grid differences...')
    
    while True:
        if not grid_set_reader1.read(in_grid_set1): # read first input grid set
            break

        if not grid_set_reader2.read(in_grid_set2): # read second input grid set
            break

        genOutputGridSet(in_grid_set1, in_grid_set2, out_grid_set)
        
        if not grid_set_writer.write(out_grid_set):
            sys.exit('Error: could not write output grids')
            
        print('.', end='')  # print a dot for each frame as progress indicator
        sys.stdout.flush()  # flush output buffer

    grid_set_writer.close()

    print('\nDone!')
      
if __name__ == '__main__':
    args = parseArguments()

    process(args.input_file[0], args.input_file[1], args.output_file[0])
