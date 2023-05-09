##
# Calculates the average of the grid point score values for a sequence of input grid sets 
#
# Author: Thomas Seidel, Department of Pharmaceutical Sciences, University of Vienna
##

import sys
import argparse

import CDPL.Grid as Grid


def parseArguments():
    parser = argparse.ArgumentParser(description='Calculates the arithmetic mean of a series of input grids.')

    parser.add_argument('-i',
                        dest='input_file',
                        required=True,
                        help='Path to the input grid set file (in CDF-format)',
                        nargs=1)

    parser.add_argument('-o',
                        dest='output_file',
                        required=True,
                        help='The output grid set file (written in CDF-format)',
                        nargs=1)

    parsed_args = parser.parse_args()

    return parsed_args

# called to initialize the output grid set
def initOutputGridSet(in_grid_set, out_grid_set):
    for in_grid in in_grid_set:
        out_grid = Grid.DRegularGrid(in_grid) # create copy of the input grid
        
        Grid.setName(out_grid, "Avg_" + Grid.getName(in_grid)) # modify name

        out_grid_set.addElement(out_grid) # append to output grid set

# called for each read input grid set
def processGridSet(in_grid_set, out_grid_set):
    for i in range(0, len(in_grid_set)):
        out_grid_set[i] += in_grid_set[i] # for each grid in the set, accumulate score values at the grid points

# main procedure
def process(input_file, output_file):
    in_grid_set = Grid.DRegularGridSet()
    out_grid_set = Grid.DRegularGridSet()
     
    grid_set_reader = Grid.FileCDFDRegularGridSetReader(input_file)
    num_frames = 0

    print('Averaging input grids...')
    
    while True:
        if not grid_set_reader.read(in_grid_set): # read next input grid set
            break

        if num_frames == 0: # initialize output grid set with a copy of the first input grid set
            initOutputGridSet(in_grid_set, out_grid_set)
        else:               # otherwise, accumulate score values
            processGridSet(in_grid_set, out_grid_set)

        num_frames += 1

    # finally, divide accumulated score values by the number of frames to obtain the arithmetic mean
    if num_frames > 1:      
         for grid in out_grid_set:
             grid /= num_frames

    # write the output grid set to disk
    grid_set_writer = Grid.FileCDFDRegularGridSetWriter(output_file)

    if not grid_set_writer.write(out_grid_set):
        sys.exit('Error: could not write output file')

    grid_set_writer.close()

    print('Done!')
      
if __name__ == '__main__':
    args = parseArguments()

    process(args.input_file[0], args.output_file[0])
