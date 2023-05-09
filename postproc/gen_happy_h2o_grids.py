##
# Calculates happy/unhappy water regions for a sequence of input grid sets
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
                        help='Path to the input grid set file (in CDF-format)',
                        nargs=1)

    parser.add_argument('-o',
                        dest='output_file',
                        required=True,
                        help='The output grid set file (written in CDF-format)',
                        nargs=1)

    parsed_args = parser.parse_args()

    return parsed_args

def getGridByName(grid_set, name):
    for grid in grid_set:
        if Grid.getName(grid) == name:
            return grid
    return None
        
# called for each read input grid set
def genOutputGridSet(in_grid_set, out_grid_set):
    h2o_grid = getGridByName(in_grid_set, 'H2O')  # get H2O atom density grid

    if h2o_grid is None:
        return False

    h_h_grid = getGridByName(in_grid_set, 'H-H')  # get grid for hydrophobic interactions

    if h_h_grid is None:
        return False

    hba_hbd_grid = getGridByName(in_grid_set, 'HBA-HBD') # get grid for lig. H-acceptor -> env. h-donor interactions

    if hba_hbd_grid is None:
        return False

    hbd_hba_grid = getGridByName(in_grid_set, 'HBD-HBA') # get grid for lig. H-donor -> env. h-acceptor interactions

    if hbd_hba_grid is None:
        return False

    happy_h2o_grid = Grid.DRegularGrid(h2o_grid)    # initialize output grid for favorable H2O regions
    unhappy_h2o_grid = Grid.DRegularGrid(h2o_grid)  # initialize output grid for unfavorable H2O regions

    Grid.setName(happy_h2o_grid, 'Happy_H2O')
    Grid.setName(unhappy_h2o_grid, 'Unhappy_H2O')

    out_grid_set.clear()                            # setup output grid set
    out_grid_set.addElement(happy_h2o_grid)
    out_grid_set.addElement(unhappy_h2o_grid)
   
    for i in range(0, h2o_grid.getSize()):          # iterate over all points
        happy_h2o_grid[i] *= max(hba_hbd_grid[i], hbd_hba_grid[i]) * (1.0 - h_h_grid[i])   # H2O,fav(i) = max(HBA-HBD(i), HBD-HBA(i)) * (1 - H-H(i))
        unhappy_h2o_grid[i] *= (1.0 - max(hba_hbd_grid[i], hbd_hba_grid[i])) * h_h_grid[i] # H2O,unfav(i) = (1 - max(HBA-HBD(i), HBD-HBA(i))) * H-H(i)

    return True

# main procedure
def process(input_file, output_file):
    in_grid_set = Grid.DRegularGridSet()
    out_grid_set = Grid.DRegularGridSet()
     
    grid_set_reader = Grid.FileCDFDRegularGridSetReader(input_file)
    grid_set_writer = Grid.FileCDFDRegularGridSetWriter(output_file)

    print('Calculating happy/unhappy water regions...')
    
    while True:
        if not grid_set_reader.read(in_grid_set): # read next input grid set
            break

        if genOutputGridSet(in_grid_set, out_grid_set):
            if not grid_set_writer.write(out_grid_set):
                sys.exit('Error: could not write output grids')
        else:
            print('Warning: generating output grid set failed')
            
        print('.', end='')  # print a dot for each frame as progress indicator
        sys.stdout.flush()  # flush output buffer

    grid_set_writer.close()

    print('\nDone!')
      
if __name__ == '__main__':
    args = parseArguments()

    process(args.input_file[0], args.output_file[0])
