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
import math
import time
import Util
import CDPL.Base as Base
import CDPL.Chem as Chem
import CDPL.Grid as Grid
import CDPL.Pharm as Pharm
import CDPL.Biomol as Biomol
import CDPL.Math as Math


GRID_STEP_SIZE = 0.75
DROP_NON_STD_RESIDUES = True

FTYPE_TO_NAME = { Pharm.FeatureType.H_BOND_ACCEPTOR : 'HBA', \
                  Pharm.FeatureType.H_BOND_DONOR : 'HBD', \
                  Pharm.FeatureType.HYDROPHOBIC : 'H', \
                  Pharm.FeatureType.POS_IONIZABLE : 'PI', \
                  Pharm.FeatureType.NEG_IONIZABLE : 'NI', \
                  Pharm.FeatureType.AROMATIC : 'AR' }

def process():
    if len(sys.argv) < 4:
        print >> sys.stderr, 'Usage:', sys.argv[0], '[input CDF-file] [output grid CDF-file] [ligand] [[residue subset]] [[chunk index] [num. chunks]]'
        sys.exit(2)

    print >> sys.stderr, '- Processing CDF-file:', sys.argv[1], '...'

    mol = Util.loadCDFMolecule(sys.argv[1])

    if not mol:
        print >> sys.stderr, '!! Could not read file'
        sys.exit(2)

    residues = Biomol.ResidueList(mol)
    residue_subset = None

    print >> sys.stderr, '- Num. residues:', residues.getSize()

    if len(sys.argv) == 5 or len(sys.argv) == 7:
        res_subset_ids = Util.toIntegerList(Util.readLines(sys.argv[4]))

        print >> sys.stderr, '- Residue subset:', res_subset_ids

        residue_subset = Chem.FragmentList(residues)

        Util.filterResidues(residue_subset, res_subset_ids)

    else:
        residue_subset = residues

    num_confs = Chem.getNumConformations(mol)

    print >> sys.stderr, '- Num. frames:', num_confs

    chunk_index = 0
    num_chunks = 1
 
    if  len(sys.argv) == 6:
        chunk_index = int(sys.argv[4])
        num_chunks = int(sys.argv[5])
       
    elif len(sys.argv) == 7:
        chunk_index = int(sys.argv[5])
        num_chunks = int(sys.argv[6])
      
    chunk_size = float(num_confs ) / num_chunks
    conf_start_idx = int(chunk_size * chunk_index)
    conf_end_idx = num_confs
    
    if chunk_index >= (num_chunks - 1):
        conf_end_idx = num_confs
    else:
        conf_end_idx = int(chunk_size * (chunk_index + 1))

    print >> sys.stderr, '- Calculating grid bounds ...'

    i = 0
    bbox_min = Math.Vector3D() 
    bbox_max = Math.Vector3D() 
    init = True

    while i < num_confs:
        coords_func = Chem.AtomConformer3DCoordinatesFunctor(i)

        for res in residue_subset:
            Chem.calcBoundingBox(res, bbox_min, bbox_max, coords_func, init)
            init = False

        i += 1
      
    grid_w = bbox_max[0] - bbox_min[0]
    grid_h = bbox_max[1] - bbox_min[1]
    grid_d = bbox_max[2] - bbox_min[2]

    num_x_grid_pts = int(math.ceil(grid_w / GRID_STEP_SIZE))
    num_y_grid_pts = int(math.ceil(grid_h / GRID_STEP_SIZE))
    num_z_grid_pts = int(math.ceil(grid_d / GRID_STEP_SIZE))

    print >> sys.stderr, '- Grid resolution =', GRID_STEP_SIZE
    print >> sys.stderr, '- Residue bounding box size =', str(str(grid_w) + ' x ' + str(grid_h) + ' x ' + str(grid_d) + ' Ang')
    print >> sys.stderr, '- Grid box size =', str(str(num_x_grid_pts) + ' x ' + str(num_y_grid_pts) + ' x ' + str(num_z_grid_pts))
    print >> sys.stderr, '- Total num. grid points =', (num_x_grid_pts * num_y_grid_pts * num_z_grid_pts)

    print >> sys.stderr, '- Extracting binding site environment residues ...'

    bsite_env = Chem.Fragment()
    water_env = Chem.Fragment()
    ligand = Chem.Fragment()
    num_res = 0

    for res in residues:
        if Biomol.getResidueCode(res) == 'HOH':
            water_env += res
            continue

        if Biomol.getResidueCode(res) == sys.argv[3]:
            ligand += res
            continue

        if DROP_NON_STD_RESIDUES and not Biomol.ResidueDictionary.isStdResidue(Biomol.getResidueCode(res)):
            continue

        i = 0

        while i < num_confs:
            coords_func = Chem.AtomConformer3DCoordinatesFunctor(i)
            
            if Chem.intersectsBoundingBox(res, bbox_min, bbox_max, coords_func):
                bsite_env += res
                num_res += 1
                break

            i += 1

    print >> sys.stderr, '-', num_res, 'Residues covered by grid calculation,', bsite_env.numAtoms, 'atoms'

    for bond in mol.bonds:
        if bsite_env.containsAtom(bond.getBegin()) and  bsite_env.containsAtom(bond.getEnd()):
            bsite_env.addBond(bond)

    Chem.perceiveSSSR(bsite_env, False)

    env_pharm = Pharm.BasicPharmacophore()
    pharm_gen = Pharm.DefaultPharmacophoreGenerator(True)
    #pharm_gen.setFeatureGenerator(Pharm.FeatureType.HYDROPHOBIC, Pharm.HydrophobicAtomFeatureGenerator()) # Atom-based hydrophobics!

    grid_ctr = 0.5 * (bbox_min + bbox_max)
    grid_xform = Math.Matrix4D(Math.DTranslationMatrix(4, grid_ctr[0], grid_ctr[1], grid_ctr[2]))

    grid_calc = Pharm.DefaultInteractionScoreGridSetCalculator(GRID_STEP_SIZE, num_x_grid_pts, num_y_grid_pts, num_z_grid_pts)
    grid_calc.setCoordinatesTransform(grid_xform)

    density_calc = Chem.AtomDensityGridCalculator()
    buriedness_calc = Chem.BuriednessGridCalculator()
    
    grid_set = Grid.DRegularGridSet()
    grid_set_writer = Grid.FileCDFDRegularGridSetWriter(sys.argv[2])

    print >> sys.stderr, '- Calculating interaction grids for frames', str('[' + str(conf_start_idx) + ', ' + str(conf_end_idx) + '['), '...' 

    water_atoms = Chem.Fragment()
    i = conf_start_idx 

    pharm_gen_times = []
    num_pharm_ftrs = []
    int_grid_calc_times = []
    water_grid_calc_times = []
    env_grid_calc_times = []
    num_water_atoms = []
    tot_calc_times = []

    while i < conf_end_idx:
        print >> sys.stderr, '- Calculating grids for frame', i, '...'
        st = time.clock()

        print >> sys.stderr, ' - Generating pharmacophore ...'

        st1 = time.clock()
        env_pharm.clear()

        coords_func = Chem.AtomConformer3DCoordinatesFunctor(i)

        pharm_gen.setAtom3DCoordinatesFunction(coords_func)
        pharm_gen.generate(bsite_env, env_pharm)

        pharm_gen_times.append(time.clock() - st1)
        num_pharm_ftrs.append(env_pharm.getNumFeatures())

        print >> sys.stderr, ' - Generated', env_pharm.getNumFeatures(), 'features ...'

        for feature in env_pharm:
            if Pharm.getType(feature) == Pharm.FeatureType.HYDROPHOBIC:
                Pharm.setWeight(feature, Pharm.getHydrophobicity(feature))

        print >> sys.stderr, ' - Calculating interaction grids ...'

        st1 = time.clock()
        grid_set.clear()
        grid_calc.calculate(env_pharm, grid_set)

        for grid in grid_set:
            int_descr = FTYPE_TO_NAME[Pharm.getFeatureType(grid)] + '-' + FTYPE_TO_NAME[Pharm.getTargetFeatureType(grid)]
            Grid.setName(grid, int_descr)

        int_grid_calc_times.append(time.clock() - st1)

        st1 = time.clock()
        env_grid = Grid.DRegularGrid(GRID_STEP_SIZE);
        env_grid.resize(num_x_grid_pts, num_y_grid_pts, num_z_grid_pts, False)
        env_grid.setCoordinatesTransform(grid_xform)
        Grid.setName(env_grid, 'ENV')

        print >> sys.stderr, ' - Calculating environment grid for', bsite_env.getNumAtoms(), 'atoms ...'

        density_calc.setAtom3DCoordinatesFunction(coords_func)
        density_calc.calculate(bsite_env, env_grid)

        env_grid_calc_times.append(time.clock() - st1)

        num_elem = env_grid.getNumElements()

        for grid in grid_set:
            j = 0

            while j < num_elem:
                grid.setElement(j, grid.getElement(j) * (1.0 - env_grid.getElement(j)))
                j += 1
            
        grid_set.addElement(env_grid)

        st1 = time.clock()
        water_atoms.clear()

        for atom in water_env.atoms:
            atom_pos = Chem.getConformer3DCoordinates(atom, i)

            if atom_pos(0) <= bbox_max(0) and atom_pos(0) >= bbox_min(0) and atom_pos(1) <= bbox_max(1) and atom_pos(1) >= bbox_min(1) and atom_pos(2) <= bbox_max(2) and atom_pos(2) >= bbox_min(2):
                water_atoms.addAtom(atom)

        num_water_atoms.append(water_atoms.getNumAtoms())

        grid = Grid.DRegularGrid(GRID_STEP_SIZE);
        grid.resize(num_x_grid_pts, num_y_grid_pts, num_z_grid_pts, False)
        grid.setCoordinatesTransform(grid_xform)
        Grid.setName(grid, 'H2O')
        grid_set.addElement(grid)

        print >> sys.stderr, ' - Calculating water grid for', water_atoms.getNumAtoms(), 'atoms ...'

        density_calc.calculate(water_atoms, grid)

        water_grid_calc_times.append(time.clock() - st1)

        if ligand.getNumAtoms() > 0:
            grid = Grid.DRegularGrid(GRID_STEP_SIZE);
            grid.resize(num_x_grid_pts, num_y_grid_pts, num_z_grid_pts, False)
            grid.setCoordinatesTransform(grid_xform)
            Grid.setName(grid, 'LIG')
            grid_set.addElement(grid)

            print >> sys.stderr, ' - Calculating ligand grid for', ligand.getNumAtoms(), 'atoms ...'

            density_calc.calculate(ligand, grid)

        if not grid_set_writer.write(grid_set):
            print >> sys.stderr, '!! Could not write to output file'
            sys.exit(2)

        tot_time = time.clock() - st
        tot_calc_times.append(tot_time)

        print >> sys.stderr, '- Calculated grid for frame', i, str('(' + str(tot_time) + ' sec)') 

        i += 1

    print >> sys.stderr, '# Statistics:'
    print >> sys.stderr, '- Avg. num. features:', str(float(sum(num_pharm_ftrs)) / len(num_pharm_ftrs))
    print >> sys.stderr, '- Avg. pharm. gen. time:', str(sum(pharm_gen_times) / len(pharm_gen_times)), 'sec'
    print >> sys.stderr, '- Avg. int. grid. calc. time:', str(sum(int_grid_calc_times) / len(int_grid_calc_times)), 'sec'
    print >> sys.stderr, '- Avg. env. grid. calc. time:', str(sum(env_grid_calc_times) / len(env_grid_calc_times)), 'sec'
    print >> sys.stderr, '- Avg. num. water atoms:', str(float(sum(num_water_atoms)) / len(num_water_atoms))
    print >> sys.stderr, '- Avg. water. grid. calc. time:', str(sum(water_grid_calc_times) / len(water_grid_calc_times)), 'sec'
    print >> sys.stderr, '- Avg. total calc. and output time:', str(sum(tot_calc_times) / len(tot_calc_times)), 'sec'


if __name__ == '__main__':
    process()
