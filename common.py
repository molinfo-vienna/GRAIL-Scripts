"""
@author: Arthur Garon
Department of Pharmaceutical Chemistry
University of Vienna
"""

import time
import math
import os
import sys

import CDPL.Base as Base
import CDPL.Biomol as Biomol
import CDPL.Chem as Chem
import CDPL.Grid as Grid
import CDPL.Math as Math
import CDPL.Pharm as Pharm
import CDPL.MolProp as MolProp
import CDPL.GRAIL as GRAIL


xrange = range

GRID_STEP_SIZE = 0.75
DROP_NON_STD_RESIDUES = True


def loadCDFMolecule(fname):
    mol = Chem.BasicMolecule()
    cdf_reader = Chem.FileCDFMoleculeReader(fname)

    try:
        if not cdf_reader.read(mol):
            sys.exit('!!! Could not load CDF file: ' + fname)
    except:
        sys.exit('!!! Could not load CDF file: ' + fname)
            
    Chem.calcImplicitHydrogenCounts(mol, False)
    Chem.perceiveHybridizationStates(mol, False)
    Chem.setAtomSymbolsFromTypes(mol, False)
    Chem.perceiveSSSR(mol, False)
    Chem.setRingFlags(mol, False)
    Chem.setAromaticityFlags(mol, False)
    MolProp.calcAtomHydrophobicities(mol, False)

    return mol

def cdfMol_pdb(pdb, output, name, write_res=True):
    initial_time = time.time()
    pdb_mol = Chem.BasicMolecule()

    pdb_str = open(pdb, 'r').read().replace('WAT', 'HOH').replace('HIE', 'HIS').replace('SPC X', 'HOH X')
    pdb_reader = Biomol.PDBMoleculeReader(Base.StringIOStream(pdb_str))

    Biomol.setPDBApplyDictAtomBondingToNonStdResiduesParameter(pdb_reader, True)

    try:
        if not pdb_reader.read(pdb_mol):
            sys.exit('!!! Could not load PDB file: ' + pdb)
    except:
        sys.exit('!!! Could not load PDB file: ' + pdb)
            
    Pharm.prepareForPharmacophoreGeneration(pdb_mol)

    for atom in pdb_mol.atoms:
        array = Math.Vector3DArray()
        array.addElement(Chem.get3DCoordinates(atom))
        Chem.set3DCoordinatesArray(atom, array)
        
    tmp_output = output + name + ".cdf"
    
    try:
        if not Chem.FileCDFMolecularGraphWriter(tmp_output).write(pdb_mol):
            sys.exit('!!! Could not write CDF file: ' + tmp_output)
    except:
        sys.exit('!!! Could not write CDF file: ' + tmp_output)

    if write_res:
        residues = Biomol.ResidueList(pdb_mol)
        tmp_output = output + name + "_aa_residue_list.txt"

        with open(tmp_output, 'w') as txt_writer:
            for res in residues:
                txt_writer.write(getResidueID(res) + '\n')

    calc_time = time.time() - initial_time
    
    print('> Cdf file generated in {}s'.format(int(calc_time)))

def getResidueID(res):
    if Biomol.getChainID(res) != ' ':
        return Biomol.getResidueCode(res) + '_' + str(Biomol.getResidueSequenceNumber(res)) + '_' + Biomol.getChainID(res)

    return Biomol.getResidueCode(res) + '_' + str(Biomol.getResidueSequenceNumber(res))

def filterResidues(residues, res_subset):
    if len(res_subset) == 0:
        return

    i = 0
    while i < len(residues):
        seq_id = getResidueID(residues[i])
        if seq_id in res_subset:
            i += 1
            continue

        del residues[i]

def getGridInfo(cdf_file, res_file):
    mol = loadCDFMolecule(cdf_file)
    residues = Biomol.ResidueList(mol)

    if res_file:
        with open(res_file) as f:
            res_subset_ids = sorted([line.strip(' ').rstrip('\n') for line in f])

            residue_subset = Chem.FragmentList(residues)
            filterResidues(residue_subset, res_subset_ids)
    else:
        residue_subset = residues

    num_confs = Chem.getNumConformations(mol)
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

    for i in range(0, 3):
        bbox_min[i] = round(bbox_min[i], 2)

    for i in range(0, 3):
        bbox_max[i] = round(bbox_max[i], 2)
        
    return bbox_min,bbox_max

def convertToBBoxMinMax(grid_spec_str):
    coords = grid_spec_str.split(',')

    bbox_min = Math.Vector3D()
    bbox_max = Math.Vector3D()

    for i in range(0, 3):
        bbox_min[i] = float(coords[i])

    for i in range(0, 3):
        bbox_max[i] = float(coords[i + 3])
    
    return bbox_min, bbox_max
    
def grailGeneration(cdf_file, lig, bbox_min, bbox_max, output, bur):
    def setHydrophicFtrWeights(pharm):
        for feature in pharm:
            if Pharm.getType(feature) == Pharm.FeatureType.HYDROPHOBIC:
                Pharm.setWeight(feature, Pharm.getHydrophobicity(feature))
                
    initial_time = time.time()

    grail_gen = GRAIL.GRAILDataSetGenerator()

    grail_gen.setGridStepSize(GRID_STEP_SIZE)
    grail_gen.setGridParamsForBoundingBox(bbox_min, bbox_max)
    grail_gen.pharmGenerator.applyConfiguration(Pharm.DefaultPharmacophoreGenerator.DEFAULT_CONFIG)
    grail_gen.pharmProcessingFunction = setHydrophicFtrWeights
    
    mol = loadCDFMolecule(cdf_file)
    num_confs = Chem.getNumConformations(mol)
  
    residues = Biomol.ResidueList(mol)
    bsite_env = Chem.Fragment()
    water_env = Chem.Fragment()
    ligand = Chem.Fragment()
    num_res = 0
    
    for res in residues:
        if Biomol.getResidueCode(res) == 'HOH':
            water_env += res
            continue

        if Biomol.getResidueCode(res) == lig:
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

    print('> Grid size (points): {} x {} x {}'.format(grail_gen.gridXSize, grail_gen.gridYSize, grail_gen.gridZSize))
    print('> Total number of grid points: {}'.format(grail_gen.gridXSize * grail_gen.gridYSize * grail_gen.gridZSize))
    print('> Spatial grid dimensions (min, max): {},{},{},{},{},{}'.format(bbox_min[0], bbox_min[1], bbox_min[2], bbox_max[0], bbox_max[1], bbox_max[2]))
    print('> Number of residues covered by grid calculation: {}'.format(num_res))

    for bond in mol.bonds:
        if bsite_env.containsAtom(bond.getBegin()) and bsite_env.containsAtom(bond.getEnd()):
            bsite_env.addBond(bond)

    Chem.perceiveSSSR(bsite_env, False)

    buriedness_calc = GRAIL.BuriednessGridCalculator()

    grid_set = Grid.DRegularGridSet()
    grid_set_writer = Grid.FileCDFDRegularGridSetWriter(
        output + '/' + os.path.basename(cdf_file)[:-4] + '_grid.cdf')

    water_atoms = Chem.Fragment()
    i = 0

    while i < num_confs:
        print('> Calculating grids for frame {} ...'.format(i))
  
        coords_func = Chem.AtomConformer3DCoordinatesFunctor(i)
        
        grail_gen.calcInteractionGrids(bsite_env, coords_func, grid_set)

        water_atoms.clear()

        for atom in water_env.atoms:
            atom_pos = Chem.getConformer3DCoordinates(atom, i)

            if bbox_min(0) <= atom_pos(0) <= bbox_max(0) and bbox_min(1) <= atom_pos(1) <= bbox_max(1) and bbox_min(
                    2) <= atom_pos(2) <= bbox_max(2):
                water_atoms.addAtom(atom)

        grid_set.addElement(grail_gen.calcAtomDensityGrid(water_atoms, coords_func, 'H2O'))

        if ligand.getNumAtoms() > 0:
            grid_set.addElement(grail_gen.calcAtomDensityGrid(ligand, coords_func, 'LIG'))

        if bur:
            bur_grid = Grid.DRegularGrid(GRID_STEP_SIZE)
            bur_grid.resize(grail_gen.gridXSize, grail_gen.gridYSize, grail_gen.gridZSize, False)
            bur_grid.setCoordinatesTransform(grail_gen.gridTransform)

            Grid.setName(bur_grid, 'BUR')

            grid_set.addElement(bur_grid)

            buriedness_calc.setAtom3DCoordinatesFunction(coords_func)
            buriedness_calc.calculate(bsite_env, bur_grid)

        if not grid_set_writer.write(grid_set):
            print('> Could not write to output file')
            return None

        i += 1
                                
    print('> Grid.cdf file generated in {}s'.format(int(time.time() - initial_time)))

def calcAtomSetCentroid(atoms, conf_idx):
    if len(atoms) == 1:
        return Chem.getConformer3DCoordinates(atoms[0], conf_idx)

    ctr = Math.Vector3D()

    for atom in atoms:
        ctr += Chem.getConformer3DCoordinates(atom, conf_idx)

    ctr /= len(atoms)
    return ctr

def kont_conversion(gcdf, index_frame, output):
    cdf_grid_reader = Grid.FileCDFDRegularGridSetReader(gcdf)
    grid_set = Grid.DRegularGridSet()

    if not cdf_grid_reader.read(index_frame, grid_set):
        print('> No grids for the frame number.')
        return

    for grid in grid_set:
        int_descr = Grid.getName(grid)
        tmp_output = output + '/' + os.path.basename(gcdf)[:-4] + '_' + int_descr + '.kont'
        with open(tmp_output, 'w') as kont_file:
            index_mol = 0
            grid_values = []

            for k in range(grid.getSize3()):
                for i in range(grid.getSize1()):
                    for j in range(grid.getSize2()):
                        index_mol += 1
                        grid_values.append(grid(i, j, k))
                        coord = Math.Vector3D()
                        grid.getCoordinates(i, j, k, coord)
                        kont_file.write(
                            "{:>7}   {:>7.3f} {:>7.3f} {:>7.3f}\n".format(index_mol, coord[0], coord[1], coord[2]))

            kont_file.write("{:>8}\n".format(int_descr))
            for txt in grid_values:
                kont_file.write("{:>8.3f}\n".format(txt))
