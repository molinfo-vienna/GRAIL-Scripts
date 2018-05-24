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
import CDPL.Biomol as Biomol
import CDPL.Math as Math
import Util
import scipy.interpolate as si
import numpy
import pyevtk


SPLINE_POINTS_PER_BB_ATOM = 10
GRID_BOX_MIN_X = -14.5 - 5
GRID_BOX_MIN_Y = 42.0 - 5
GRID_BOX_MIN_Z = 4.26 - 5
GRID_BOX_MAX_X = -4.5 + 5
GRID_BOX_MAX_Y = 53.5 + 5
GRID_BOX_MAX_Z = 18.26 + 5


def insideBBox(bbox_min, bbox_max, x, y, z):
    return (x <= bbox_max(0) and x >= bbox_min(0) and y <= bbox_max(1) and y >= bbox_min(1) and z <= bbox_max(2) and z >= bbox_min(2))

def spline(cv, npts):
    tck, u_ = si.splprep(cv.T,s=0.0)
    return si.splev(numpy.linspace(0, 1, npts), tck)

def process():
    if len(sys.argv) < 3:
        #print >> sys.stderr, 'Usage:', sys.argv[0], '[input.cdf] [output directory] [[ligand] [ligand ctab]]'
        print >> sys.stderr, 'Usage:', sys.argv[0], '[input.cdf] [output directory]'
        sys.exit(2)

    in_fname = path.splitext(path.basename(sys.argv[1]))[0]
    mol = Chem.BasicMolecule()
    cdf_reader = Chem.FileCDFMoleculeReader(sys.argv[1])
    pvd_file = open(path.join(sys.argv[2], in_fname + '.pvd'), 'w')

    Util.writePVDHeader(pvd_file)

    print >> sys.stderr, '- Processing CDF-file:', sys.argv[1], '...'

    if not cdf_reader.read(mol):
        print '!! Could not read file'
        sys.exit(2)

    backbone_atoms = []

    if False: #len(sys.argv) > 3:
        #lig_atoms = {}
        
        for atom in mol.atoms:
            if Biomol.isPDBBackboneAtom(atom) and Biomol.getResidueAtomName(atom) == 'C':
                backbone_atoms.append(atom)

            #if Biomol.getResidueCode(atom) == sys.argv[3]:
            #    lig_atoms[Biomol.getResidueAtomName(atom)] = atom
                
            #    while atom.getNumBonds() > 0:
            #        mol.removeBond(atom.bonds[0].getIndex())

        #f = file(sys.argv[4], 'r')
        #lines = []

        #for l in f.readlines():
        #    l = l.strip()

        #    if not l:
        #        continue
            
        #    lines.append(l)

        #num_atoms = int(lines[0])
        #i = 1

        #while i <= num_atoms:
        #    tokens = lines[i].split()

        #    if tokens[0] in lig_atoms:
        #        Chem.setType(lig_atoms[tokens[0]], Chem.AtomDictionary.getType(tokens[1]))

        #    i += 1

        #while i < len(lines):
        #    tokens = lines[i].split()
        #    i += 1

        #    if tokens[0] in lig_atoms and tokens[1] in lig_atoms:
        #        Chem.setOrder(mol.addBond(lig_atoms[tokens[0]].getIndex(), lig_atoms[tokens[1]].getIndex()), int(tokens[2]))

        #for atom in lig_atoms.values():
        #    Chem.setImplicitHydrogenCount(atom, Chem.calcImplicitHydrogenCount(atom, mol))
        #    Chem.setHybridizationState(atom, Chem.perceiveHybridizationState(atom, mol))

    else:
        for atom in mol.atoms:
            if Biomol.isPDBBackboneAtom(atom) and Biomol.getResidueAtomName(atom) == 'C':
                backbone_atoms.append(atom)

    #print 'num backbone atoms =', len(backbone_atoms)

    bond_list = []

    for bond in mol.bonds:
        if Biomol.getResidueCode(bond.getAtom(0)) == 'HOH' or Biomol.getResidueCode(bond.getAtom(1)) == 'HOH':
            continue
 
        if Chem.getType(bond.getAtom(0)) == Chem.AtomType.H or Chem.getType(bond.getAtom(1)) == Chem.AtomType.H:
            continue
        
        bond_list.append(bond)

    num_confs = Chem.getNumConformations(mol)

    num_coords = len(bond_list) * 4 + (len(backbone_atoms) * SPLINE_POINTS_PER_BB_ATOM - 1) * 2
    bond_ctr = Math.Vector3D()
    i = 0

    bbox_min = Math.Vector3D() 
    bbox_max = Math.Vector3D() 

    bbox_min[0] = GRID_BOX_MIN_X
    bbox_min[1] = GRID_BOX_MIN_Y
    bbox_min[2] = GRID_BOX_MIN_Z
 
    bbox_max[0] = GRID_BOX_MAX_X
    bbox_max[1] = GRID_BOX_MAX_Y
    bbox_max[2] = GRID_BOX_MAX_Z

    while i < num_confs:
        line_x_coords = numpy.ndarray(num_coords, numpy.float32)
        line_y_coords = numpy.ndarray(num_coords, numpy.float32)
        line_z_coords = numpy.ndarray(num_coords, numpy.float32)
        atom_types = numpy.ndarray(num_coords, numpy.uint32)

        spline_ctrl_points = numpy.ndarray((len(backbone_atoms), 3), numpy.float32)
        j = 0

        for atom in backbone_atoms:
            atom_pos = Chem.getConformer3DCoordinates(atom, i)

            spline_ctrl_points[j, 0] = atom_pos(0)
            spline_ctrl_points[j, 1] = atom_pos(1)
            spline_ctrl_points[j, 2] = atom_pos(2)
            j += 1

        spline_pts = spline(spline_ctrl_points, len(backbone_atoms) * SPLINE_POINTS_PER_BB_ATOM)
        j = 0
        k = 0

        while k < (len(backbone_atoms) * SPLINE_POINTS_PER_BB_ATOM - 1):
            #if not insideBBox(bbox_min, bbox_max, spline_pts[0][k], spline_pts[1][k], spline_pts[2][k]) and not insideBBox(bbox_min, bbox_max, spline_pts[0][k + 1], spline_pts[1][k + 1], spline_pts[2][k + 1]):
            #    k += 1
            #    continue

            line_x_coords[j] = spline_pts[0][k]
            line_y_coords[j] = spline_pts[1][k]
            line_z_coords[j] = spline_pts[2][k]
            atom_types[j] = 0
            j += 1
    
            line_x_coords[j] = spline_pts[0][k + 1]
            line_y_coords[j] = spline_pts[1][k + 1]
            line_z_coords[j] = spline_pts[2][k + 1]
            atom_types[j] = 0
            j += 1
            k += 1

        for bond in bond_list:
            atom1 = bond.getAtom(0)
            atom2 = bond.getAtom(1)
           
            atom1_pos = Chem.getConformer3DCoordinates(atom1, i)
            atom2_pos = Chem.getConformer3DCoordinates(atom2, i)

            #if not insideBBox(bbox_min, bbox_max, atom1_pos(0), atom1_pos(1), atom1_pos(2)) and not insideBBox(bbox_min, bbox_max, atom2_pos(0), atom2_pos(1), atom2_pos(2)):
            #    continue

            atom1_type = Chem.getType(atom1)
            atom2_type = Chem.getType(atom2)

            bond_ctr.assign(atom1_pos)
            bond_ctr += atom2_pos
            bond_ctr *= 0.5

            line_x_coords[j] = atom1_pos(0)
            line_y_coords[j] = atom1_pos(1)
            line_z_coords[j] = atom1_pos(2)
            atom_types[j] = atom1_type
            j += 1

            line_x_coords[j] = bond_ctr(0)
            line_y_coords[j] = bond_ctr(1)
            line_z_coords[j] = bond_ctr(2)
            atom_types[j] = atom1_type
            j += 1
 
            line_x_coords[j] = bond_ctr(0)
            line_y_coords[j] = bond_ctr(1)
            line_z_coords[j] = bond_ctr(2)
            atom_types[j] = atom2_type
            j += 1

            line_x_coords[j] = atom2_pos(0)
            line_y_coords[j] = atom2_pos(1)
            line_z_coords[j] = atom2_pos(2)
            atom_types[j] = atom2_type
            j += 1

        line_x_coords.resize(j)
        line_y_coords.resize(j)
        line_z_coords.resize(j)
        atom_types.resize(j)

        out_fname = in_fname + '_frame_no_' + str(i)
        out_path = path.join(sys.argv[2], out_fname)
        line_data = { 'atom_type' : atom_types }

        print >> sys.stderr, '- Writing structure data for frame', i, '...'

        if not pyevtk.hl.linesToVTK(out_path, line_x_coords, line_y_coords, line_z_coords, pointData = line_data):
            print '!! Could not write structure output file'
            sys.exit(2)

        Util.writePVDEntry(pvd_file, i, out_fname, 'vtu')

        i += 1
        #if i == 4:
        #    break

    Util.writePVDFooter(pvd_file)

if __name__ == '__main__':
    process()
