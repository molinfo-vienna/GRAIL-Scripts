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

import MDAnalysis
import CDPL.Chem as Chem
import CDPL.Biomol as Biomol
import CDPL.Math as Math
import numpy as np
import sys
import math


def process():
    if len(sys.argv) < 4:
        print >> sys.stderr, 'Usage:', sys.argv[0], '[input topology-file] [input coordinates-file] [output CDF-file]'
        sys.exit(2)

    print >> sys.stderr, '- Processing topology-file', sys.argv[1], 'and coordinates-file', sys.argv[2], '...'

    u = MDAnalysis.Universe(sys.argv[1], sys.argv[2])
    cdf_mol = Chem.BasicMolecule()

    cdf_mol.reserveMemoryForAtoms(len(u.atoms))
    cdf_mol.reserveMemoryForBonds(len(u.bonds))

    print >> sys.stderr, '- Num. atoms:', len(u.atoms) 
    print >> sys.stderr, '- Num. bonds:', len(u.bonds) 

    num_frames = len(u.trajectory) 

    print >> sys.stderr, '- Num. frames:', num_frames

    # construct atoms

    print >> sys.stderr, '- Building atoms ...'

    waters = {}
    i = 0

    for md_atom in u.atoms:
        atom = cdf_mol.addAtom()
        sym = MDAnalysis.topology.guessers.guess_atom_element(md_atom.name)

        Chem.setSymbol(atom, sym.title())
        Chem.setImplicitHydrogenCount(atom, 0)
        Biomol.setChainID(atom, md_atom.segid)
        
        if md_atom.resname == 'WAT':
            Biomol.setResidueCode(atom, 'HOH')
        else:
            Biomol.setResidueCode(atom, md_atom.resname)

        if Biomol.getResidueCode(atom) == 'HOH':
            if md_atom.resid in waters:
                waters[md_atom.resid].append(i)
            else:
                waters[md_atom.resid] = [ i ]

        Biomol.setResidueSequenceNumber(atom, int(md_atom.resid))
        Biomol.setResidueAtomName(atom, md_atom.name)

        # fix positive charge on arginin nitrogen
        if md_atom.resname == 'ARG' and md_atom.name == 'NH2':
            Chem.setFormalCharge(atom, 1)
    
        coords = []
        for coord in md_atom.position:
            coords.append(float(coord))
        
        Chem.set3DCoordinates(atom, coords)

        coords_array = Math.Vector3DArray()
        coords_array.reserve(num_frames)

        Chem.set3DCoordinatesArray(atom, coords_array)
        Chem.setPEOECharge(atom, float(md_atom.charge))
    
        i += 1

    Chem.setAtomTypesFromSymbols(cdf_mol, True)

    # construct bonds

    print >> sys.stderr, '- Building bonds ...'

    for md_bond in u.bonds:
        cdf_mol.addBond(int(md_bond.atoms[0].index), int(md_bond.atoms[1].index))

    print >> sys.stderr, '- Building water atom bonds ...'

    for water in waters.values():
        if len(water) < 2:
            continue

        for atom_idx in water:
            if Chem.getType(cdf_mol.atoms[atom_idx]) == Chem.AtomType.O:
                if atom.numBonds > 1:
                    break
                    
                for atom_idx2 in water:
                    if Chem.getType(cdf_mol.atoms[atom_idx2]) == Chem.AtomType.H:
                        cdf_mol.addBond(atom_idx, atom_idx2)

                break

    # make sane biomolecule

    Chem.perceiveSSSR(cdf_mol, True)
    Chem.setRingFlags(cdf_mol, True)
    Chem.perceiveBondOrders(cdf_mol, True)
    Chem.perceiveHybridizationStates(cdf_mol, True)
    Chem.setAromaticityFlags(cdf_mol, True)
    Chem.calcFormalCharges(cdf_mol, True)

    # read timsteps and write cdf

    print >> sys.stderr, '- Importing coordinates ...'

    i = 0
    traj_coords = []
    atom_coords = Math.Vector3D()

    for ts in u.trajectory:
        print >> sys.stderr, '- Processing time step', i, '...'

        for md_atom in u.atoms:
            del traj_coords[:]

            for coord in md_atom.position:
                traj_coords.append(float(coord))

            coords_array = Chem.get3DCoordinatesArray(cdf_mol.getAtom(int(md_atom.index)))        
           
            atom_coords[0] = traj_coords[0]
            atom_coords[1] = traj_coords[1]
            atom_coords[2] = traj_coords[2]

            coords_array.addElement(atom_coords)
           
        i += 1

    print >> sys.stderr, '- Writing output file:'

    if not Chem.FileCDFMolecularGraphWriter(sys.argv[3]).write(cdf_mol):
        print >> sys.stderr, '!! Could not write file'
        sys.exit(2)


if __name__ == '__main__':
    process()
