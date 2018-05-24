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
import Util
import CDPL.Base as Base
import CDPL.Chem as Chem
import CDPL.Biomol as Biomol
import CDPL.Math as Math


def process():
    if len(sys.argv) < 4:
        print >> sys.stderr, 'Usage:', sys.argv[0], '[threshold] [input CDF-file] [output CDF-file] [[residue subset]]'
        sys.exit(2)

    print '- Processing CDF-file:', sys.argv[2], '...'

    threshold = float(sys.argv[1])
    mol = Util.loadCDFMolecule(sys.argv[2])

    if not mol:
        print '!! Could not read file'
        sys.exit(2)

    #old_seq_nos, old_to_new_seq_no_map = Util.fixResidueSeqNumbers(mol)
    residues = Biomol.ResidueList(mol)

    print '- Num. residues:', residues.getSize()

    if len(sys.argv) > 4:
        res_subset = Util.toIntegerList(Util.readLines(sys.argv[4]))

        print '- Residue subset:', res_subset

        #for i in range(len(res_subset)):
        #    res_subset[i] = old_to_new_seq_no_map[res_subset[i]]

        Util.filterResidues(residues, res_subset)

    num_confs = Chem.getNumConformations(mol)

    print '- Num. frames:', num_confs
    print '- Calculating residue mobilities...'

    saved_mean_positions = []
    saved_res_positions = []

    for res in residues:
        atoms = Util.getBackboneAtoms(res)
        positions = []
        i = 0

        while i < num_confs:
            positions.append(Util.calcAtomSetCentroid(atoms, i))
            i += 1
        
        mean_pos = Util.calcMeanPos(positions)
        fluct = Util.calcMaxFluct(positions, mean_pos)

        if fluct <= threshold:
            saved_mean_positions.append(mean_pos)
            saved_res_positions.append(positions)

            print 'Using residue: ' + Util.getResidueID(res) + ', fluct: ' + str(fluct)

    if len(saved_mean_positions) < 3:
        print '!! Not enough reference positions for alignment'
        sys.exit(2)     

    alignment = Math.DKabschAlgorithm()
    al_ref_positions = Math.DMatrix(3, len(saved_mean_positions))
    al_positions = Math.DMatrix(3, len(saved_mean_positions))
    i = 0

    print '- Aligning frames...'
   
    while i < len(saved_mean_positions):
        pos = saved_mean_positions[i]

        al_ref_positions.setElement(0, i, pos[0])
        al_ref_positions.setElement(1, i, pos[1])
        al_ref_positions.setElement(2, i, pos[2])
        i += 1

    i = 0
    xform = Math.Matrix4D()

    while i < num_confs:
        j = 0

        while j < len(saved_mean_positions):
            pos = saved_res_positions[j][i]

            al_positions.setElement(0, j, pos[0])
            al_positions.setElement(1, j, pos[1])
            al_positions.setElement(2, j, pos[2])
            j += 1

        if not alignment.align(al_positions, al_ref_positions):
            print '!! Could not align frame', i

        else:
            #print 'Aligning frame', i
            xform.assign(alignment.getTransform())
            Chem.transformConformation(mol, i, xform)

        i += 1

    #Util.setResidueSeqNumbers(mol, old_seq_nos)

    if not Util.saveCDFMolecule(sys.argv[3], mol):
        print '!! Could not write file'
        sys.exit(2)


if __name__ == '__main__':
    process()
