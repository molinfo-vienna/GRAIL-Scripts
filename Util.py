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


import math
import CDPL.Chem as Chem
import CDPL.Biomol as Biomol
import CDPL.Math as Math


def getAlphaCAtom(res):
    for atom in res.atoms:
        if Biomol.getResidueAtomName(atom) == 'CA':
            return [atom]

    return []

def getBackboneAtoms(res):
    atoms = []

    for atom in res.atoms:
        if  Biomol.isPDBBackboneAtom(atom):
            atoms.append(atom)

    return atoms

def getAllHeavyAtoms(res):
    atoms = []

    for atom in res.atoms:
        if  Chem.getType(atom) != Chem.AtomType.H:
            atoms.append(atom)

    return atoms

def calcAtomSetCentroid(atoms, conf_idx):
    if len(atoms) == 1:
        return Chem.getConformer3DCoordinates(atoms[0], conf_idx)

    ctr = Math.Vector3D()

    for atom in atoms:
        ctr += Chem.getConformer3DCoordinates(atom, conf_idx)

    ctr /= len(atoms)
    return ctr

def calcMeanPos(positions):
    mean_pos = Math.Vector3D()

    for pos in positions:
        mean_pos += pos

    mean_pos /= len(positions)
    return mean_pos


def calcMSF(positions, mean_pos):
    msf = 0.0
    tmp = Math.Vector3D()

    for pos in positions:
        tmp.assign(pos)
        tmp -= mean_pos
        msf += Math.innerProd(tmp, tmp)

    msf = msf / len(positions)
    return msf

def calcRMSF(positions, mean_pos):
    msf = calcMSF(positions, mean_pos)

    return math.sqrt(msf)

def calcMaxFluct(positions, mean_pos):
    max_fluct = 0.0
    tmp = Math.Vector3D()

    for pos in positions:
        tmp.assign(pos)
        tmp -= mean_pos
        dev = Math.norm2(tmp)

        if dev > max_fluct:
            max_fluct = dev

    return max_fluct

def getResidueID(res):
    if Biomol.getChainID(res) != ' ':
        return Biomol.getResidueCode(res) + '_' + str(Biomol.getResidueSequenceNumber(res)) + '_' + Biomol.getChainID(res)

    return Biomol.getResidueCode(res) + '_' + str(Biomol.getResidueSequenceNumber(res))

def loadCDFMolecule(fname):
    mol = Chem.BasicMolecule()
    cdf_reader = Chem.FileCDFMoleculeReader(fname)
      
    if not cdf_reader.read(mol):
        return None

    Chem.calcImplicitHydrogenCounts(mol, False)
    Chem.perceiveHybridizationStates(mol, False)        
    Chem.setAtomSymbolsFromTypes(mol, False)        
    Chem.perceiveSSSR(mol, False)
    Chem.setRingFlags(mol, False)
    Chem.setAromaticityFlags(mol, False)

    return mol

def saveCDFMolecule(fname, mol):
    cdf_writer = Chem.FileCDFMolecularGraphWriter(fname)
      
    if not cdf_writer.write(mol):
        return None

    return mol

def readLines(fname):
    f = file(fname, 'r')
    lines = []

    for l in f.readlines():
        l = l.strip()

        if l: lines.append(l.strip())

    return lines

def toIntegerList(str_list):
    ints = []

    for l in str_list:
        ints.append(int(l))

    return ints

def fixResidueSeqNumbers(mol):
    curr_seq_no = 0
    res_start = 0
    res_end = 0
    num_atoms = mol.getNumAtoms()
    old_res_seq_nos = []
    old_to_new_seq_no_map = {}

    for atom in mol.atoms:
         old_res_seq_nos.append(Biomol.getResidueSequenceNumber(atom))

    while res_start < num_atoms:
        while res_end < num_atoms and (old_res_seq_nos[res_start] == old_res_seq_nos[res_end]):
            res_end += 1

        if not old_res_seq_nos[res_start] in old_to_new_seq_no_map:
            old_to_new_seq_no_map[old_res_seq_nos[res_start]] = curr_seq_no

        while res_start < res_end:
            Biomol.setResidueSequenceNumber(mol.getAtom(res_start), curr_seq_no)
            res_start += 1

        curr_seq_no += 1
        
    return old_res_seq_nos, old_to_new_seq_no_map

def setResidueSeqNumbers(mol, seq_nos):
    i = 0

    while i < len(seq_nos):
        Biomol.setResidueSequenceNumber(mol.getAtom(i), seq_nos[i])
        i += 1

def filterResidues(residues, res_subset):
    if len(res_subset) == 0:
        return

    i = 0
    
    while i < len(residues):
        seq_no = Biomol.getResidueSequenceNumber(residues[i])

        if seq_no in res_subset:
            i += 1
            continue

        del residues[i]

def writePVDHeader(pvd_file):
    pvd_file.write('<?xml version="1.0"?>\n<VTKFile type="Collection" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">\n<Collection>\n')

def writePVDFooter(pvd_file):
    pvd_file.write('</Collection>\n</VTKFile>\n')
    pvd_file.close()

def writePVDEntry(pvd_file, frame_no, fname, suffix):
    pvd_file.write('<DataSet timestep="' + str(frame_no) + '" group="" part="0" file="' + fname + '.' + suffix + '"/>\n')
