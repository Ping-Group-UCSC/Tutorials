#!/usr/bin/env python
#Given two POSCAR, convert the atoms in the first one to the same as second one
from argparse import ArgumentParser
import os
import sys
import numpy as np
from numpy.linalg import norm, inv

import copy
from io_package import read_cell_and_pos_poscar, write_poscar


parser = ArgumentParser(description="Re-order atoms in first structure according to second structure")
parser.add_argument('-i', dest="input", help='POSCAR type file to convert')
parser.add_argument('-r', dest="ref", help='POSCAR type file as reference')
parser.add_argument('-o', dest="output", help='Filename to output')
args = parser.parse_args()

#Main program
if (not os.path.exists(args.input) or not os.path.exists(args.ref)):
    print("Cannot find files")
    sys.exit(1)

vecR, list_atom_base  = read_cell_and_pos_poscar(args.input)
vecR2, list_atom2  = read_cell_and_pos_poscar(args.ref)

if ((np.abs(vecR - vecR2) > 1e-4).any()):
    print("Different crystals")
    sys.exit(2)


list_atom = []

for atom2 in list_atom2:
    pos2 = atom2["pos"]
    b_add = False
    for atom in list_atom_base:
        pos1 = atom["pos"]
        delta = pos1 - pos2
        for i in range(len(delta)):
            delta[i] -= round(delta[i])
        delta = np.dot(vecR, delta)
        if ((np.abs(delta) < 0.2).all()):
            list_atom.append(atom)
            b_add = True
            break
    if (not b_add):
        raise ValueError("Did not find any atom match for %s" % pos2)

#Calculate difference
tot = 0
for atom1, atom2 in zip(list_atom, list_atom2):
    delta = atom1["pos"] - atom2["pos"]
    for i in range(len(delta)):
        delta[i] -= round(delta[i])
    tot += norm(np.dot(vecR, delta))


print("Total difference: %.4f Ang" % tot)
write_poscar(args.output, vecR, list_atom)

