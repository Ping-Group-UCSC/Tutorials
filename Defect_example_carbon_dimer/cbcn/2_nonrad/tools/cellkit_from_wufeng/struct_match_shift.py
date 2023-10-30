#!/usr/bin/env python
#Same as struct_match, match two atom list, and also shift atoms to make sure they are close in coordinates (for example, 0.0 and 0.99 will be shift to 0.0 and -0.01)
from argparse import ArgumentParser
import os
import sys
import numpy as np
from numpy.linalg import norm, inv

import copy
import itertools
from collections import Counter
from io_package import write_cell_and_pos_auto, read_cell_and_pos_auto

parser = ArgumentParser(description="Re-order atoms in first structure according to second structure")
parser.add_argument('-i', dest="input", help='POSCAR-type or QE prefix to convert (read prefix.in and optimized .out)')
parser.add_argument('-r', dest="ref", help='POSCAR-type or QE prefix as reference (read prefix.in and optimized .out')
parser.add_argument('-o', dest="output", help='Filename to output, using same format as -i')
args = parser.parse_args()

#Main program
#if (not os.path.exists(args.input) or not os.path.exists(args.ref)):
#    print("Cannot find files")
#    sys.exit(1)

(vecR, list_atom_base), package1 = read_cell_and_pos_auto(args.input)
(vecR2, list_atom2), package2  = read_cell_and_pos_auto(args.ref)

if ((np.abs(vecR - vecR2) > 1e-4).any()):
    print("Different crystals")
    sys.exit(2)

if (len(list_atom_base) != len(list_atom2)):
    print("Different number of atoms")
print("Map %i atoms: " % len(list_atom2))

list_atom = []

for atom2 in list_atom2:
    pos2 = atom2["pos"]
    b_add = False
    for atom in list_atom_base:
        pos1 = atom["pos"]
        delta = pos1 - pos2
        for i in range(len(delta)):
            delta[i] -= round(delta[i])
        delta_cart = np.dot(vecR, delta)
        if ((np.abs(delta_cart) < 0.2).all()):
            atom["pos"]  = pos2 + delta
            list_atom.append(atom)
            b_add = True
#           print("Map %s to %s" % (pos2, atom["pos"]))
            break
    if (not b_add):
        raise ValueError("Did not find any atom match for %s" % pos2)

#Calculate difference
list_diff = []
for atom1, atom2 in zip(list_atom, list_atom2):
    delta = atom1["pos"] - atom2["pos"]
#   for i in range(len(delta)):
#       delta[i] -= round(delta[i])
    list_diff.append(norm(np.dot(vecR, delta)))

list_diff.sort(reverse=True)
    
tot = sum(list_diff)

#print(list_diff)

print("Total difference: %.4f Ang" % tot)
write_cell_and_pos_auto(package1, args.output, vecR, list_atom)

