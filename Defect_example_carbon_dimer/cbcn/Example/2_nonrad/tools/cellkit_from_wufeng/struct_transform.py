#!/usr/bin/env python
#Transform a structure based on matrix and shift
from argparse import ArgumentParser
import os
import numpy as np
from numpy.linalg import norm, inv
import copy
import sys

from io_package import write_cell_and_pos_auto, read_cell_and_pos_auto

parser = ArgumentParser(description="Transform atoms based on rotation matrix and shift vector. New atom positions would be R.atom + S. Cell vectors are kept as before. Note if cell is assumed to be R.C, each row in C is a lattice vector, the structure is not changed")
parser.add_argument('-i', dest="input", help='POSCAR type file to convert')
parser.add_argument('-o', dest="output", help='Filename to output')
parser.add_argument('-r', dest="rotate", help="Rotation matrix, 9 numbers as m[1,1], m[1,2]. New position for atom (a,b,c) is a' = m[1,1] * a + m[1,2] * b + m[1,3] * c")
parser.add_argument('-s', dest="shift" , help="Shift vector, 3 numbers")
args = parser.parse_args()

rot = np.asarray(map(float, args.rotate.split())).reshape((3,3))
shift = np.asarray(map(float, args.shift.split()))

(vecR, list_atom_base), package  = read_cell_and_pos_auto(args.input)
list_atom = copy.deepcopy(list_atom_base)
for x in list_atom:
    x["pos"] = np.dot(rot, x["pos"]) + shift

write_cell_and_pos_auto(package, args.output, vecR, list_atom)

