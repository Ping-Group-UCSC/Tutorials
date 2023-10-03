#!/usr/bin/env python
#Transform a structure based on matrix and shift
from argparse import ArgumentParser
import os
import numpy as np
from numpy.linalg import norm, inv
import copy
import sys
from io_package import write_pos_keep_qe, read_pos_only_qe

parser = ArgumentParser(description="Expand the tructrure with a given fraction, calculate Cartesian coordinates after. Note if the unit is crystal we do not need that")
parser.add_argument('-i', dest="input", required=True, help='QE input file to convert')
parser.add_argument('-o', dest="output", required=True, help='Filename to output')
parser.add_argument('-s', dest="scale", type=float, default=None, required=True, help="The coefficient to expand, in the format of 0.99, 1.01")
args = parser.parse_args()


vecR, list_atom, unit, lines_before, lines_after = read_pos_only_qe(args.input)

if (unit == "crystal"):
    print("Crystal coordinates, should not be expanded")

for atom in list_atom:
    atom["pos"] *= args.scale

write_pos_keep_qe(args.output, vecR, list_atom, unit, lines_before, lines_after)
print("Warning: the lattice is not changed, please modify it manully")

