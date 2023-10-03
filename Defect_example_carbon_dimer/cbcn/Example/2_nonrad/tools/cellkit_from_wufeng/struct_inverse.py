#!/usr/bin/env python
#Inverse the POSCAR structure
from argparse import ArgumentParser
import os
import numpy as np
from numpy.linalg import norm, inv

import copy
import itertools
from collections import Counter

from io_package import write_cell_and_pos_auto, read_cell_and_pos_auto

parser = ArgumentParser(description="Inverse the structure")
parser.add_argument('-i', dest="input", help='POSCAR type or QE in/out file to convert')
parser.add_argument('-o', dest="output", help='Filename to output; same as what read')
args = parser.parse_args()

#Main program
(vecR, list_atom_base), package  = read_cell_and_pos_auto(args.input)
list_atom = copy.deepcopy(list_atom_base)
for x in list_atom:
    x["pos"] = 1 - x["pos"]

write_cell_and_pos_auto(package, args.output, vecR, list_atom)

