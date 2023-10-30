#!/usr/bin/env python
#Mix two structures from QE or vasp results
from argparse import ArgumentParser
import os
import numpy as np
from string import digits, ascii_letters
import re
import sys
import shutil

from io_package import read_cell_and_pos_auto, write_cell_and_pos_auto


def mix_structure(vecR, list_pos1, list_pos2, ratio, head, tail, filename_out, package):
    '''
    Mix two structures
    If ratio is 0, as pos1, 1 as pos2
    Otherwise use linear extrapolation

    Also check atom names (except number suffix, which can be different as polaron positions)

    If two atoms are too far away, treated as cross-boundary ; will mix the minimum difference (i.e. 0.05 and 0.95 will be mixed the same as 0.05 and -0.05)

    '''
    if (len(list_pos1) != len(list_pos2)):
        raise ValueError("Different number of atoms")


    l3 = []
    suffix1 = None
    suffix2 = None
    for p1, p2 in zip(list_pos1, list_pos2):
        element1 = p1["species"].translate(None, digits)
        element2 = p2["species"].translate(None, digits)
        if (element1 != element2):
            raise ValueError("Different atom found : %s %s" % (p1["species"], p2["species"]))
#Match atom positions to the nearest position (0.01 and 0.99 => 0.01 and -0.01)
        vec2 = p2["pos"].copy()
        for i in range(len(vec2)):
            if (abs(vec2[i] - p1["pos"][i]) > 0.5):
                vec2[i] -= round(vec2[i] - p1["pos"][i], 0)

        pos3 = p1["pos"] * (1-ratio) + vec2 * ratio
        name = p1["species"]
        l3.append({"species": name, "pos" : pos3})

    if package == "qe":
        with open(filename_out, 'w') as f:
            for line in head:
                f.write(line)
            for atom in l3:
                f.write("%-4s %.10f %.10f %.10f\n" % (atom["species"], atom["pos"][0], atom["pos"][1], atom["pos"][2]))
            for line in tail:
                f.write(line)
    elif package == "vasp":
        write_cell_and_pos_auto("vasp", filename_out, vecR, l3)

    return

parser = ArgumentParser(description="Linear interpolate two structures. A series from structures will be A template.in as QE header must present to generate QE output files. If occ.in presents it will be appended as QE OCCUPATIONS card")
parser.add_argument('filenames', nargs=2, help='POSCAR type or QE in/out file to convert. For QE, specify "xx" for an input "xx.in" and an output "xx.out", the relaxed structure in .out will be used.')
parser.add_argument('-n', dest="n_ratio", type=int, default=None, help='The number of structures generated  between two structures. 3 means mixing ratio 0.25, 0.5 and 0.75. Only one of -r and -n must be specified.')
parser.add_argument('-r', dest="list_ratio", default=None, help='A list of ratio to generate structures, like "0.25 0.50 0.75". Only one of -r and -n must be specified.')
parser.add_argument('--package', dest="package", default="qe", help='Possible values: qe/vasp, for QE write scf.in and for VASP write POSCAR')
args = parser.parse_args()

if (args.n_ratio is None and args.list_ratio is None):
    raise ValueError("One of -r and -n must be specified.")

if (args.list_ratio is not None):
    list_ratio = map(float, args.list_ratio.split())
else:
    list_ratio = np.linspace(0, 1, args.n_ratio + 1, endpoint=False)[1:]


if args.package == "qe":
    filename_scf = "template.in" 
    with open(filename_scf, 'r') as f:
       head = f.readlines() 

    filename_occ = "occ.in" 
    if (os.path.exists(filename_occ)):
        with open(filename_occ, 'r') as f:
           tail = f.readlines() 
    else:
        tail = []
else:
    head = []
    tail = []

(vecR, list_pos1), prog = read_cell_and_pos_auto(args.filenames[0])
(vecR2, list_pos2), prog = read_cell_and_pos_auto(args.filenames[1])

filename_new = {"qe" : "scf.in" , "vasp" : "POSCAR"}[args.package]

if ((vecR != vecR2).any()):
    print(vecR, vecR2)
    raise ValueError("Inconsistent lattice vectors!")
for ratio in list_ratio:
    dirname = "ratio-%.4f" % ratio
    if (not os.path.exists(dirname)):
        os.makedirs(dirname)
    if (os.path.exists(os.path.join(dirname, filename_new))):
        continue
    mix_structure(vecR, list_pos1, list_pos2, ratio, head, tail, os.path.join(dirname, filename_new), args.package)

