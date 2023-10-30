#!/usr/bin/env python
#Given two structures, convert the atoms in the first one to the same as second one
from argparse import ArgumentParser
import os
import sys
import numpy as np
from numpy.linalg import norm, inv

import copy
from io_package import write_cell_and_pos_auto, read_cell_and_pos_auto


parser = ArgumentParser(description="Re-order atoms in first structure according to second structure")
parser.add_argument('-i', dest="input", help='POSCAR type file to convert')
parser.add_argument('-r', dest="ref", help='POSCAR type file as reference')
parser.add_argument('-o', dest="output", help='Filename to output')
parser.add_argument('-t', dest="threshold", type=float, default=0.2, help='Threshold of same distance between atoms')
parser.add_argument('-c', dest="b_continue_notfound", default=False, action="store_true", help='Do not stop if one atom cannot be matched within threshold. Will take the nearest instead.')
args = parser.parse_args()

(vecR, list_atom_base), package1  = read_cell_and_pos_auto(args.input)
(vecR2, list_atom2), package2  = read_cell_and_pos_auto(args.ref)

if ((np.abs(vecR - vecR2) > 1e-1).any()):
    print("Error: Different crystals, diff > 1e-1")
    sys.exit(2)
if ((np.abs(vecR - vecR2) > 1e-4).any()):
    print("Warning: Different crystals, diff > 1e-4")
if len(list_atom_base) != len(list_atom2):
    print("Different number of atoms")


list_atom = []

#Old code : no duplicate check
#for atom2 in list_atom2:
#    pos2 = atom2["pos"]
#    b_add = False
#    for atom in list_atom_base:
#        pos1 = atom["pos"]
#        delta = pos1 - pos2
#        for i in range(len(delta)):
#            delta[i] -= round(delta[i])
#        delta = np.dot(vecR, delta)
#        if ((np.abs(delta) < args.threshold).all()):
#            list_atom.append(atom)
#            b_add = True
#            break
#    if (not b_add):
#        raise ValueError("Did not find any atom match for %s" % pos2)
set_matched = set(range(len(list_atom_base)))
for atom2 in list_atom2:
    pos2 = atom2["pos"]
    list_dist = []
    for ix_atom in set_matched:
        atom = list_atom_base[ix_atom]
        if (atom2["species"] != atom["species"]):
            continue
        pos1 = atom["pos"]
        delta = pos1 - pos2
#Shift back to the same cell; always assume the atom does not move > 0.5
        for i in range(len(delta)):
            delta[i] -= round(delta[i])
        delta = np.dot(vecR, delta)
        list_dist.append([ix_atom, norm(delta)])

    if (len(list_dist) == 0):
        raise ValueError("Did not find any atom match for %s" % pos2)

    list_dist.sort(key=lambda x:x[1])
    ix_atom, dist_min = list_dist[0]

    print("Atom %.3f %.3f %.3f : Dist %f -> Ref %i" % (
        atom2["pos"][0], atom2["pos"][1], atom2["pos"][2], dist_min, ix_atom+1))
    if (dist_min < args.threshold or args.b_continue_notfound):
        list_atom.append(list_atom_base[ix_atom])
        set_matched.remove(ix_atom)
    else:
        raise ValueError("No match within threshold: %s" % pos2)

#Calculate difference
tot = 0
for atom1, atom2 in zip(list_atom, list_atom2):
    delta = atom1["pos"] - atom2["pos"]
    for i in range(len(delta)):
        delta[i] -= round(delta[i])
    tot += norm(np.dot(vecR, delta))


print("Total difference: %.4f Ang" % tot)
write_cell_and_pos_auto(package1, args.output, vecR, list_atom)

