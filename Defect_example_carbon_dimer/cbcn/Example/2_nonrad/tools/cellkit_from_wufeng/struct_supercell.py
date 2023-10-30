#!/usr/bin/env python
#Transform a structure based on matrix and shift
from math import ceil, floor
from argparse import ArgumentParser
import os
import numpy as np
from numpy.linalg import norm, inv
import copy
import sys
import itertools

from io_package import write_cell_and_pos_auto, read_cell_and_pos_auto

parser = ArgumentParser(description="Transform cell parameters based on rotation matrix and shift vector. New atom positions would be R.atom + S. Matrix can be a supercell, new atoms will be added")
parser.add_argument('-i', dest="input", help='POSCAR type file / QE prefix to convert')
parser.add_argument('-o', dest="output", help='Filename to output')
parser.add_argument('-r', dest="rotate", help="Rotation matrix, 9 numbers as m[1,1], m[1,2]. New position for atom (a,b,c) is a' = m[1,1] * a + m[1,2] * b + m[1,3] * c")
parser.add_argument('-s', dest="shift" , help="Shift vector, 3 numbers")
args = parser.parse_args()

rot = np.asarray(map(float, args.rotate.split())).reshape((3,3))
shift = np.asarray(map(float, args.shift.split()))

(vecR, list_atom_base), package  = read_cell_and_pos_auto(args.input)
list_atom = copy.deepcopy(list_atom_base)
print(vecR)

vecRI = inv(vecR)
vecR2 = np.dot(rot, vecR.T).T
print(vecR2)
vecR2I = inv(vecR2)

#Find the maximum of cells we must search
#Find positions of 8 verticies
list_vunit = []
for vertex_abc in (
        (0,0,0), (0,1,0), (0,0,1), (1,0,0), (1,1,0), (1,0,1), (0,1,1), (1,1,1)
        ):
#Cartesian
    v1 = np.dot(vecR2, np.asarray(vertex_abc).reshape(3,1))
#In unit of unit cell
    list_vunit.append(np.dot(vecRI, v1))

ar_vunit = np.asarray(list_vunit)
list_r = []
for i in range(3):
    list_r.append(range(int(floor(ar_vunit[:,i].min()-1)),int(ceil(ar_vunit[:,i].max()+1+1))))

tol = 1e-4

list_atom2 = []
for vcell in itertools.product(*list_r):
    for x in list_atom:
        pos = np.dot(vecR, np.asarray(vcell) + x["pos"] + shift)
        pos2 = np.dot(vecR2I, pos)
#       print(pos2)
        if ((pos2 > -tol).all() and (pos2 < 1 + tol).all()):
            atom2 = copy.deepcopy(x)
            atom2["pos"] = pos2
            list_atom2.append(atom2)
#           print(pos2)

def check_duplicate(atoms, tol=1e-4):
    '''
    Check if atoms are too close, later one will be removed

    tol in crystal coordinate
    '''
    l1 = []
    set_remove = set()
    for i0, atom1 in enumerate(atoms):
        if (i0 in set_remove):
            continue
        for i in range(len(atoms)):
            if (i0 == i):
                continue
            diff = np.abs(atom1["pos"] - atoms[i]["pos"])
            diff = diff - np.floor_divide(diff, 1)
            if ( (diff < tol).all()):
                set_remove.add(i)

    for i in range(len(atoms)):
        if (i not in set_remove):
            l1.append(atoms[i])

    print("Remove %i duplicated atoms, %i left" % (len(set_remove), len(l1)))

    return l1
        

#Duplicate check
list_atom2 = check_duplicate(list_atom2)


write_cell_and_pos_auto(package, args.output, vecR2, list_atom2)

