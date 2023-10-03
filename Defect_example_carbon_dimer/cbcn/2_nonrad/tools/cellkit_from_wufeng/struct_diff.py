#!/usr/bin/env python
#Given two POSCAR, measure the difference between all atoms and find largest differences
from argparse import ArgumentParser
import os
import sys
import numpy as np
from numpy.linalg import norm, inv

import copy
import itertools
from collections import Counter

cell_neighbor = [np.asarray(x) for x in itertools.product((-1,0,1),(-1,0,1),(-1,0,1))]
def find_distance(v1, v2, vecR):
    '''
    Find nearest distance between two atoms, couting all cells
    '''
    list_delta = []
    list_v = []
    for cellshift in cell_neighbor:
        list_v.append(np.dot(vecR , (cellshift + v1 - v2)))
        list_delta.append(norm(list_v[-1]))
    delta = min(list_delta)
    ix = list_delta.index(delta)
    return (delta, cell_neighbor[ix], list_v[ix])


def read_cell_and_pos_qe(prefix):
    '''
    Read the cell and positions from .in file
    Note all numbers after species (like V1, Fe1) will be removed :  to (V,Fe)
    '''
#Read cell parameters
    with open(prefix + ".in", 'r') as f:
        lines = f.readlines()

    list_pos = []
    for i, line in enumerate(lines):
        if ("nat" in line):
            ar = line.split(",")
            for st2 in ar:
                if ("nat" in st2):
                    ar2 = st2.split("=")
                    nat = int(ar2[-1])
        elif ("CELL_PARAMETER" in line):
            if ("ang" not in line):
                raise ValueError("Only angstrom unit is supported")
            vecR = np.asarray([[float(x) for x in line.split()] for line in lines[i+1:i+4]]).T
        elif ("ATOMIC_POSITIONS" in line):
            if ("crystal" not in line):
                raise ValueError("Only crystal coordinate is supported")
            for line2 in lines[i+1:i+nat+1]:
                ar = line2.strip().split()
                list_pos.append({"species": filter(lambda x:x.isalpha(), ar[0]), "pos" : np.asarray([float(x) for x in ar[1:]])})

    return vecR, list_pos

def read_cell_and_pos_poscar(filename):
    '''
    Read the cell and positions
    '''
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    vecR = np.asarray([[float(x) for x in line.split()] for line in lines[2:5]]).T
    atoms = lines[5].split()
    natoms = [int(x) for x in lines[6].split()]
    list_pos = []
    i = 8
    for atom, n in zip(atoms, natoms):
        for j in range(n):
            list_pos.append({"species": atom, "pos" : np.asarray([float(x) for x in lines[i].split()])})
            i += 1

    return vecR, list_pos


list_site_species = ["Fe"]

def is_polaron_site(atom):
    '''
    Determine whether a site belongs to a polaron
    '''
    return atom["species"] in list_site_species


def find_center_group(vecR, list_pos, list_center_species, list_neighbour_species):
    '''
    Find all connected atoms in list_neighbour_species to species in list_site_species
    Generally find XOn unit in oxides
    Note one atom could be in more than one group
    '''
    list_v = [x for x in list_pos if x["species"] in list_center_species]
    list_o = [(i, x) for i, x in enumerate(list_pos) if x["species"] in list_neighbour_species]

    list_group = []
    for i in range(len(list_v)):
        atom = list_v[i]
        group = [atom]
        for j, atom2 in list_o:
            dist, shift, v = find_distance(atom["pos"], atom2["pos"], vecR)
            if (dist < 2.5):
                group.append(j)
        list_group.append(group)

    return list_group

def write_poscar(filename, vecR, list_pos):
    '''
    Write to POSCAR file
    '''
    with open(filename, 'w') as f:
        f.write("Polaron structure\n1.0\n")
        for i in range(3):
            f.write("%.10f %.10f %.10f\n" % tuple(vecR[:,i]))

        list_species = [x["species"] for x in list_pos]
        dic_species = Counter(list_species)
        for species, val in dic_species.items():
            f.write("%4s " % species)
        f.write("\n")
        for species, val in dic_species.items():
            f.write("%4i " % val)
        f.write("\n")
        f.write("Direct\n")
        for species, val in dic_species.items():
            for atom in list_pos:
                if (atom["species"] == species):
                    f.write("%.10f %.10f %.10f\n" % tuple(atom["pos"]))

parser = ArgumentParser(description="Re-order atoms in first structure according to second structure")
parser.add_argument('-i', dest="input", help='POSCAR type file to convert')
parser.add_argument('-r', dest="ref", help='POSCAR type file as reference')
#parser.add_argument('-o', dest="output", help='Filename to output')
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
list_delta = []

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
            list_delta.append((len(list_atom), norm(delta), delta))
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


list_delta.sort(key=lambda x:-x[1])
print("Largest different atoms: ")
for i in range(5):
    print("%4i %.4f %s" % list_delta[i])
