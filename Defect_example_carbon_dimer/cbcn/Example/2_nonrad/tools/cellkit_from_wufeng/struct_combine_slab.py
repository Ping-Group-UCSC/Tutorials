#!/usr/bin/env python
#Combine multiple layered structure
from argparse import ArgumentParser
import os
import numpy as np
from numpy.linalg import norm, inv
import copy
import sys

from io_package import write_cell_and_pos_auto, read_cell_and_pos_auto

parser = ArgumentParser(description="Combine two slab structures together to create a new slab structure.")
parser.add_argument('-i1', dest="input1", required=True, help='First layerd structure (placed close to 0)')
parser.add_argument('-i2', dest="input2", required=True, help='Second layered structure (placed between first top and 1')
parser.add_argument('-o', dest="output", help='Filename to output')
parser.add_argument('--ab', dest="vec_source", required=True, type=int, help='Can be 1 or 2, Use input1 or input2 ab axis as the final ab axis. Even from different lattice vectors, the a/b crystal coordinate of the other will be copied as is to final output.')
parser.add_argument('--distance', dest="distance", type=float, required=True, help='The distance between structure 1 top layer and structure 2 bottom layer. Unit in Ang.')
parser.add_argument('--vac', dest="vacuum", type=float, required=True, help='The vacuum to separaete the whole 1+2 slab. Unit in Ang.')
parser.add_argument('--forceperpc', default=True, dest="b_forceperpc", help='The c-axis must be exactly z and the a/b must be in xy plane, 0.0001 style error will be enforced to be zero. This cannot be turned off now.')

args = parser.parse_args()

def check_ortho_z(vecR):
    """
    Check if the system ab-axis are in plane and c-axis are along z

    :return: the vecR with ab-axis z-component and c-axis ab component cleared
    """
    tol = 1e-2
    b_notortho = (np.abs(vecR[2,0:2]) > tol).any() or (np.abs(vecR[0:2,2]) > tol).any()
    if (b_notortho):
        raise ValueError("ab-axis must be in xy plane and c-axis must be along z-axis")

    vecR[2,0:2] = 0
    vecR[0:2,2] = 0

    return vecR

def shift_layer_to_0(ar_z_old):
    '''
    Shift the layer out of bound back
    Put it in a region at 0 to x
    '''
#Convert to positive values first 
    ar_z = np.mod(ar_z_old, 1)

#Test where is the vacuum (largest seperation)
    ar_z.sort()
    diff = (ar_z[1:] - ar_z[:-1]).tolist()
    diff.append(ar_z[0] - ar_z[-1]+1)
    diff = np.array(diff)
    ix = diff.argmax() + 1
    if (ix == len(diff)):
        ix = 0
    shift = -ar_z[ix]
    return np.mod(ar_z_old + shift, 1)

    

def find_layers(ar_z, tol=0.05):
    """
    Find layers from a list of positions with specific tol
    """
    ar_ix = np.argsort(ar_z)
    list_layer = []
    layer = []
    last_z = -100
    for ix in ar_ix:
        if (ar_z[ix] > last_z + tol):
            list_layer.append(layer)
            layer = [ix]
        else:
            layer.append(ix)

    list_layer.append(layer)
    del list_layer[0]
    return list_layer


(vecR1, list_atom1), package  = read_cell_and_pos_auto(args.input1)
(vecR2, list_atom2), package  = read_cell_and_pos_auto(args.input2)

vecR1 = check_ortho_z(vecR1)
vecR2 = check_ortho_z(vecR2)

if (args.vec_source == 1):
    vecR = vecR1.copy()
elif (args.vec_source == 2):
    vecR = vecR2.copy()
else:
    raise ValueError("Unrecognized option %s" % args.vec_source)

#Shift layers to 0 and convert to Cartesian coord
ar_z1 = vecR1[2,2] * shift_layer_to_0(np.array([x["pos"][2] for x in list_atom1]))
ar_z2 = vecR2[2,2] * shift_layer_to_0(np.array([x["pos"][2] for x in list_atom2]))

#Find layers
#layers1 = find_layers(ar_z1)
#layers2 = find_layers(ar_z2)

#Compute the cell size
#list_thick = [ lz[ll[-1][-1]] - lz[ll[0][0]] for lz, ll in zip([ar_z1, ar_z2], [layers1, layers2])]
#print(list_thick)
shift2 = ar_z1.max() + args.distance
thick = ar_z1.max() + args.vacuum + ar_z2.max() + args.distance
vecR[2,2] = thick

print("New cell z : %.3f" % thick)

#Shift z
#Rule : 
#Structure1 : smallest to 0, 
#Structure2 : smallest to Structure1 largest + distance
list_atom = []
for ix, atom in enumerate(list_atom1):
    atom["pos"][2] = ar_z1[ix] / thick
    list_atom.append(atom)

for ix, atom in enumerate(list_atom2):
    atom["pos"][2] = (ar_z2[ix] + shift2) / thick
    list_atom.append(atom)

write_cell_and_pos_auto("vasp", args.output, vecR, list_atom)


