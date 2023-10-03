#!/usr/bin/env python
#This module include functions calculate distance of two atoms under periodical boundary condition
import os
import numpy as np
from numpy.linalg import norm, inv

import copy
import itertools

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


def find_center_group(vecR, list_pos, list_center_species, list_neighbour_species, bondlength_max=2.5):
    '''
    Find all connected atoms in list_neighbour_species to species in list_site_species
    Generally find XOn unit in oxides
    Note one atom could be in more than one group

    :return: Two lists, one of atoms in the group : [ [atomcenter, atom1, atom2, ...], [...], []], another of dist [ [ center-1, center-2 ,..], [] ,[]..]
    '''
    list_v = [x for x in list_pos if x["species"] in list_center_species]
    list_o = [(i, x) for i, x in enumerate(list_pos) if x["species"] in list_neighbour_species]

    list_group = []
    list_dist = []
    for i in range(len(list_v)):
        atom = list_v[i]
        group = [atom]
        group_dist = []
        for j, atom2 in list_o:
            dist, shift, v = find_distance(atom["pos"], atom2["pos"], vecR)
            if (dist < bondlength_max):
                group.append(j)
                group_dist.append(dist)
        list_group.append(group)
        list_dist.append(group_dist)

    return list_group, list_dist




