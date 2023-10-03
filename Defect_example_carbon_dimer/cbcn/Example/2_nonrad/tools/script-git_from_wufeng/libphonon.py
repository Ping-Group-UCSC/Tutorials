#!/usr/bin/env python
#Require poscar to have species line
import os
import itertools
import numpy as np
from math import *
from constant import *

elements = ['H' , 'He', 'Li', 'Be', 'B',  'C',  'N',  'O' , 'F' , 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P' , 'S',  'Cl', 'Ar', 'K' , 'Ca',\
            'Sc', 'Ti', 'V' , 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',\
            'Rb', 'Sr', 'Y' , 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I' , 'Xe',\
            'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu',\
                  'Hf', 'Ta', 'W' , 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn',\
            'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U' , 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr',\
                  'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mg', 'Ds', 'Rg']

dic_element = dict(zip(elements, range(1, len(elements)+1)))

def plot_phonon_mode_multi(dir1, vecR, list_pos, list_eigenvector_phonon, **kwargs):
    '''
    Plot all phonon normal mode together in a folder
    '''
#   dir1 = "phonon-axsf"
    if (not os.path.exists(dir1)):
        os.mkdir(dir1)

    for ix in range(len(list_eigenvector_phonon)):
        plot_phonon_mode(os.path.join(dir1, "%s" % (ix+1)), 
                vecR, list_pos, list_eigenvector_phonon[ix], **kwargs)

    return

def plot_phonon_mode(prefix, vecR, list_pos, eigenvector_phonon, unit_pos="crystal", unit_eigen="angstrom", scale=1, list_dataformat=["axsf", "xsf"]):
    '''
    Plot phonon normal mode into axsf for visualization

    All input variable must be in unit of angstrom, unless marked

    :param prefix: The prefix, add ".***" with data format as the filename

    '''
#Number of snapshots in animation
#Note the last will be a copy of the first and not included in this number
    n_step = 23

    v_st =  "".join(["%.7f %.7f %.7f\n" % tuple(vecR[:,i]) for i in range(3)])
    list_atom = [dic_element[x["species"]] for x in list_pos]
    n_atom = len(list_atom)
    ar_pos = np.asarray([x["pos"] for x in list_pos])
    ar_force = eigenvector_phonon * scale

    if (unit_pos == "crystal"):
        ar_pos = np.dot(vecR, ar_pos.T).T
    elif (unit_pos == "bohr"):
        ar_pos = ar_pos * Bohr2Ang

    if (unit_eigen == "crystal"):
        ar_force = np.dot(vecR, ar_force.T).T
    elif (unit_eigen == "bohr"):
        ar_force = ar_force * Bohr2Ang

#Save for different format
    for dataformat in list_dataformat:
        filename = "%s.%s" % (prefix, dataformat)
        if (dataformat == "axsf"):
            with open(filename, 'w') as f:
                f.write("ANIMSTEPS %i\nCRYSTAL\nPRIMVEC\n%sCONVVEC\n%s" % (n_step + 1, v_st, v_st))
                for ix_step in range(n_step+1):
                    f.write("PRIMCOORD %i\n%i %i\n" % (ix_step + 1, n_atom, 1))
                    shift = sin(2*pi*ix_step / n_step + pi/2)
                    ar_pos2 = ar_pos + shift * ar_force
                    for ix_atom in range(len(list_atom)):
                        f.write("%i %s\n" % (list_atom[ix_atom],
                            " ".join(["%.6f" % x for x in ar_pos2[ix_atom,:]])))
        elif (dataformat == "xsf"):
            ar_posforce = np.concatenate((ar_pos, ar_force), axis=1)
            with open(filename, 'w') as f:
                f.write("CRYSTAL\nPRIMVEC\n%sCONVVEC\n%s" % (v_st, v_st))
                f.write("PRIMCOORD\n%i %i\n" % (n_atom, 1))
                for ix_atom in range(len(list_atom)):
                    f.write("%i %s\n" % (list_atom[ix_atom],
                        " ".join(["%.6f" % x for x in ar_posforce[ix_atom,:]])))

