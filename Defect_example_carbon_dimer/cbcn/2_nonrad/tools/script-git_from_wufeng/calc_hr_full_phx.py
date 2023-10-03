#!/usr/bin/env python
#This script calculate HR factor for all phonons according to Alkauskas New J. Phys. 16 (2014) 073026
#Support only one phonon and ph.x output
#Output format :

#
#case = initph or finalph ( calculated from initial or final structure)
#["%s-Sk" % case]  : all S_k for each phonon mode
#["%s-S" % case] : summation of S_k
#["%s-IndexLargeS" % case] : band index (1-based) sorted in descending order

import sys
import os
import yaml
import numpy as np
from constant import *
from numpy.linalg import norm
from libphonon import plot_phonon_mode_multi, plot_phonon_mode
from io_package import read_cell_and_pos_auto
np.seterr(all="log")

def read_yaml(filename):
    with open(filename, 'r') as f:
        d = yaml.load(f.read())
    return d

def saveinput():
    with open(file_input, 'w') as f:
        yaml.dump(dinput, f, default_flow_style=False, width=200)

def save():
    with open(file_store, 'w') as f:
        yaml.dump(data, f, default_flow_style=False, width=200)

    with open(file_store2, 'w') as f:
        yaml.dump(data2, f, default_flow_style=False, width=200)

def read_dynmat_mass(filename):
    '''
    Read Mass from dynmat file
    Convert to amu
    '''
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    for i, line in enumerate(lines):
        if ("Basis" in line):
            ix_start = i+4

    i = ix_start
    dic_mass = {}
    while (lines[i][4] == " "):
        line = lines[i]
        dic_mass[line[15:19].strip()] = float(line[21:]) / AMU2me
        i += 1

    return dic_mass


def read_phonon_gamma_mold(filename, dic_mass):
    '''
    Read eigenvalues and eigenvectors from dynmat.mold output
    Note: the atom list is abandoned
    '''
    with open(filename, 'r') as f:
        lines = f.readlines()

    state = "none"
#Freq here is cm-1
    list_freq = []
    list_atom = []
    for i, line in enumerate(lines):
        if (line.startswith("[")):
            state = line.strip()
        else:
#Convert unit 
            if state == "[FREQ]":
                list_freq.append(float(line) / THz2cmm1)
            elif state == "[FR-COORD]":
                list_atom.append(line)
            elif state == "[FR-NORM-COORD]":
                ix_start = i
                break

    n_atom = len(list_atom)

    list_eigenvec = []
    for ix_v in range(n_atom*3):
        eigenvec = []
        for ix_a in range(n_atom):
            ix = ix_start + ix_v * (n_atom+1) + ix_a + 1
            line = lines[ix]
            eigenvec.append(np.asarray(map(float, line.split())))
        list_eigenvec.append(eigenvec)

    list_data = []
    for eigenvec, freq in zip(list_eigenvec, list_freq):
        list_data.append({"frequency":freq, 
            "eigenvector": eigenvec,
            })

    return list_data

file_input = "input.yaml"
file_store = "output.yaml"
file_store2 = "detail.yaml"

if (not os.path.exists(file_input)):
    dinput = yaml.load('''
title: HR 
prefix_init_relax : The folder name and prefix of relaxation or scf with initial state geometry
prefix_final_relax : The folder name and prefix of relaxation or scf with final state geometry
''')
else:
    dinput = read_yaml(file_input)

if ("The " in dinput["prefix_init_relax"]):
    print("Please modify %s according to instruction inside" % file_input)
    saveinput()
    sys.exit(1)

if (not os.path.exists(file_store)):
    data = {}
    data2 = {}
else:
    data = read_yaml(file_store)
    data2 = read_yaml(file_store2)


#Compuate Huang-Rhys factor for each mode
if (not "HR_k" in data):
#Load structures
#Note all units in Bohr
    (vecRi, list_atom_init), package = read_cell_and_pos_auto(dinput["prefix_init_relax"])
    (vecRf, list_atom_final), package = read_cell_and_pos_auto(dinput["prefix_final_relax"])
#Convert to Bohr
    vecRi *= Ang2Bohr
    vecRf *= Ang2Bohr

#This is to read the shared information (mass and species)
    list_atom_common = list_atom_init

    print(vecRi, vecRf)

    if  (np.abs(vecRi -vecRf) > 1e-5).any():
        raise ValueError("Lattice parameter of initial and final structure not consistent")
    vecR = vecRi

    dic_mass = read_dynmat_mass(dinput["filename_phonon_dynmat"])
#Append mass information to list_atom
    for x in list_atom_common:
        x["mass"] = dic_mass[x["species"]]

    print("Loading state information...")
    phonon = read_phonon_gamma_mold(dinput["filename_phonon_mold"], dic_mass)


#   phonon_init = read_yaml(os.path.join(dinput["folder_init_state"], "band.yaml"))
#   phonon_final = read_yaml(os.path.join(dinput["folder_final_state"], "band.yaml"))

#Compute eq(6) m^1/2 (R-R)
    list_diffR = []
    for i in range(len(list_atom_common)):
#       v1 = np.dot(vecR, list_atom_init[i]["pos"] - list_atom_final[i]["pos"])  * list_atom_init[i]["mass"]
        v1 = np.dot(vecR, list_atom_init[i]["pos"] - list_atom_final[i]["pos"])
        list_diffR.append(v1)
            

#phonons are assumed to be the same in initial and final states 
    
    for case, phonon in [
            ["ph", phonon]
            ]:
        list_Sk = []
        list_dQk2 = []
#Check normalization of eigenvectors
        q = 0
        for band in phonon:
            for shift in band["eigenvector"]:
                q += np.dot(shift, shift)
        q = int(q)
        deg = len(list_atom_common) * 3
        print("\sum_{ik}q_{ik}q_{ik} = %i for %i degreee of freedom" % (q, deg))

        if (q != deg):
            raise ValueError("Eigenvectors not normalzied")
        else:
            print("Eigenvectors are normalized")

        for band in phonon:
#Here we need the angular frequency ; what calculated is ordinary frequency
            freq = band["frequency"] * 2 * pi
            qk2 = 0
#Note eigenvector is m^0.5*dR
            for ix_atom, shift_mr in enumerate(band["eigenvector"]):
                qk2 += list_atom_common[ix_atom]["mass"] * ((list_diffR[ix_atom] * shift_mr)**2).sum()
#               qk2 += list_atom_common[ix_atom]["mass"] * ((list_diffR[ix_atom] * np.array([1.0/np.sqrt(deg), 1.0/np.sqrt(deg), 1.0/np.sqrt(deg)]))**2).sum()
            Sk  = freq*qk2/2
            dQk2 = qk2
#Convert unit: 
#qk^2: AMU * Bohr^2 
#Freq : THz (ordinary)
#           Sk = Sk * AMU2me * 1e12 / second2au
            Sk = Sk * AMU2kg * Bohr2m**2 * 1e12 / hbar_Js
            list_Sk.append(Sk)
            list_dQk2.append(dQk2)

        ar_Sk = np.asarray(list_Sk)
        ar_dQk2 = np.asarray(list_dQk2)
        S = ar_Sk.sum()
        dQk2 = ar_dQk2.sum()
#       data2["%s-Sk" % case] = ["%.8e" % x for x in list_Sk]
#       data["%s-S" % case] = "%.8e" % S
#       data2["%s-dQk2" % case] = ["%.8e" % x for x in list_dQk2]
#       data["%s-dQ2" % case] = "%.8e" % (ar_dQk2.sum()) 
        data2["%s-Sk" % case] = [float(x) for x in list_Sk]
        data["%s-S" % case] =  float(S)
        data2["%s-dQk2-Bohr2" % case] = [float(x) for x in list_dQk2]
        data["%s-dQ2-Bohr2" % case] =  float(ar_dQk2.sum()) 
        data2["%s-FreqOrd-meV" % case] = [float(band["frequency"] * THzOrd2meV) for band in phonon]

        ar_ix = np.flip(np.argsort(ar_Sk))
        data2["%s-IndexLargeS" % case] =  [[int(ix), float(v)] for ix, v in zip(ar_ix + 1, ar_Sk[ar_ix])]

#Plot for final phonon
        if (case == "finalph"):
            plot_phonon_mode_multi("phonon-final", vecR*Bohr2Ang, list_atom_final, [band["displacement"] for band in phonon], unit_pos="crystal", unit_eigen="bohr", scale=20)

#Plot the phonon between two structures
    plot_phonon_mode("phonon-1dcoord", vecR*Bohr2Ang, list_atom_final, np.asarray([x["pos"] - x2["pos"] for x, x2 in zip(list_atom_init, list_atom_final)]), unit_pos="crystal", unit_eigen="crystal", scale=4)

    list_dQ2 = [np.dot(vecR, x["pos"] - x2["pos"]) for x, x2 in zip(list_atom_init, list_atom_final)]
    list_dQ2 = [np.dot(x, x) * list_atom_common[i]["mass"] for i, x in enumerate(list_dQ2)]
    data["dQ2-Bohr2"] = float(sum(list_dQ2))

save()


