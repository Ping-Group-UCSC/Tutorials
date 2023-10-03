#!/usr/bin/env python
#This script calculate HR factor for all phonons according to Alkauskas New J. Phys. 16 (2014) 073026
#Support  pw.x + ph.x output
#Or Phonopy output
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
from math import sqrt
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

def is_option_set(term):
    '''
    Check if an options is set
    This is just a rough test with string
    '''
    return not term.startswith("The ") and not term.strip() == ""

def read_dynmat_mass(filename):
    '''
    Read Mass from dynmat file
    Convert to amu
    '''
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    ix_start = 0+3
    for i, line in enumerate(lines):
        if ("Basis" in line):
            ix_start = i+4

    i = ix_start
    dic_mass = {}
    while (lines[i][4] == " "):
        line = lines[i]
#Note the unit in ph.x is Ry-based Amu
        dic_mass[line[15:19].strip()] = float(line[21:]) / AMU2me * 2
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
            "eigenvector": np.array(eigenvec),
            })

    return list_data


def read_phonopy_struct(filename):
    '''
    Read structure from phonopy.yaml

    vecR in unit bohr
    pos in unit crystal

    Note format is converted to "pos" and "species" 
    '''
    data = read_yaml(filename)

    vecR = np.asarray(data["primitive_cell"]["lattice"]).T
    list_pos = data["primitive_cell"]["points"]
    for x in list_pos:
        x["pos"] = np.asarray(x["coordinates"])
        x["species"] = x["symbol"]

    return vecR, list_pos

def read_phonon_gamma_phonopy(filename, list_mass):
    '''
    A (much) fast version to read phonon information than YAML
    '''
    with open(filename, 'r') as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        if ("q-position" in line):
            if ("[    0.0000000,    0.0000000,    0.0000000 ]" in line):
                i0 = i
                break
            else:
                raise ValueError("First q-point is not Gamma point")

    list_data = []

    i = i0 + 4
    while ("frequency" in lines[i]):
        freq = float(lines[i].split()[-1])
        l1 = []
        l2 = []
        i += 2
        while ("atom" in lines[i]):
#Must be real; read only real parts
            l1.append([float(line.split()[2][:-1]) for line in lines[i+1:i+4]])
            l2.append([x / list_mass[len(l1)-1] for x in l1[-1]])
            i += 4
        ar1 = np.asarray(l1)
        ar2 = np.asarray(l2)
        list_data.append({"frequency":freq, 
            "eigenvector": ar1,
            "displacement" : ar2
            })
        i += 1

    return list_data


file_input = "input.yaml"
file_store = "output.yaml"
file_store2 = "detail.yaml"

if (not os.path.exists(file_input)):
    dinput = yaml.load('''
title: HR 
pwx_prefix_init_relax : The folder name and prefix of relaxation or scf with initial state geometry
pwx_prefix_final_relax : The folder name and prefix of relaxation or scf with final state geometry
phonopy_folder : The folder name of phonopy calculation
phx_filename_phonon_dynmat:  The filename of dynamt file generated by ph.x
phx_filename_phonon_mold: The filename of .mold file generated by ph.x
''')
else:
    dinput = read_yaml(file_input)

if ("The " in dinput["pwx_prefix_init_relax"]):
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
    (vecRi, list_atom_init), package = read_cell_and_pos_auto(dinput["pwx_prefix_init_relax"])
    (vecRf, list_atom_final), package = read_cell_and_pos_auto(dinput["pwx_prefix_final_relax"])
#Convert to Bohr
    vecRi *= Ang2Bohr
    vecRf *= Ang2Bohr

#This is to read the shared information (mass and species)
    list_atom_common = list_atom_init

    print(vecRi, vecRf)

    if  (np.abs(vecRi -vecRf) > 1e-5).any():
        raise ValueError("Lattice parameter of initial and final structure not consistent")
    vecR = vecRi

#Try to read Phonopy calculation
    if (is_option_set(dinput["phonopy_folder"])):
        print("Read Phonopy calculation...")
        vecR_phonopy, list_atom_phonopy = read_phonopy_struct(os.path.join(dinput["phonopy_folder"], "phonopy.yaml"))
        phonon = read_phonon_gamma_phonopy(os.path.join(dinput["phonopy_folder"], "band.yaml"), [x["mass"] for x in list_atom_phonopy])
#Append mass information to list_atom
        for x, x2 in zip(list_atom_common, list_atom_phonopy):
            x["mass"] = x2["mass"]
    else:
#Try to read ph.x calculation
        print("Read ph.x calculation...")
        dic_mass = read_dynmat_mass(dinput["phx_filename_phonon_dynmat"])
#Append mass information to list_atom
        for x in list_atom_common:
            x["mass"] = dic_mass[x["species"]]

        print("Loading state information...")
        phonon = read_phonon_gamma_mold(dinput["phx_filename_phonon_mold"], dic_mass)


#   phonon_init = read_yaml(os.path.join(dinput["folder_init_state"], "band.yaml"))
#   phonon_final = read_yaml(os.path.join(dinput["folder_final_state"], "band.yaml"))

#Compute eq(6) m^1/2 (R-R)
    list_diffR = []
    for i in range(len(list_atom_common)):
#       v1 = np.dot(vecR, list_atom_init[i]["pos"] - list_atom_final[i]["pos"])  * list_atom_init[i]["mass"]
        v1 = np.dot(vecR, list_atom_init[i]["pos"] - list_atom_final[i]["pos"])
        list_diffR.append(v1)
            
#Save 1D-effective
    list_dR2 = [np.dot(vecR, x["pos"] - x2["pos"]) for x, x2 in zip(list_atom_init, list_atom_final)]
    data2["1D-dR-Bohr"] = np.array(list_dR2).flatten().tolist()
    list_dQ2 = [np.dot(x, x) * list_atom_common[i]["mass"] for i, x in enumerate(list_dR2)]
    data["1D-dQ2-Bohr2"] = float(sum(list_dQ2))

#phonons are assumed to be the same in initial and final states 
    
    for case, phonon in [
            ["ph", phonon]
            ]:
        list_Sk = []
        list_dQk = []
        list_dQk2 = []
#Check normalization of eigenvectors
#Two orthogonal and normalization condition:
#\sum_ai \q_(k,ai)q(k2,ai) = delta k1,k2
#\sum_k q_(ai1) q_(ai2) = delta ai1, ai2

        mat_q_ki = np.stack([x["eigenvector"].flatten() for x in phonon])
#       q = 0
#       for band in phonon:
#           for shift in band["eigenvector"]:
#               q += np.dot(shift, shift)
#       q = int(q)
#       deg = len(list_atom_common) * 3

#Check normalization between phonons
        n1 = np.matmul(mat_q_ki, mat_q_ki.T)
        error = np.abs(n1 - np.eye(n1.shape[0]))
        error_diag = np.diag(error)
        print("Normalization of phonons error :  %f" % np.max(error_diag))
        print("Orthongalnality between phonons error :  %f" % np.max(error))

#       n1 = np.matmul(mat_q_ki.T, mat_q_ki)
#       error = np.abs(n1 - np.eye(n1.shape[0]))
#       error_diag = np.diag(error)
#       print(error_diag)
#       print("Normalization of phonons error :  %f" % np.max(error_diag))
#       print("Orthongalnality between atomxyz error :  %f" % np.max(error))

        if (np.max(error_diag) > 1e-4):
            raise ValueError("Eigenvectors not normalzied")
        else:
            print("Eigenvectors are normalized")

        for ix_phonon, band in enumerate(phonon):
#Here we need the angular frequency ; what calculated is ordinary frequency
            freq = band["frequency"] * 2 * pi
#Note eigenvector is m^0.5*dR
            qk2 = 0
            for ix_atom, shift_mr in enumerate(band["eigenvector"]):
#Square before sum  (which is wrong)
#               qk2 += list_atom_common[ix_atom]["mass"] * ((list_diffR[ix_atom] * shift_mr)**2).sum()

#Square after sum
                qk2 += list_atom_common[ix_atom]["mass"]**0.5 * ((list_diffR[ix_atom] * shift_mr)).sum()
            
            qk2 = qk2**2

#               qk2 += list_atom_common[ix_atom]["mass"] * ((list_diffR[ix_atom] * np.array([1.0/np.sqrt(deg), 1.0/np.sqrt(deg), 1.0/np.sqrt(deg)]))**2).sum()
            Sk  = freq*qk2/2
#Convert unit: 
#qk^2: AMU * Bohr^2 
#Freq : THz (ordinary)
#           Sk = Sk * AMU2me * 1e12 / second2au
            Sk = Sk * AMU2kg * Bohr2m**2 * 1e12 / hbar_Js
            list_Sk.append(Sk)
            list_dQk.append(sqrt(qk2))
            list_dQk2.append(qk2)

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
        data2["%s-dQk-Bohr" % case] = [float(x) for x in list_dQk]
        data2["%s-dQk2-Bohr2" % case] = [float(x) for x in list_dQk2]
        data["%s-total-dQk2-Bohr2" % case] =  float(ar_dQk2.sum()) 
        data2["%s-FreqOrd-meV" % case] = [float(band["frequency"] * THzOrd2meV) for band in phonon]

        ar_ix = np.flip(np.argsort(ar_Sk))
        data2["%s-IndexLargeS" % case] =  [[int(ix), float(v)] for ix, v in zip(ar_ix + 1, ar_Sk[ar_ix])]

#Plot for final phonon
        if (case == "finalph"):
            plot_phonon_mode_multi("phonon-final", vecR*Bohr2Ang, list_atom_final, [band["displacement"] for band in phonon], unit_pos="crystal", unit_eigen="bohr", scale=20)

#Plot the phonon between two structures
    plot_phonon_mode("phonon-1dcoord", vecR*Bohr2Ang, list_atom_final, np.asarray([x["pos"] - x2["pos"] for x, x2 in zip(list_atom_init, list_atom_final)]), unit_pos="crystal", unit_eigen="crystal", scale=4)


save()


