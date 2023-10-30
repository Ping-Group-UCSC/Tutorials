#!/usr/bin/env python
#Read from different packages
import os
import numpy as np
from numpy.linalg import norm, inv

import copy
import itertools
from collections import Counter

Ang2m = 1e-10
Bohr2m = 0.52917720859E-10

Ang2Bohr = Ang2m / Bohr2m
Bohr2Ang = Bohr2m / Ang2m


def read_cell_and_pos_auto(arg):
    '''
    Detect file type and read lattice information
    All returned positions are in unit of angstrom
    '''
    if (os.path.exists(arg + ".in")):
        return read_cell_and_pos_qe(arg), "qe"
    elif ("POSCAR" in arg or "CONTCAR" in arg or arg.endswith(".vasp")):
        return read_cell_and_pos_poscar(arg), "vasp"
    else:
        return None, None

def write_cell_and_pos_auto(package, *args, **kwargs):
    '''
    Write into given format
    '''
    if (package == "vasp"):
        write_poscar(*args, **kwargs)
    elif (package == "qe"):
        write_cell_and_pos_qe(*args, **kwargs)
    else:
        raise ValueError("Unknown package %i" % package)
    return


def read_cell_and_pos_qe(prefix):
    '''
    Read the cell and positions from .in file and .out file (if exists)
    Note all numbers after species (like V1, Fe1) will be removed :  to (V,Fe)
    The speceis + number is save in speciesfull

    vecR must be read from CELL_PARAMETERS from input; if that does not exists, read output instead

    :return: vecR (lattice vector in columns),  list_atom(dictionary with species and pos)
    '''
#Read cell parameters
    with open(prefix + ".in", 'r') as f:
        print("Read from QE input")
        lines = f.readlines()

    vecR = None
    list_pos = []
    for i, line in enumerate(lines):
        if ("nat" in line):
            ar = line.split(",")
            for st2 in ar:
                if ("nat" in st2):
                    ar2 = st2.split("=")
                    nat = int(ar2[-1])
        elif ("CELL_PARAMETER" in line):
            scale = 1
            if ("ang" in line.lower()):
                pass
            elif ("bohr" in line.lower()):
               scale = Bohr2Ang  
            else:
                raise ValueError("Only angstrom unit is supported")
            vecR = np.asarray([[float(x) for x in line.split()] for line in lines[i+1:i+4]]).T * scale
            vecRi = inv(vecR)
            break

    for i, line in enumerate(lines):
        if ("ATOMIC_POSITIONS" in line):
            unit_coord = line.split()[-1].replace("{","").replace("}","").replace("(","").replace(")","").lower()
            if ("crystal" not in line and vecR is None):
                print("Only crystal coordinate is supported in QE input without CELL_PARAMETERS")
                break
            for line2 in lines[i+1:i+nat+1]:
                ar = line2.strip().split()
                pos = np.asarray([float(x) for x in ar[1:]])
                if (unit_coord == "crystal"):
                    pass
                elif (unit_coord == "angstrom"):
                    pos = np.dot(vecRi, pos)
                elif (unit_coord == "bohr"):
                    pos = np.dot(vecRi, pos) / Bohr2Ang
                else:
                    raise ValueError("Unsupported unit %s" % unit_coord)
                list_pos.append({"species": filter(lambda x:x.isalpha(), ar[0]), "pos" : pos, "speciesfull" : ar[0]})

    if (os.path.exists(prefix + ".out")):
        print("Read from QE output")
        with open(prefix + ".out", 'r') as f:
            lines = f.readlines()

        if (vecR is None):#Not found in input ; read initial structure in output
            for i, line in enumerate(lines):
                if ("celldm(1)" in line):
                    i += 4
                    break
            celldm1 = float(line[15:26])
            vecR = np.asarray([[float(x)*celldm1*Bohr2Ang for x in line[23:56].split()] for line in lines[i:i+3]]).T

        if (vecR is None):
            raise ValueError("Cannot find cell parameters information or output")

        i_start = None
        for i, line in enumerate(lines):
            if (line.startswith("Begin final coordinates")):
                i_start = i + 1
                break

        if (i_start is None):
#Do not read .out if .in is read
            if (len(list_pos) > 0):
                print("Not relax calculation, skip .out file")
            else: #Try read .out initial coordinate 
                for i, line in enumerate(lines):
                    if ("Crystallographic axes" in line):
                        i_start =  i + 3
                        break
                i = i_start
                while (len(lines[i]) > 1):
#                   print(lines[i])
                    ar = lines[i].split()
                    pos = np.array([float(x) for x in ar[-4:-1]])
#                   print(pos)
                    species = ar[1]
                    list_pos.append({"species": filter(lambda x:x.isalpha(), species), "pos" : pos, "speciesfull" : species})
                    i += 1
        else:
        #Check if it is vc-relax
#Read cell parameters
            if ("CELL" in lines[i_start+3]):#vc-relax
                vecR = np.asarray([[float(x) for x in line.split()] for line in lines[i_start+4:i_start+7]])
                if ("Ang" not in lines[i_start+3]):#Convert unit
                    vecR *= 0.529177249
                i = i_start + 8
            else: #relax
#Read from start
                for i, line in enumerate(lines):
                    if ("celldm" in line):
                        alat = float(line.split()[1]) * 0.529177249
                        break
                vecR = np.asarray([[float(x) for x in line.split()[3:6]] for line in lines[i+4:i+7]]) * alat
                #Recover i
                i = i_start + 1
            vecR = vecR.T
            vecRi = inv(vecR)

#Read atoms
            unit_coord = lines[i].split()[-1].replace("{","").replace("}","").replace("(","").replace(")","").lower()

            list_pos = []
            while (True):
                i += 1
                line = lines[i]
                if (line.startswith("End")):
                    break
                ar = line.split()
                pos = np.asarray([float(x) for x in ar[1:]])
                if (unit_coord == "crystal"):
                    pass
                elif (unit_coord == "angstrom"):
                    pos = np.dot(vecRi, pos)
                else:
                    raise ValueError("Unsupported unit %s" % unit_coord)
                list_pos.append({"species": filter(lambda x:x.isalpha(), ar[0]), "pos" : pos, "speciesfull" : ar[0]})
#               print(list_pos[-1])


    return vecR, list_pos

def read_cell_and_pos_poscar(filename):
    '''
    Read the cell and positions
    '''
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    vecR = np.asarray([[float(x) for x in line.split()] for line in lines[2:5]]).T
    vecRi = np.linalg.inv(vecR)
    atoms = lines[5].split()
    natoms = [int(x) for x in lines[6].split()]
    unit = lines[7][0]
    i = 8
    if (unit == "S"):#Selective Dynamics; read next line
        unit = lines[8][0]
        i += 1

    list_pos = []

    for atom, n in zip(atoms, natoms):
        for j in range(n):
            pos = np.asarray([float(x) for x in lines[i].split()[:3]])
            if (unit == "C"):#Convert unit to crystal coordinate
                pos = np.dot(vecRi, pos)
            list_pos.append({"species": atom, "pos" : pos})
            i += 1

    return vecR, list_pos

def write_cell_and_pos_qe(filename, vecR, list_pos):
    '''
    Write into QE atomic position crystal format
    '''
    with open(filename, 'w') as f:
        f.write("! nat = %i\n" % len(list_pos))
        f.write("CELL_PARAMETERS (angstrom)\n")
        for i in range(3):
            f.write("%.10f %.10f %.10f\n" % tuple(vecR[:,i]))
        f.write("ATOMIC_POSITIONS (crystal)\n")
        for atom in list_pos:
            f.write("%-5s %14.9f %14.9f %14.9f\n" % (atom["speciesfull"],
                atom["pos"][0], atom["pos"][1], atom["pos"][2]))
    return


def write_poscar(filename, vecR, list_pos):
    '''
    Write to POSCAR file
    '''
    with open(filename, 'w') as f:
        f.write("auto structure\n1.0\n")
        for i in range(3):
            f.write("%.10f %.10f %.10f\n" % tuple(vecR[:,i]))

        list_species = [x["species"] for x in list_pos]
        dic_species = Counter(list_species)
#Sort dic_species as appearance order in list_species
        list_count = []
        for x in list_species:
            if (x in dic_species):
                list_count.append((x, dic_species[x]))
                del dic_species[x]

#       print(list_count)
        for species, val in list_count:
            f.write("%4s " % species)
        f.write("\n")
        for species, val in list_count:
            f.write("%4i " % val)
        f.write("\n")
        f.write("Direct\n")
        for species, val in list_count:
            for atom in list_pos:
                if (atom["species"] == species):
                    f.write("%20.16f %20.16f %20.16f\n" % tuple(atom["pos"]))

def read_pos_only_qe(filename):
    '''
    Read atomic positions in QE
    Only take care of atomic position, no unit conversion did
    '''
#Read cell parameters
    with open(filename, 'r') as f:
        lines = f.readlines()

    lines_before = []
    lines_after = []

    vecR = None
    list_atom = []
    for i, line in enumerate(lines):
        if ("nat" in line):
            ar = line.split(",")
            for st2 in ar:
                if ("nat" in st2):
                    ar2 = st2.split("=")
                    nat = int(ar2[-1])
        elif ("CELL_PARAMETER" in line):
            scale = 1
            if ("ang" in line.lower()):
                pass
            elif ("bohr" in line.lower()):
               scale = Bohr2Ang  
            else:
                raise ValueError("Only angstrom unit is supported")
            vecR = np.asarray([[float(x) for x in line.split()] for line in lines[i+1:i+4]]).T * scale
            vecRi = inv(vecR)
            break

    for i, line in enumerate(lines):
        if ("ATOMIC_POSITIONS" in line):
            lines_before = lines[:i]
            unit_coord = line.split()[-1].replace("{","").replace("}","").replace("(","").replace(")","").lower()
            for line2 in lines[i+1:i+nat+1]:
                ar = line2.strip().split()
                pos = np.asarray([float(x) for x in ar[1:]])
                if (unit_coord == "crystal"):
                    pass
                elif (unit_coord == "angstrom"):
                    pass
                elif (unit_coord == "bohr"):
                    pass
                else:
                    raise ValueError("Unsupported unit %s" % unit_coord)
                list_atom.append({"species": filter(lambda x:x.isalpha(), ar[0]), "pos" : pos, "speciesfull" : ar[0]})
            lines_after = lines[i+nat+1:]

    return vecR, list_atom, unit_coord, lines_before, lines_after

def write_pos_keep_qe(filename, vecR, list_atom, unit, lines_before, lines_after):
    '''
    Write QE with all other parts kept, only atomic positions user defined
    To be used with read_pos_only_qe
    '''
    with open(filename, 'w') as f:
        for line in lines_before:
            f.write(line)
        f.write("ATOMIC_POSITIONS (%s)\n" % unit)
        for atom in list_atom:
            f.write("%-5s %14.9f %14.9f %14.9f\n" % (atom["speciesfull"],
                atom["pos"][0], atom["pos"][1], atom["pos"][2]))
        for line in lines_after:
            f.write(line)
    return
