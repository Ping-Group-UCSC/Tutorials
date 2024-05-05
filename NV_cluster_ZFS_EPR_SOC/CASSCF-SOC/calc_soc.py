#!/usr/bin/env python3
import numpy as np
import os
import re
import sys

cwd = os.getcwd()
print("Enter file name to read: ")
f = open(os.path.join(cwd, sys.argv[1]), "r")
print("Reading {}".format(sys.argv[1]))
lines = f.readlines()
num_lines = len(lines)

def get_soc(line):
    cm2GHz = 29979245.8e-6
    re_soc = float(line.split()[8])
    im_soc = float(line.split()[9])
    soc = np.sqrt(re_soc**2 + im_soc**2)
    soc *= cm2GHz
    return soc



for i, line in enumerate(lines):
# Block 0 for triplet, Block 1 for singlet
# Root starts from 0 (ground state)
#              Bra                       Ket
#   <Block Root  S    Ms  | HSOC |  Block Root  S    Ms>
    if "0    2  1.0  0.0              0     1  1.0  0.0" in line: # within 3E(ms=0), lambda z
        soc = get_soc(line)
        print("3E(ms=0) (z): ", np.round(soc, 3), "GHz")
    if "0    2  1.0 -1.0              0     1  1.0  1.0" in line: # within 3E(ms=-/+1), lambda z
        soc = get_soc(line)
        print("3E(ms=-/+1) (z): ", np.round(soc, 3), "GHz")
    if "0    2  1.0  1.0              0     1  1.0 -1.0" in line: # within 3E(ms=+/-1), lambda z
        soc = get_soc(line)
        print("3E(ms=+/-1) (z): ", np.round(soc, 3), "GHz")
    if "0    2  1.0  1.0              0     1  1.0  1.0" in line: # within 3E(ms=1), lambda z
        soc = get_soc(line)
        print("3E(ms=1) (z): ", np.round(soc, 3), "GHz")
    if "0    2  1.0 -1.0              0     1  1.0 -1.0" in line: # within 3E(ms=-1), lambda z
        soc = get_soc(line)
        print("3E(ms=-1) (z): ", np.round(soc, 3), "GHz")
    # ---------------------------------------------------------------- #
#              Bra                       Ket
#   <Block Root  S    Ms  | HSOC |  Block Root  S    Ms>
    if "1    3  0.0  0.0              0     1  1.0  0.0" in line: # 1E'(root3)->3E(root1), lambda z
        soc = get_soc(line)
        print("1E'(root3)->3E(root1) (z): ", np.round(soc, 3), "GHz")
    if "1    3  0.0  0.0              0     1  1.0  1.0" in line: # 1E'(root3)->3E(root1), lambda perp
        soc = get_soc(line)
        print("1E'(root3)->3E(root1) (⟂): ", np.round(soc, 3), "GHz")
    if "1    4  0.0  0.0              0     1  1.0  0.0" in line: # 1E'(root4)->3E(root1), lambda z
        soc = get_soc(line)                                                              
        print("1E'(root4)->3E(root1) (z): ", np.round(soc, 3), "GHz")
    if "1    4  0.0  0.0              0     1  1.0  1.0" in line: # 1E'(root4)->3E(root1), lambda perp
        soc = get_soc(line)
        print("1E'(root4)->3E(root1) (⟂): ", np.round(soc, 3), "GHz")
    # ---------------------------------------------------------------- #
#              Bra                       Ket
#   <Block Root  S    Ms  | HSOC |  Block Root  S    Ms>
    if "1    3  0.0  0.0              0     2  1.0  0.0" in line: # 1E'(root3)->3E(root2), lambda z
        soc = get_soc(line)
        print("1E'(root3)->3E(root2) (z): ", np.round(soc, 3), "GHz")
    if "1    3  0.0  0.0              0     2  1.0  1.0" in line: # 1E'(root3)->3E(root2), lambda perp
        soc = get_soc(line)
        print("1E'(root3)->3E(root2) (⟂): ", np.round(soc, 3), "GHz")
    if "1    4  0.0  0.0              0     2  1.0  0.0" in line: # 1E'(root4)->3E(root2), lambda z
        soc = get_soc(line)                                                             
        print("1E'(root4)->3E(root2) (z): ", np.round(soc, 3), "GHz")
    if "1    4  0.0  0.0              0     2  1.0  1.0" in line: # 1E'(root4)->3E(root2), lambda perp
        soc = get_soc(line)
        print("1E'(root4)->3E(root2) (⟂): ", np.round(soc, 3), "GHz")
    # ---------------------------------------------------------------- #
#              Bra                       Ket
#   <Block Root  S    Ms  | HSOC |  Block Root  S    Ms>
    if "1    2  0.0  0.0              0     1  1.0  0.0" in line: # 3E(root1)->1A1, lambda z
        soc = get_soc(line)
        print("3E(root1)->1A1 (z): ", np.round(soc, 3), "GHz")
    if "1    2  0.0  0.0              0     1  1.0  1.0" in line: # 3E(root1)->1A1, lambda perp
        soc = get_soc(line)
        print("3E(root1)->1A1 (⟂): ", np.round(soc, 3), "GHz")
    if "1    2  0.0  0.0              0     2  1.0  0.0" in line: # 3E(root2)->1A1, lambda z
        soc = get_soc(line)
        print("3E(root2)->1A1 (z): ", np.round(soc, 3), "GHz")
    if "1    2  0.0  0.0              0     2  1.0  1.0" in line: # 3E(root2)->1A1, lambda perp
        soc = get_soc(line)
        print("3E(root2)->1A1 (⟂): ", np.round(soc, 3), "GHz")
    # ---------------------------------------------------------------- #
#              Bra                       Ket
#   <Block Root  S    Ms  | HSOC |  Block Root  S    Ms>
    if "1    0  0.0  0.0              0     1  1.0  0.0" in line: # 3E(root1)->1E(root0), lambda z
        soc = get_soc(line)                                                             
        print("3E(root1)->1E(root0) (z): ", np.round(soc, 3), "GHz")                           
    if "1    0  0.0  0.0              0     1  1.0  1.0" in line: # 3E(root1)->1E(root0), lambda perp
        soc = get_soc(line)                                                                    
        print("3E(root1)->1E(root0) (⟂): ", np.round(soc, 3), "GHz")                           
    if "1    1  0.0  0.0              0     1  1.0  0.0" in line: # 3E(root1)->1E(root1), lambda z
        soc = get_soc(line)                                                                    
        print("3E(root1)->1E(root1) (z): ", np.round(soc, 3), "GHz")                           
    if "1    1  0.0  0.0              0     1  1.0  1.0" in line: # 3E(root1)->1E(root1), lambda perp
        soc = get_soc(line)        
        print("3E(root1)->1E(root1) (⟂): ", np.round(soc, 3), "GHz")
    # ---------------------------------------------------------------- #
#              Bra                       Ket
#   <Block Root  S    Ms  | HSOC |  Block Root  S    Ms>
    if "1    0  0.0  0.0              0     2  1.0  0.0" in line: # 3E(root2)->1E(root0), lambda z
        soc = get_soc(line)
        print("3E(root2)->1E(root0) (z): ", np.round(soc, 3), "GHz")
    if "1    0  0.0  0.0              0     2  1.0  1.0" in line: # 3E(root2)->1E(root0), lambda perp
        soc = get_soc(line)
        print("3E(root2)->1E(root0) (⟂): ", np.round(soc, 3), "GHz")
    if "1    1  0.0  0.0              0     2  1.0  0.0" in line: # 3E(root2)->1E(root1), lambda z
        soc = get_soc(line)
        print("3E(root2)->1E(root1) (z): ", np.round(soc, 3), "GHz")
    if "1    1  0.0  0.0              0     2  1.0  1.0" in line: # 3E(root2)->1E(root1), lambda perp
        soc = get_soc(line)
        print("3E(root2)->1E(root1) (⟂): ", np.round(soc, 3), "GHz")
    # ---------------------------------------------------------------- #
#              Bra                       Ket
#   <Block Root  S    Ms  | HSOC |  Block Root  S    Ms>
    if "1    2  0.0  0.0              0     0  1.0  0.0" in line: # 1A1->3A2, lambda z
        soc = get_soc(line)
        print("1A1->3A2 (z): ", np.round(soc, 3), "GHz")
    if "1    2  0.0  0.0              0     0  1.0  1.0" in line: # 1A1->3A2, lambda perp
        soc = get_soc(line)
        print("1A1->3A2 (⟂): ", np.round(soc, 3), "GHz")
    # ---------------------------------------------------------------- #
#              Bra                       Ket
#   <Block Root  S    Ms  | HSOC |  Block Root  S    Ms>
    if "1    0  0.0  0.0              0     0  1.0  0.0" in line: # 1E(root0)->3A2, lambda z
        soc = get_soc(line)
        print("1E(root0)->3A2 (z): ", np.round(soc, 3), "GHz")
    if "1    0  0.0  0.0              0     0  1.0  1.0" in line: # 1E(root0)->3A2, lambda perp
        soc = get_soc(line)
        print("1E(root0)->3A2 (⟂): ", np.round(soc, 3), "GHz")
    if "1    1  0.0  0.0              0     0  1.0  0.0" in line: # 1E(root1)->3A2, lambda z
        soc = get_soc(line)
        print("1E(root1)->3A2 (z): ", np.round(soc, 3), "GHz")
    if "1    1  0.0  0.0              0     0  1.0  1.0" in line: # 1E(root1)->3A2, lambda perp
        soc = get_soc(line)
        print("1E(root1)->3A2 (⟂): ", np.round(soc, 3), "GHz")

