#!/usr/bin/env python
from argparse import ArgumentParser
import os
import sys
import re

parser = ArgumentParser(description="Copy finished NEB calculation as a new neb input with all intermidate geomtetry copied")
parser.add_argument('-d', dest="dirneb", required=True, help='NEB directory with *_1/PW.out as output')
parser.add_argument('-r', dest="ref", default=None, help='NEB reference input file, geometry will be replaced; if not specified, use neb.in in -d folder')
parser.add_argument('-i', dest="input", required=True, help='NEB input file without geometry')
args = parser.parse_args()

if (args.ref is None):
    args.ref = os.path.join(args.dirneb, "neb.in")

dic_geometry = {}
for dirroot, dirnames, filenames in os.walk(args.dirneb):
    for dirname in dirnames:
        filename = os.path.join(dirroot, dirname, "PW.out")
        if (os.path.exists(filename)):
            n = int(re.match(".+_(\\d+)", dirname).group(1))
            with open(filename, 'r') as f:
                lines = f.readlines()

            for i in range(len(lines)-1, -1, -1):
                if ("ATOMIC_POSITIONS" in lines[i]):
                    break
            list_pos = []
            while (lines[i].strip() != ""):
                list_pos.append(lines[i])
                i += 1
            dic_geometry[n] = list_pos

list_geometry = list(dic_geometry.items())
list_geometry.sort(key=lambda x:x[0])
if (len(list_geometry) != list_geometry[-1][0]):
    raise ValueError("Missed geometry, only found : %s" % ([x[0] for x in list_geometry]))
list_geometry = [x[1] for x in list_geometry]
                
with open(args.ref, 'r') as f:
    lines = f.readlines()

for i, line in enumerate(lines):
    if ("BEGIN_POSITIONS" in line):
        ix_start = i
    elif ("END_POSITIONS" in line):
        ix_end = i
        break

if (os.path.exists(args.input)):
    raise ValueError("%s exists, will not overwrite" % args.input)

with open(args.input, 'w') as f:
    for i in range(0, ix_start+1):
        f.write(lines[i])
    for ix, list_pos in enumerate(list_geometry):
        if (ix == 0):
            f.write("FIRST_IMAGE\n")
        elif (ix == len(list_geometry) - 1):
            f.write("LAST_IMAGE\n")
        else:
            f.write("INTERMEDIATE_IMAGE\n")
        print(list_pos)
        for line in list_pos:
            f.write(line)
    for i in range(ix_end, len(lines)):
        f.write(lines[i])



