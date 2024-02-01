#!/usr/bin/env python3

import sys

try:
    filename = sys.argv[1]
except IndexError:
    print('Please provide path to orca output')
    sys.exit(0)

# cm-1 to GHz
cm_to_GHz = 29.9702547

with open(filename, 'r') as f:
    for line in f:
        if "D   =" in line:
            D = float(line.split()[2]) * cm_to_GHz
        elif "E/D =" in line:
            E = float(line.split()[2]) * D * cm_to_GHz
        elif "SPIN-SPIN " in line:
            Dss = float(line.split()[2]) * cm_to_GHz
            Ess = float(line.split()[3]) * cm_to_GHz
        elif line.startswith("SPIN-ORBIT "):
            Dso = float(line.split()[2]) * cm_to_GHz
            Eso = float(line.split()[3]) * cm_to_GHz


print("             D(GHz)           E(GHz)")
print("SPIN-SPIN  = {:10.6f}       {:10.6f}".format(Dss,Ess))
print("SPIN-ORBIT = {:10.6f}       {:10.6f}".format(Dso,Eso))
print("----------------------------------------")
print("TOTAL      = {:10.6f}       {:10.6f}".format(D,E))


