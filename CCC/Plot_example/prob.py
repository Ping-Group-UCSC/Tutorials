#!/usr/bin/env python3

import math
kT = 0.0257
def bolt(e):
    return math.exp(-e/kT)

eList = [0, 0.099428537, 0.438443348]
Z = sum([bolt(e) for e in eList])

for e in eList:
    print("State with e = {:.4f} has probability = {}".format(e, bolt(e)/Z))
