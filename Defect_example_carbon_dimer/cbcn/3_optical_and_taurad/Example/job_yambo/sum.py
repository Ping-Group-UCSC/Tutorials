#!/usr/bin/env python3
import numpy as np
import os


dirc = "z"
files = ["o-"+dirc+".exc_weights_at_1", "o-"+dirc+".exc_weights_at_2", 
    "o-"+dirc+".exc_weights_at_3", "o-"+dirc+".exc_weights_at_4", 
    "o-"+dirc+".exc_weights_at_5", "o-"+dirc+".exc_weights_at_6", "o-"+dirc+".exc_weights_at_7"]
for i in range(len(files)):
    dir = os.path.join("./", files[i])
    data = np.genfromtxt(dir, dtype=float)
    print(type(data))
    try:
        print(i)
        print("\n", data[:, 0])
        print(data[:, 1])
        print(data[:, 4])
        weight = sum(data[:, 5])
        print("ENERGY: ", data[:, -1])
        print("weight: ", weight)
    except: continue



