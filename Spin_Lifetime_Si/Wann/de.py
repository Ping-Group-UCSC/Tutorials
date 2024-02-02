#!/usr/bin/env python
import numpy as np

c1 = np.loadtxt("wannier.eigenvals", dtype=np.float64, usecols=(8))
c2 = np.loadtxt("wannier.eigenvals", dtype=np.float64, usecols=(9))
dc = c2 - c1
print("from k50 to k70 around CBM\n",dc[50:70])
