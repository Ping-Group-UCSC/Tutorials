#!/usr/bin/env python
from math import pi
#import pandas as pd
import numpy as np
from numpy.linalg import norm
from StringIO import StringIO

Bohr2Ang = 0.52917720859
Ry2eV = 13.6056917
Ha2eV = Ry2eV * 2
SpeedOfLight = 137
au2second = 2.418884326505E-017

dim = 0

excE = np.loadtxt("./o.exc_E_sorted_abs")

vecR = np.array([
[  56.85708,  0.000000,  0.000000 ],
[ -28.42854,  49.23968,  0.000000 ],
[  0.000000,  0.000000,  20.00000 ],
]).T

LxLy = np.cross(vecR[:,0], vecR[:,1])
Lz = norm(vecR[:,2])
vol = np.dot(vecR[:,2], LxLy)
LxLy = norm(LxLy)

f = open("taurad.dat","w")
f.write("%8s %7s %8s %8s %8s %8s\n" % ("Index", "Volume", "Peak", "Amp", "mu^2", "Lifetime"))
for ix_exc in range(1, 10+1):
    case = ix_exc
    peak = excE[ix_exc-1, 0] / Ha2eV
    amp = excE[ix_exc-1, 4]
    mu2 = vol * amp / 8 / pi

    if (dim == 0):
        gamma0 = 4.0 / 3.0 / SpeedOfLight ** 3 *  peak ** 3* mu2
    elif (dim == 1):
        gamma0 = 2 * pi / SpeedOfLight ** 2 / Lz * peak ** 2 * mu2
    elif (dim == 2):
        gamma0 = 4 * pi / SpeedOfLight ** 1 / LxLy * peak ** 1 * mu2

    lifetime = 1 / gamma0 * au2second * 1e12
#   l0.append((case, vol, peak, amp, mu2, lifetime))

    f.write("%8s %7.0f %.3f %.2e %.2e %s ps\n" % (case, vol, peak, amp, mu2, lifetime))
