#!/usr/bin/env python
from math import pi
#import pandas as pd
import numpy as np
from numpy.linalg import norm
#from StringIO import StringIO

Bohr2Ang = 0.52917720859
Ry2eV = 13.6056917
Ha2eV = Ry2eV * 2
SpeedOfLight = 137
au2second = 2.418884326505E-017

dim = 0

excE = np.loadtxt("./o-x.exc_E_sorted_abs")

vecR = np.array([
[  28.48680253,  0.000000,  0.000000 ],
[ -14.24340068,  24.67030124,  0.000000 ],
[  0.000000,  0.000000,  28.34588983 ],
]).T

LxLy = np.cross(vecR[:,0], vecR[:,1])
Lz = norm(vecR[:,2])
vol = np.dot(vecR[:,2], LxLy)
LxLy = norm(LxLy)

f = open("taurad.dat","w")
f.write("%8s %7s %8s %8s %8s %8s\n" % ("Index", "Volume", "Peak", "Amp", "mu^2", "Lifetime"))
for ix_exc in range(1, 9+1):
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
    B = gamma0 / au2second * LxLy * Bohr2Ang ** 2 * 1e-16
#   l0.append((case, vol, peak, amp, mu2, lifetime))

    f.write("%8s %7.0f %.3f %.2e %.2e %s ps %.2e cm^2/s\n" % (case, vol, peak*Ha2eV, amp, mu2, lifetime, B))
#    f.write("%8s %7.0f %.3f %.2e %.2e %s ps\n" % (case, vol, peak*Ha2eV, amp, mu2, lifetime))
