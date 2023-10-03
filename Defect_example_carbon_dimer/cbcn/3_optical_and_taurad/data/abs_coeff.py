#!/usr/bin/env python
import numpy as np
import os
import matplotlib.pyplot as plt

plt.rcParams["figure.titlesize"] = 18
plt.rcParams["lines.linewidth"] = 2
plt.rcParams["xtick.labelsize"] = 18
plt.rcParams["ytick.labelsize"] = 18
plt.rcParams["font.size"] = 18
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = "Arial"
plt.rcParams["legend.fontsize"] = 18
plt.rcParams["figure.figsize"] = (8, 6)
plt.rcParams["figure.dpi"] = 200


data = np.loadtxt("o-x.eps_q1_diago_bse")

omega = data[:, 0]
eps_1 = data[:, 2]
eps_2 = data[:, 1]

def calc_alpha(omega, eps_1, eps_2):
    y = omega * eps_2 / (np.sqrt((eps_1 + np.sqrt(eps_1**2 + eps_2**2))/2.0))
    return y

unit_conv = 2.0 * np.pi * 8065.6e-3

y = calc_alpha(omega, eps_1, eps_2) * unit_conv

fig, ax = plt.subplots(nrows=1, ncols=1, constrained_layout=True)

ax.plot(omega, eps_1, label="Re")
ax.plot(omega, eps_2, label="Im")
ax.legend()
ax.set_xlim(4, 8)
ax.set_xlabel("Energy (eV)")
#ax.set_ylabel("$\mathrm{\\alpha~(cm^{-1})}$")
ax.set_ylabel("$\mathrm{\\varepsilon~(a.u.)}$")
plt.savefig("abs_coeff.png")

