#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

################################################################
direction = "z"
carrier = "e"
b2_e = [0.375, 0.310, 0.347] # from fitQuality
b2_h = [0.023, 0.027, 0.023]
with_ee_scattering = True
################################################################


plt.rcParams["figure.titlesize"] = 22
plt.rcParams["lines.linewidth"] = 2
plt.rcParams["xtick.labelsize"] = 22
plt.rcParams["ytick.labelsize"] = 22
plt.rcParams["font.size"] = 22
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = "Arial"
plt.rcParams["legend.fontsize"] = 18
plt.rcParams["figure.figsize"] = (8, 6)
plt.rcParams["figure.dpi"] = 200

# combine the name of the file to read with the current work direcotry
cwd = os.getcwd()
dir_f = os.path.join(cwd, "spin_lifetime.xlsx")

# read calc. data from the sheets
ws = pd.ExcelFile(dir_f).parse("tau_s_pristine")
T = ws["T"]
elec_tau_sx = ws["elec_tau_sx"]
elec_tau_sy = ws["elec_tau_sy"]
elec_tau_sz = ws["elec_tau_sz"]
hole_tau_sx = ws["hole_tau_sx"]
hole_tau_sy = ws["hole_tau_sy"]
hole_tau_sz = ws["hole_tau_sz"]
elec_tau_carrier_rate_avg_eph = ws["elec_tau_carrier_rate_avg_eph"]
elec_tau_carrier_rate_avg_eph_ee = ws["elec_tau_carrier_rate_avg_eph_ee"]
hole_tau_carrier_rate_avg_eph = ws["hole_tau_carrier_rate_avg_eph"]
hole_tau_carrier_rate_avg_eph_ee = ws["hole_tau_carrier_rate_avg_eph_ee"]

b2_e = np.asarray(b2_e)
b2_h = np.asarray(b2_h)

elec_tau_s_eph = np.zeros((3, len(elec_tau_carrier_rate_avg_eph)))
elec_tau_s_eph_ee = np.zeros((3, len(elec_tau_carrier_rate_avg_eph_ee)))
hole_tau_s_eph = np.zeros((3, len(hole_tau_carrier_rate_avg_eph)))
hole_tau_s_eph_ee = np.zeros((3, len(hole_tau_carrier_rate_avg_eph_ee)))

for i in range(3):
    elec_tau_s_eph_ee[i] = elec_tau_carrier_rate_avg_eph_ee/4/b2_e[i]
    hole_tau_s_eph_ee[i] = hole_tau_carrier_rate_avg_eph_ee/4/b2_h[i]


fig, ax = plt.subplots(nrows=1, ncols=1, constrained_layout=True)

#### spin lifetime calc. by spin mixing ####
if carrier == "h":
    if direction == "x":
        ax.plot(
            T, hole_tau_s_eph_ee[0], 
            marker="D", markersize=10, markerfacecolor="white", color="tab:red", label="hole (EY)"
        )
    elif direction == "y":
        ax.plot(
            T, hole_tau_s_eph_ee[1], 
            marker="D", markersize=10, markerfacecolor="white", color="tab:red", label="hole (EY)"
        )
    else:
        ax.plot(
            T, hole_tau_s_eph_ee[2], 
            marker="D", markersize=10, markerfacecolor="white", color="tab:red", label="hole (EY)"
        )
else:
    if direction == "x":
        ax.plot(
            T, elec_tau_s_eph_ee[0], 
            marker="D", markersize=10, markerfacecolor="white", color="tab:blue", label="elec. (EY)"
        )
    elif direction == "y":
        ax.plot(
            T, elec_tau_s_eph_ee[1], 
            marker="D", markersize=10, markerfacecolor="white", color="tab:blue", label="elec. (EY)"
        )
    else:
        ax.plot(
            T, elec_tau_s_eph_ee[2], 
            marker="D", markersize=10, markerfacecolor="white", color="tab:blue", label="elec. (EY)"
        )


#### spin lifetime calc. by DMD ####
if carrier == "h":
    if direction == "x":
        ax.plot(
            T, hole_tau_sx, 
            marker="o", markersize=10, color="tab:red", label="hole (calc. by DMD)"
        )
    elif direction == "y":
        ax.plot(
            T, hole_tau_sz, # lattice constant z of MAPbBr3 is the second longest so set it as y
            marker="o", markersize=10, color="tab:red", label="hole (calc. by DMD)"
        )
    else:
        ax.plot(
            T, hole_tau_sy, # lattice constant y of MAPbBr3 is the longest so set it as z 
            marker="o", markersize=10, color="tab:red", label="hole (calc. by DMD)"
        )
else:
    if direction == "x":
        ax.plot(
            T, elec_tau_sx, 
            marker="o", markersize=10, color="tab:blue", label="elec. (calc. by DMD)"
        )
    elif direction == "y":
        ax.plot(
            T, elec_tau_sz, # lattice constant z of MAPbBr3 is the second longest so set it as y 
            marker="o", markersize=10, color="tab:blue", label="elec. (calc. by DMD)"
        )
    else:
        ax.plot(
            T, elec_tau_sy, # lattice constant y of MAPbBr3 is the longest so set it as z  
            marker="o", markersize=10, color="tab:blue", label="elec. (calc. by DMD)"
        )


ax.set_yscale("log")
ax.legend()
ax.set_xlabel("Temperature (K)")
ax.set_ylabel("$\mathrm{\\tau_s~(ps)}$")
plt.savefig("tau_s_vs_temp.png")
