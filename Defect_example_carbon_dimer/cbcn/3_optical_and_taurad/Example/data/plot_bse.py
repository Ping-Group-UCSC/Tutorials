#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import os
import re

plt.rcParams["lines.linewidth"] = 2
plt.rcParams["xtick.labelsize"] = 18
plt.rcParams["ytick.labelsize"] = 18
plt.rcParams["font.size"] = 18
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = "Arial"
plt.rcParams["legend.fontsize"] = 18
plt.rcParams["figure.figsize"] = (8, 6)
plt.rcParams["figure.dpi"] = 200

fig, ax = plt.subplots(nrows=1, ncols=1, constrained_layout=True)

cwd = os.getcwd()
file_to_plot = []
labels = []
colors = {"x":"tab:blue", "y":"tab:red", "z":"tab:green"}
z_length=30 # bohr
scale_factor = 1/z_length*4*np.pi

for f in os.listdir(cwd):
    if "o" in f and "eps_" in f:
        file_to_plot.append(os.path.join(cwd, f))
        labels.append(re.findall("x|y|z", f)[0])
print(file_to_plot)
for key in colors:
    for i, f in enumerate(file_to_plot):
        if labels[i] == key:
            data = np.loadtxt(f)
            ax.plot(
                data[:, 0], data[:, 1]*scale_factor, 
                color=colors[labels[i]], 
                label="along " + labels[i] + "-axis"
            )

ax.legend()
ax.set_xlabel("$\mathrm{\u03c9}$ [eV]")
ax.set_ylabel("Im[$\mathrm{\u03B5_M(\u03c9)}$]")
ax.set_xlim(1, 7)
plt.savefig("bse.png")
