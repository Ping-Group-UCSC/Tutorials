#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt
from parse_soc import df_dict_to_array
import numpy as np


xls = pd.ExcelFile('socme_theta.xlsx')
df_dict = {
    sheet_name: xls.parse(sheet_name)
    for sheet_name in xls.sheet_names
}
angles = np.array(xls.sheet_names, dtype=float)

df_simple = pd.DataFrame(df_dict_to_array(df_dict, 2, 1))
df_simple.columns = ['X', 'Y', 'Z']
df_simple['perp'] = (df_simple['X']**2 + df_simple['Y']**2)**0.5 * (1/2)**0.5

plt.style.use('seaborn-talk')
for col in ['X', 'Y', 'Z']:
    plt.scatter(angles, np.abs(df_simple[col]), label=col)
plt.scatter(angles, df_simple['perp'], label=r'$\perp$')
plt.legend()
plt.ylabel(r'SOC (cm$^{{-1}}$) [T={},S={}]'.format(2, 1))
plt.xlabel(r'$\{}$'.format('theta'))
plt.xlim(0, np.pi)
plt.xticks(
    ticks=[0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi], labels=['0', r'$\pi/4$', r'$\pi/2$', r'$3\pi/4$', r'$\pi$']
)
plt.savefig('compare.png')
plt.close()


# print(df_dict['0.0000']['T'].dtype)
# def df_dict_to_array(socme_df_dict, T, S):
