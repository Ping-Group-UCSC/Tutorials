#!/usr/bin/env python

import numpy as np
import pandas as pd
import glob
import matplotlib.pyplot as plt
import sys


def read_socme_as_df(filepath):
    ''' read socme from orca tddft output to pandas dataframe '''
    with open(filepath) as f:
        # look for SOCME in file
        for line in f:
            if 'SOCME' in line:
                for _ in range(4):
                    next(f)
                socme = []
                while True:
                    line = next(f)
                    if not line.strip():
                        break
                    split_line = np.array(line.split())[[0, 1, 5, 10, 15]]
                    socme.append(split_line)
                # return df
                socme_df = pd.DataFrame(socme)
                socme_df.columns = ['T', 'S', 'Z', 'X', 'Y']
                return socme_df
    raise ValueError('SOCME not found in filepath: {}'.format(filepath))


def read_zfs_evec(filepath):
    with open(filepath) as f:
        for line in f:
            # look for diagonalized in file
            if 'diagonalized' in line:
                for _ in range(2):
                    next(f)
                evec = []
                while True:
                    line = next(f)
                    if not line.strip():
                        break
                    evec.append(line.split())
                # return numpy array
                evec = np.array(evec, dtype=float)
                return evec
    raise ValueError('diagonalized not found in filepath: {}'.format(filepath))


def find_index(socme_df, T, S):
    filt = (socme_df['T'] == T) & (socme_df['S'] == S)
    return socme_df[filt].index[0]


def df_dict_to_array(socme_df_dict, T, S):
    ''' return numpy array from socme_df_dict for particular T and S '''
    T_S_index = find_index(socme_df_dict['0.0000'], T, S)
    socme_array = np.array([
        socme_df_dict[theta].loc[T_S_index, ['X', 'Y', 'Z']]
        for theta in socme_df_dict
    ], dtype=float)
    return socme_array


def plot_soc(angle_name, T, S):
    fname = 'socme_{}.xlsx'.format(angle_name)
    with pd.ExcelWriter(fname) as writer:
        socme_df_dict = {}
        for angle_folder in sorted(glob.glob('calcs/{}_*'.format(angle_name))):
            angle = angle_folder.split('_')[-1]
            filepath = angle_folder + '/tddft.out'
            socme_df = read_socme_as_df(filepath)
            socme_df.to_excel(writer, sheet_name=angle, index=False)
            socme_df_dict[angle] = socme_df
    angles = np.array(list(socme_df_dict.keys()), dtype=float)
    socme_array = df_dict_to_array(socme_df_dict, T, S)
    labels = ['X', 'Y', 'Z']
    plt.style.use('seaborn-talk')
    for i in range(3):
        plt.scatter(angles, np.abs(socme_array[:, i]), label=labels[i])
    plt.legend()
    plt.ylabel(r'SOC (cm$^{{-1}}$) [T={},S={}]'.format(T, S))
    plt.xlabel(r'\{}'.format(angle_name))
    plt.xlim(0, np.pi)
    plt.xticks(
        ticks=[0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi], labels=['0', r'$\pi/4$', r'$\pi/2$', r'$3\pi/4$', r'$\pi$']
    )
    plt.savefig('socme_{}.png'.format(angle_name))
    plt.close()


def plot_zfs(angle_name, ivec=0):
    angles = []
    evecs = []
    for angle_folder in sorted(glob.glob('calcs/{}_*'.format(angle_name))):
        angle = angle_folder.split('_')[-1]
        filepath = angle_folder + '/zfs.out'
        evec = read_zfs_evec(filepath)
        evecs.append(evec[ivec])
        angles.append(angle)
    angles = np.array(angles, dtype=float)
    evecs = np.array(evecs, dtype=float)
    labels = ['X', 'Y', 'Z']
    plt.style.use('seaborn-talk')
    print(evecs)
    print()
    print(evecs[:, 0])
    for i in range(3):
        plt.scatter(angles, np.abs(evecs[:, i]), label=labels[i])
    plt.legend()
    plt.ylabel(r'Vec')
    plt.xlabel(r'\{}'.format(angle_name))
    plt.xlim(0, np.pi)
    plt.xticks(
        ticks=[0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi], labels=['0', r'$\pi/4$', r'$\pi/2$', r'$3\pi/4$', r'$\pi$']
    )
    plt.savefig('zfs_{}.png'.format(angle_name))
    plt.close()


if __name__ == '__main__':
    # angle_name = 'phi'
    angle_name = 'theta'
    try:
        calc = sys.argv[1]
    except IndexError:
        print("Please provide calc: 'zfs' or 'soc'")
        sys.exit(0)
    assert calc in ['zfs', 'z', 'soc', 's'], "Invalid calc: {}".format(calc)
    if calc == 'zfs' or calc == 'z':
        plot_zfs(angle_name)
    elif calc == 'soc' or calc == 's':
        plot_soc(angle_name, '2', '1')
