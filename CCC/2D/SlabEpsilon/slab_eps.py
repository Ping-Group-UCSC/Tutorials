#!/usr/bin/env python3

import os
import sys
import numpy as np
from scipy.integrate import simps
import matplotlib.pyplot as plt



def readInpFile():
    '''
    read input file returning read values as tuple
    '''
    # default calculation input (ci)
    ci = {
        'path2qe'       : "QE",
        'path2jdftx'    : "JDFTx",
        'outfile'       : "epsilon_z.dat"
    }
    # read inpFile and assign values in ci
    try:
        with open(sys.argv[1]) as f:
            for line in f:
                tLine = line.split('#')[0]
                try:
                    key, value = tuple( map(str.strip, tLine.split('=')) )
                    ci[key] =  value
                except:
                    None    # ignore all lines without a clear assignment statement
    except IndexError:
        None # just use the defaults

    print("Input read as: ", ci)
    return tuple(ci.values())


def readDielectricMatrix(QE, dynmat_file="dynmat.out"):
    '''
    read dielectric matrix from QE output
    '''
    dynmat_file = os.path.join(QE, dynmat_file)

    with open(dynmat_file) as f:
        lines  = iter(f.readlines())

    for line in lines:
        if "... with zone-center polar mode contributions" in line:
            return np.array([ next(lines).split()[0:3] for _ in range(3) ], dtype=np.float64)
    
    sys.stderr.write("dielectric matrix not found in {}\n".format(dynmat_file))
    
    return None


def readDielectricPerp(JDFTx, dump_file="plus_dump.slabEpsilon"):
    '''
    read perpendicular dielectric profile from JDFTx output
    '''
    return np.loadtxt(os.path.join(JDFTx, dump_file))


def calcDielectricPara(eps_z_perp, eps_mat):
    '''
    calculate parallel dielectric profile from the perpendicular component and the dielectric matrix
    '''
    # average value of eps_perp(z)
    eps_avg_perp = simps(eps_z_perp[:,1],eps_z_perp[:,0])/max(eps_z_perp[:,0])

    # averaged parallel components of the dielectric matrix
    eps_mat_para = (eps_mat[0,0] + eps_mat[1,1]) / 2

    # prefactor for eps_z_para
    prefactor = (eps_mat_para - 1) / ( eps_avg_perp - 1)

    # return eps_z_para = prefactor * (eps_z_perp - 1) + 1
    return prefactor * (eps_z_perp[:,1] - 1) + 1


def eps_plot(eps_z, plt_file='epsilon_z.eps'):
    '''
    plot dielectric function
    '''
    print("Plotting data to {}".format(plt_file))

    plt.plot(eps_z[:,0], eps_z[:,1], 'o-')
    plt.plot(eps_z[:,0], eps_z[:,2], 'o-')
    plt.xlabel(r'$z$')
    plt.ylabel(r'$\epsilon (z)$')
    plt.savefig(plt_file, format='eps')
    plt.close()
    
    return None


def slab_eps(lprint=False, lplot=True):
    '''
    This program is used to calculate the dielectric profile from JDFTx and QE calculations
    '''

    help_message = "\
    This program is used to calculate the dielectric profile from JDFTx and QE calculations.\n\
        Usage: {0} [options] <input file>\n\
        Options:\n\
            -h | --help     display this help menu and quit\n\
    ".format(os.path.basename(sys.argv[0]))
    
    '''
    Handle wrong usage or help
    '''
    if "-h" in sys.argv or "--help" in sys.argv:
        print(help_message)
        sys.exit(0)

    '''
    Begin code
    '''

    # read input file
    try:
        path2qe, path2jdftx, outfile = readInpFile()
    except:
        path2qe = "QE"
        path2jdftx = "JDFTx"
        outfile = "epsilon_z.dat"
    
    # read QE dielectric matrix
    eps_mat = readDielectricMatrix(path2qe)
    if lprint:
        print("Dielectric Matrix read as:")
        for eps in eps_mat:
            print("    {:16.9f}  {:16.9f}  {:16.9f}".format(eps[0], eps[1], eps[2]))

    # read dielectric profile (perpendicular)
    eps_z_perp = readDielectricPerp(path2jdftx)
    if lprint:
        print("Dielectric profile (perpendicular) read as:")
        for i, eps in enumerate(eps_z_perp):
            print("    {:16.9f}  {:16.9f}".format(eps[0], eps[1]))
    
    # calculate dielectric profile (parallel)
    eps_z_para = calcDielectricPara(eps_z_perp, eps_mat)
    if lprint:
        print("Dielectric profile (parallel) read as:")
        for i, eps in enumerate(eps_z_para):
            print("    {:16.9f}".format(eps))


    # combine eps_z_perp and eps_z_para to one array
    eps_z = np.stack((eps_z_perp[:,0], eps_z_perp[:,1], eps_z_para), axis=-1)
    
    # save to txt
    print("Writing data to {}".format(outfile))
    np.savetxt(outfile, eps_z, \
        fmt='%.18e', header='{:23s}{:25s}{:25s}'.format('z','e_perp','e_para'))
    
    # create plot
    if lplot:
        eps_plot(eps_z)
    

    return eps_z


if __name__ == "__main__":
    slab_eps()