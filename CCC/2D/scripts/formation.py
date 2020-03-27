#!/usr/bin/env python3

import json
import mendeleev as melv
import numpy as np
import os
import pandas as pd
import re
import signal
import sys

# # personal modules
# import tjs_resource as tjs
# import pw2py as pw

ry2ev = 13.605662285137
ha2ev = 2 * ry2ev


def compoundToDict(compound):
    '''
    convert supplied compound to dictionary
    '''
    compound = compound.split()
    out = {}
    for i in range(0, len(compound) // 2 + 2, 2):
        out[compound[i]] = int(compound[i+1])
    
    return out


def grepFinal(filename, scale=0, search="Final"):
    '''
    grep final energy from filename, scale result by scale/nat if scale != 0
    '''
    with open(filename) as f:
        for line in f:
            if "number of atoms/cell" in line:
                nat = int(line.split()[-1])
            elif search in line:
                if scale == 0:
                    return float(line.split()[-2]) * ry2ev
                return float(line.split()[-2]) / nat * scale * ry2ev
    
    # tjs.warn("No '{}' energy in '{}'".format(search, filename))
    return None


def grepCCC(filename):
    '''
    grep charge cell correction from filename
    '''
    with open(filename) as f:
        for line in f:
            if "Net" in line:
                return float(line.split()[2]) * ha2ev

    # tjs.warn("No Net CCC in {}".format(filename))
    return None


def gatherPrist():
    '''
    return dictionary object with total energy(E), valence band maximum (VBM), 
        conduction band minimum (CBM), and band gap (GAP)
    '''
    # prist
    prist = {}
    prist['E'] = grepFinal("Prist/relax.out")
    prist['VBM'] = -3.6165
    prist['CBM'] = 1.0890
    prist['GAP'] =  prist['CBM'] - prist['VBM']

    # # chem
    # chem = {}
    # for ion in ['BN']:
    #     chem[ion] = {}
    # return prist, chem

    return prist


def FE(FE0, q, ef):
    '''
    FE_X(q,ef) = FE_X(q,ef=0) + q * ef
    '''
    return FE0 + q * ef


def formationLine(data, prist):
    '''
    calculation formation Line data (return numpy array)
    '''
    minCharge = min( map(int, data.keys()) )
    maxCharge = max( map(int, data.keys()) )

    line = []
    ie = {}
    ie['CBM'], ie['VBM'] = {}, {}
    line.append(np.array([0.0, data[str(maxCharge)]]))
    for charge in range(maxCharge, minCharge, -1):
        ef = data[str(charge - 1)] - data[str(charge)]
        line.append([ef, FE(data[str(charge)], charge, ef)])
        ie['CBM']["{}/{}".format(str(charge - 1), str(charge))] = prist['GAP'] - ef
        ie['VBM']["{}/{}".format(str(charge - 1), str(charge))] = ef
    
    line.append(np.array([prist['GAP'], FE(data[str(minCharge)], minCharge, ef)]))

    return {'line' : np.array(line), 'CBM' : ie['CBM'], 'VBM' : ie['VBM'] }


def readInpFile(inpFile):
    '''
    read input file returning read values as tuple
    '''
    # default calculation input (ci)
    ci = {
        'dopants'       : None,
        'conditions'    : None
    }
    # read inpFile and assign values in ci
    checkFile(inpFile)
    with open(inpFile) as f:
        for line in f:
            tLine = line.split('#')[0]
            try:
                key, value = tuple( map(str.strip, tLine.split('=')) )
                if key not in ci.keys():
                    tjs.die("In input file, key '{}' not recognized".format(key))
                ci[key] =  value
            except:
                None    # ignore all lines without a clear assignment statement

    print("Input read as: ", ci)
    return tuple(ci.values())


def formation():
    '''
    Main program, collect total energies and generate data folder w/ data files and gnuplot scripts
    '''

    help_message = "\
    This program collect total energies and generate data folder w/ data files and gnuplot scripts.\n\
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

    dumpDir = "FormData"
    if not os.path.exists(dumpDir):
        os.mkdir(dumpDir)

    # read input file
    try:
        dopants, conditions = readInpFile(sys.argv[1])
    except:
        dopants = "Ti"
        # conditions = ["rich", "poor"]

    dopants = dopants.split()

    # prist, chem = gatherPrist()
    prist = gatherPrist()

    # # gather chemical potentials
    # for ion in (dopants):
    #     chem[ion] = {}
    #     chem[ion]['poor'] = grepFinal(os.path.join(ion, \
    #         "ChemPot/Oxide-{}O2/vc.out".format(ion)), scale=3) - 2 * chem['O']['poor']
    #     chem[ion]['rich'] = grepFinal(os.path.join(ion, \
    #         "ChemPot/Oxide-{}O2/vc.out".format(ion)), scale=3) - 2 * chem['O']['rich']


    '''
    Calculate formation energy at ef=0 (VBM)
        Ef_X(q,ef=0) = E_X(q) - E_prist + mu_Fe - mu_X + delta_q + q * vbm
        prist + X -> doped + Fe
    '''
    # large data structure for all dopants and calcs
    data = {}
    # for condition in conditions:
    #     data[condition] = {}
    for dopant in dopants:
        data[dopant] = {}

        '''
        This part may be moved to a funciton as it can be considered seperate
        '''
        # gather total energies (w/ CCC when charged)
        for path in os.listdir(dopant):
            charge = str(int(re.sub('Q', '', path)))
            # total energy
            data[dopant][charge] = grepFinal(os.path.join(dopant, path, "relax.out"))
            # charge cell correction
            if int(charge) != 0:
                data[dopant][charge] += grepCCC(os.path.join(dopant, path, "CCC/completed/main.out"))
            # shift by q * vbm
            data[dopant][charge] += float(charge) * prist['VBM']
            # shift by prist energy
            data[dopant][charge] -= prist['E']
            # shift by chempot
            # data[dopant][calc][charge] += chem['BN'][condition] - chem[dopant][condition]

        
        # reorder keys
        for k in range(min( map(int, data[dopant].keys() ) ), \
            max( map(int, data[dopant].keys() ) ) + 1 ):
            data[dopant][str(k)] = data[dopant].pop(str(k))



    '''
    Calculate formation energy data: line and ionization energies (with different references)
    '''
    formData = {}
    for dopant in dopants:
        formData[dopant] = formationLine(data[dopant], prist)
        # save data
        np.savetxt(os.path.join(dumpDir, "fe_{}.dat".format(dopant)), \
            formData[dopant]['line'], delimiter=' ', \
                fmt='%19.12f', header="{:19s}{:19s}".format("e_f(eV)", "FE(eV)"))
            


    return prist, data, chem, formData


if __name__ == "__main__":
    formation()
