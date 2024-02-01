#!/usr/bin/env python

import pw2py as pw
import numpy as np
from scipy.spatial.transform import Rotation as R
from copy import deepcopy
import os
import shutil


xyz_file = 'nv.xyz'
geo = pw.atomgeo.from_file(os.path.join('common', xyz_file))
og = deepcopy(geo)

calc_folder = 'calcs'
if not os.path.exists(calc_folder):
    os.mkdir(calc_folder)


def make_rotated_calcs(axis):
    if axis == 'z':
        angle_name = 'theta'
    elif axis == 'y':
        angle_name = 'phi'
    else:
        raise ValueError("Axis must be 'y' or 'z', passed: {}".format(axis))
    for angle in np.linspace(0, np.pi, num=33):
        print('{} = {:.4f} ... '.format(angle_name, angle), end='')
        # mkdir
        folder = os.path.join(
            calc_folder, '{}_{:.4f}'.format(angle_name, angle)
        )
        if not os.path.exists(folder):
            os.mkdir(folder)
            shutil.copy(os.path.join('common', 'tddft.in'), folder)
            shutil.copy(os.path.join('common', 'zfs.in'), folder)
            # calc geo
            r = R.from_euler(axis, angle).as_matrix()
            geo.pos = r.dot(og.pos.T).T
            # write file
            geo.write_file(os.path.join(folder, xyz_file))
            print('done')
        else:
            print('skipped')


make_rotated_calcs('z')
