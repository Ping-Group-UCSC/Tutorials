#!/usr/bin/env python

import numpy as np


def read_zfs_evecs(filename):
    with open(filename) as f:
        for line in f:
            if line.startswith('diagonalized'):
                for _ in range(2):
                    next(f)
                evecs = np.array([
                    np.fromstring(next(f), sep=' ', dtype=float)
                    for _ in range(3)
                ])
                return evecs.T
    raise ValueError("Did not find keyworkd 'diagonalized' in filename: {}".format(filename))


if __name__ == '__main__':
    evecs = read_zfs_evecs('zfs.out')
    print(evecs)
