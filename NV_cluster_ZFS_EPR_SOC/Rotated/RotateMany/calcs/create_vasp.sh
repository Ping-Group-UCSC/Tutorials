#!/usr/bin/env bash

for d in */ ; do
    cd $d
    printf "$d ... "
    conv_geo.py vasp nv.xyz > nv.vasp
    echo "done"
    cd ..
done
