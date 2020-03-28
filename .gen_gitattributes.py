#!/usr/bin/env python3

import os
import sys
import subprocess as sp

skip = ["py", "html", "gnu", "md", "sh", "git"]

out = sp.check_output("find . -name '*.*'", shell=True).splitlines()

list_suffix = []

for v in out:
    suffix = str(v).split('.')[-1].split("'")[0]
    if suffix not in list_suffix:
        list_suffix.append(suffix)

for suffix in list_suffix:
    if suffix not in skip + [""] :
        print("*.{} -text".format(suffix))
