#!/usr/bin/env python3

import os
import subprocess as sp

skip = ["", ".py", ".html", ".gnu", ".md", ".sh", ".git"]

filenames = sp.check_output("find . -name '*.*'", shell=True).decode('ascii').splitlines()

extensions = []
for filename in filenames:
    _, extension = os.path.splitext(filename)
    if extension not in extensions + skip:
        extensions.append(extension)

# generate .gitattributes file
print("$ cat .gitattributes")
for extension in extensions:
    print("*{} linguist-documentation".format(extension))
