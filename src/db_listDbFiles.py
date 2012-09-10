#!/usr/bin/python

import os, optparse
usage = "%prog > filelist"
description="List all executable files provided as part of this software package"
parser = optparse.OptionParser(usage=usage, description=description)
parser.parse_args()

fileloc = os.path.abspath(__file__)
path = os.path.dirname(fileloc)
ls = os.listdir(path)

ln = "\n".join(s for s in sorted(ls))
print ln
