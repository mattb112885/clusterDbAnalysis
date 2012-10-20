#!/usr/bin/python

import os, optparse, pprint, sys

usage = "%prog [options] > filelist"
description="List all files in the src/ directory provided as part of this software package"
parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-w", "--maxw", help="Maximum number of characters wide (D=Print each file on its own line)", action="store", type="int", dest="maxw", default=1)
(options, args) = parser.parse_args()

fileloc = os.path.abspath(__file__)
path = os.path.dirname(fileloc)
ls = os.listdir(path)

# Pretty-print into columns
maxw = options.maxw
maxl = max( [ len(s) for s in ls ] )
maxn = int(maxw/maxl)
c = 0
for s in sorted(ls):
    if c >= maxn:
        sys.stdout.write("\n")
        c = 0
    sys.stdout.write(s.ljust(maxl) + "\t")
    c += 1

sys.stdout.write("\n")
