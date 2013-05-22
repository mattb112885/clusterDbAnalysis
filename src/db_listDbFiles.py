#!/usr/bin/python

from FileLocator import *
import os, optparse, sys

usage = "%prog [options] [searchstring] > filelist"
description="List all files in the ITEP directories provided as part of this software package"
parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-w", "--maxw", 
                  help="Maximum number of characters wide (D=Print each file on its own line)", 
                  action="store", type="int", dest="maxw", default=1)
(options, args) = parser.parse_args()

subpaths = ['src','scripts','src/utilities']
paths = [os.path.join(locateRootDirectory(), sub) for sub in subpaths]

def printcol(ls):
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

for path in paths:
    ls = os.listdir(path)
    if args==[]:
        printcol(ls)
    else:
        query= set([q.lower() for q in args])
        def hits(h): return any([h.lower().find(q)>-1 for q in query])
        ls = [h for h in ls if hits(h)]
        if len(ls) > 0:
            print('\nPrograms found in '+path+':')
            printcol(["\t"+l for l in ls])
        else:
            print('\nNo programs found in '+path+'.')
print("\n")

