#!/usr/bin/python

# Because Unix's join function is picky EVEN THOUGH I SORTED THE FUCKIGN ROWS
# I'm writing this instead.

import optparse
import sys

usage="%prog -1 [col1] -2 [col2] (options) File1 File2"
description = "To get around some of the obnoxiousness of the UNIX join command..."

parser=optparse.OptionParser(usage=usage, description=description)
parser.add_option("-1", help="Column in first file starting from 1 (required)", action="store", type="int", dest="firstcol", default=None)
parser.add_option("-2", help="Column in second file starting from 1 (required)", action="store", type="int", dest="secondcol", default=None)
parser.add_option("-d", help="Delimiter (default: Tab)", action="store", type="str", dest="delimiter", default="\t")
parser.add_option("--k1", help="Keep non-matching rows in the first file (D: Only print matching rows)", action="store_true", dest="keep1", default=False)
parser.add_option("--k2", help="Keep non-matching rows in the second file (D: Only print matching rows)", action="store_true", dest="keep2", default=False)
(options, args) = parser.parse_args()

if options.firstcol is None or options.secondcol is None:
    sys.stderr.write("ERROR: Both -1 and -2 arguments are required\n")
    exit(2)
if len(args) != 2:
    sys.stderr.write("Must provide exactly two file names\n")
    exit(2)

firstcol = options.firstcol - 1
secondcol = options.secondcol - 1

firstKeyToRows = {}
firstFileNumCols = 0
for line in open(args[0], "r"):
    spl = line.strip("\r\n").split(options.delimiter)
    firstFileNumCols = len(spl)
    key = spl[firstcol]
    if key in firstKeyToRows:
        firstKeyToRows[key].append(spl)
    else:
        firstKeyToRows[key] = [ spl ]

secondKeyToRows = {}
secondFileNumCols = 0
for line in open(args[1], "r"):
    spl = line.strip("\r\n").split(options.delimiter)
    secondFileNumCols = len(spl)
    key = spl[secondcol]
    if key in secondKeyToRows:
        secondKeyToRows[key].append(spl)
    else:
        secondKeyToRows[key] = [ spl ]

for key in firstKeyToRows:
    if key not in secondKeyToRows:
        if options.keep1:
            for spl in firstKeyToRows[key]:
                st = "\t".join(spl) + "\t"*secondFileNumCols
                print st
    else:
        for spl in firstKeyToRows[key]:
            for spl2 in secondKeyToRows[key]:
                st = "\t".join(spl) + "\t" + "\t".join(spl2)
                print st

for key in secondKeyToRows:
    if key not in firstKeyToRows:
        if options.keep2:
            for spl2 in secondKeyToRows[key]:
                st = "\t"*firstFileNumCols + "\t".join(spl2)
                print st
