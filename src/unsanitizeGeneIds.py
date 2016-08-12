#!/usr/bin/env python

from __future__ import print_function
import fileinput, optparse, re

usage = "%prog < infile > outfile"
description="""Replace the gene IDs in infile from fig_\d+_\d+_peg_\d+ to fig|\d+.\d+.peg.\d+ 
(the former comes from some of the sanitation scripts)"""
parser = optparse.OptionParser(usage=usage, description=description)
parser.parse_args()

for line in fileinput.input("-"):
    line = line.strip("\r\n")
    # The r is needed becuase otherwise the python interpreter interprets \1 as 1 ...
    ln = re.sub("fig_(\d+)_(\d+)_peg_(\d+)", r"fig|\1.\2.peg.\3", line)
    print(ln)
