#!/usr/bin/python

# Import the %ID cutoff from the tblastn
#
# This is important to prevent us from having HUGE tables full of garbage
# to join when looking for pseudogenes...

import fileinput, sys

idCutoff = float(sys.argv[1])

for line in fileinput.input("-"):
    spl = line.strip('\r\n').split("\t")
    if(float(spl[2]) > idCutoff):
        print "\t".join(spl)
