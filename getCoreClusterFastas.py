#!/usr/bin/python

# Will generate a FASTA file for every core gene...
# I'm a fucking moron. I'm sorry this took so long.
#
# Requres an info file (below)
# Files outputted with first 10 digits of run ID followed by the cluster ID
# Cluster ID MUST BE SORTED in the input file (this shouldn't be a problem
# but if something weird happens it isn't my fault!!)
#

import sys
infofile = sys.argv[1]

s = open(infofile, "r")

lastone = None
fid = -1
for line in s:
    spl = line.strip().split("\t")
    if not spl[1] == lastone:
        if not fid == -1:
            fid.close()
        fid = open(spl[0][1:10] + "_" + spl[1] + ".fasta", "w")
        lastone = spl[1]
    
    fid.write(">" + spl[2] + " " + spl[4] + "\n")
    fid.write(spl[5] + "\n")
    
