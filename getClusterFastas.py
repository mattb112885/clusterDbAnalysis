#!/usr/bin/python

# Will generate a FASTA file for every different cluster
# in the results of "db_getClusterGeneInfo.py" ("infofile").
# The results are placed in "outputfolder"
#
# Files outputted with first 10 digits of run ID followed by the cluster ID
#
##############
# WARNING!!!!!
##############
# Cluster ID MUST BE SORTED in the input file (this shouldn't be a problem
# but if something weird happens it isn't my fault!!)
#

import sys, os

if not len(sys.argv) == 3:
    print "Usage: ./getCoreClusterFastas.py [infofile] [outputfolder]"

infofile = sys.argv[1]
outputfolder = sys.argv[2]

s = open(infofile, "r")

lastone = None
fid = -1
for line in s:
    spl = line.strip().split("\t")
    myrunid = spl[0]
    if not spl[1] == lastone:
        if not fid == -1:
            fid.close()
        fid = open(os.path.join(outputfolder, myrunid[0:9] + "_" + spl[1] + ".fasta"), "w")
        lastone = spl[1]
    
    fid.write(">" + spl[2] + " " + spl[4] + "\n")
    fid.write(spl[5] + "\n")
    
