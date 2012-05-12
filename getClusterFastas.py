#!/usr/bin/python

# Will generate a FASTA file for every different cluster
# in the results of "db_getClusterGeneInformation.py" ("infofile").
# The fasta files are placed in "outputfolder"
#
# Organism name and annotation are placed in "labelfolder"
#
# Files outputted with run ID followed by the cluster ID
#
##############
# WARNING!!!!!
##############
# Cluster ID MUST BE SORTED in the input file (this shouldn't be a problem
# but if something weird happens it isn't my fault!!)
#

import sys, os

if not len(sys.argv) == 4:
    print "Usage: ./getCoreClusterFastas.py [infofile] [fastafolder] [labelfolder]"
    exit(2)

infofile = sys.argv[1]
outputfolder = sys.argv[2]
labelfolder = sys.argv[3]

s = open(infofile, "r")

lastone = None
fid = -1
labelfid = -1
for line in s:
    spl = line.strip().split("\t")
    myrunid = spl[0]
    if not spl[1] == lastone:
        if not fid == -1:
            fid.close()
        if not labelfid == -1:
            labelfid.close()
        
        fid = open(os.path.join(outputfolder, myrunid + "_" + spl[1] + ".fasta"), "w")
        labelfid = open(os.path.join(labelfolder, myrunid + "_" + spl[1] + ".fasta_aln_newick.info"), "w")
        lastone = spl[1]
    
    fid.write(">" + spl[2] + " " + spl[4] + "\n")
    fid.write(spl[5] + "\n")
    labelfid.write(spl[2] + "\t" + spl[3] + " " + spl[4] + "\n")
