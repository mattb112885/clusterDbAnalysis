#!/usr/bin/python

# This is a pipe command.
# Pipe in the "geneinfo" file.
#
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
import fileinput
import optparse

usage = "%prog [fastafolder] < cluster_info_file"
description="Generate a FASTA file for each cluster present in the specified clusterinfo file (as generated from e.g. db_getClusterGeneInfo). The file must be in order by cluster"
parser = optparse.OptionParser(usage=usage, description=description)
(options, args) = parser.parse_args()

if not len(args) == 1:
    sys.stderr.write("ERROR: Fastafolder is a required argument (see -h for usage details)\n")
    exit(2)

outputfolder = args[0]

lastone = None
fid = -1
for line in fileinput.input("-"):
    spl = line.strip('\r\n').split("\t")
    myrunid = spl[0]
    if not spl[1] == lastone:
        if not fid == -1:
            fid.close()
        fname = "%s_%s.fasta" %(myrunid, spl[1])
        fid = open(os.path.join(outputfolder, fname), "w")
        lastone = spl[1]
    # Write to fasta
    ln = ">%s %s\n%s\n" %(spl[2], spl[4], spl[5])
    fid.write(ln)

