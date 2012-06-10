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

usage = "%prog [fastafolder] [labelfolder] < cluster_info_file"
description="Generate a FASTA file for each cluster present in the specified clusterinfo file (as generated from e.g. db_getClusterGeneInfo). The file must be in order by cluster"
parser = optparse.OptionParser(usage=usage, description=description)
(options, args) = parser.parse_args()

if not len(args) == 2:
    sys.stderr.write("ERROR: Both fastafolder and labelfolder are required (see -h for usage details)\n")
    exit(2)

outputfolder = args[0]
labelfolder = args[1]

lastone = None
fid = -1
labelfid = -1
for line in fileinput.input("-"):
    spl = line.strip('\r\n').split("\t")
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
