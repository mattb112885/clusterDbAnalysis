#!/usr/bin/env python

# This is a pipe command.
# Pipe in the "geneinfo" file.
#
# Will generate a FASTA file for every different cluster
# in the results of "db_getClusterGeneInformation.py" ("infofile").
# The fasta files are placed in "outputfolder"
#
# Files outputted with run ID followed by the cluster ID
#
#

import sys, os
import fileinput
import optparse

usage = "%prog [options] fastafolder < cluster_info_file"
description="""Generate a FASTA file for each cluster present in the specified clusterinfo file 
(as generated from e.g. db_getClusterGeneInformation.py)

If cluster info is not available (i.e. the file was generated from db_getGeneInformation.py),
it will generate a fasta file with the name "NOCLUSTER_.fasta" with all of the sequences
in the file. """
parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-n", "--nucleotides", help="Export nucleotide fasta files, not protein (D: Protein fasta files)", action="store_true", dest="nuc", default=False)
(options, args) = parser.parse_args()

if not len(args) == 1:
    sys.stderr.write("ERROR: Fastafolder is a required argument (see -h for usage details)\n")
    exit(2)

outputfolder = args[0]

# Make it not fail if the directory doesn't exist already.
if not os.path.isdir(outputfolder):
    if os.path.exists(outputfolder):
        sys.stderr.write("ERROR: Specified output folder name already exists as an ordinary file\n")
        exit(2)
    os.mkdir(outputfolder)

# Make a dictionary that will hold our outfile names and handles
# This is used to remove any dependency on order of rows in the file and also
# make it non-ambiguous if multiple clusters have the same cluster IDs but different run IDs
outfiles = {}
for line in fileinput.input("-"):
    spl = line.strip('\r\n').split("\t")
    if len(spl) < 14:
        myrunid="NOCLUSTER"
        myclustid = ""
    else:
        myrunid = spl[12]
        myclustid = spl[13]
    fname = "%s_%s.fasta" %(myrunid, myclustid)
    #if we haven't seen this cluster run and ID yet, add it's handle to out dictionary, and make it the current fid
    if fname not in outfiles.keys(): 
        fid = open(os.path.join(outputfolder, fname), "w")
        outfiles[fname] = fid
    if options.nuc:
        mysequence = spl[10]
    else:
        mysequence = spl[11]
    mygeneid = spl[0]
    myannotation = spl[9]
    ln = ">%s %s\n%s\n" %(mygeneid, myannotation, mysequence)
    outfiles[fname].write(ln)

#close all files
for fhandle in outfiles.values():
    fhandle.close()


