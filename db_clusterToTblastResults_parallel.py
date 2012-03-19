#!/usr/bin/python

# This is a pipe command. Pipe in the run ID(s) you wish to use to run the command.
# 
# Command line options:
# -n --ncores : Number of cores. Default: 1
# -r --rcolumn : Run ID column in stdin (count from 1). Default: 1
#
# Prints the results for particular clustter to files in the "tblastn_by_cluster" folder
#
# Performs the same function as db_clusterToTblastResults (actually, it calls that function with different parameters
# in parllel)
#
# However, Rather than lumping all of the results into one giant file, this function will instead
# generate separate files for each of the separate numbers of clusters from 1 to the specified number
# and call db_clusterToTblastResults separately for each of them.
#
# Unlike the direct db_clusterToTblastResults call, this one is intended to be used when you want ALL
# of the possible combinations at a time, not just those for a particular subset of interesting
# clusters. It is much faster than the monolothic version.
#
# Note - sqlite is threadsafe so making the temporary table doesn't crash this.

import os, fileinput, optparse
from ruffus import *

# Max number of genes in a cluster (note - there isn't too much of a time penalty if this is an overestimate
# but if it is an underestimate we could miss the biggest clusters)
maxn = 1000

# Get input arguments
usage = "%prog [options]"
description="Given a set of run IDs for which to extract TBLAST results (from stdin), extract TBLASTN results within each cluster (break the task up into parallel parts for speed). Resutls are stored in the tblastn_by_cluster folder"
parser = optparse.OptionParser(usage=usage)
parser.add_option("-r", "--rcolumn", help="Column number (start from 1) for run ID (defualt = 1)", action="store", type="int", dest="runcolumn", default=1)
parser.add_option("-n", "--ncores", help="Number of cores to use in processing (default = 1)", action="store", type="int", dest="NCORES", default=1)
(options, args) = parser.parse_args()

rc = options.runcolumn - 1 # Convert to Pythonic indexes                                                                                                                                                      
NUMCORES = options.NCORES

# Read run IDs from stdin
s = []
for line in fileinput.input("-"):
    ln = line.strip().split("\t")
    s.append(ln[0])
unqruns = set(s) # Unique run IDs

# Set up command to run in parallel (basically, get the clusters with 1 gene, then 2, then ... and then do each one of those separately
# and save them each to different files).
params = []
for runid  in unqruns:
    for i in range(1, maxn):
        fname = "tblastn_by_cluster/RUN_" + runid[1:10] + "_N_" + str(i) + "_tblastn"
        cmd = "echo \"" + runid + "\" | python src/db_getClustersWithNumGenes.py " + str(i) + " | python src/db_clusterToTblastResults.py > " + fname
        params.append([cmd,])

# Setup for ruffus
@parallel(params)
def runcmd(command):
    print command
    os.system(command)
    return

# Run in parallel
pipeline_run([runcmd], multiprocess=int(NUMCORES))
