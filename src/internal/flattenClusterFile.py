#!/usr/bin/env python

# This is a pipe command.
# 
# Pipe in the cluster file (from MCL)
# and the flattened file comes out in stdout
# You must provide this function what parameters you used to generate
# the clusters:
#
# - inflation
# - distance metric
# - cutoff
#
# All of these will go into the run ID.
#
# The "flattened file" contains only three columns:
# - A run ID (which should be unique)
# - Cluster number 
# - Gene ID
#
# This will make it easy to generate SQL queries for particular clusters
# (e.g. give me the blast results for all of the genes listed in cluster 1
# against anything else)

import fileinput
import optparse
import sys

usage = "%prog -n runname < MCL_cluster_file > flattened_file"
description = "Turn MCL cluster file into a more database-friendly format, the run is given its own ID (input to this function) and each cluster within it is given a cluster ID"
parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-n", "--name", help="Run name (should have descriptions of Inflation, method and cutoff)", action="store", type="string", dest="name", default=None)
(options, args) = parser.parse_args()

if options.name == None:
    sys.stderr.write("ERROR in flattenClusterFile: -n is required")
    exit(2)

runid = options.name
if len(runid) > 255:
    sys.stderr.write("WARNING: Resulting clusterID is too long - will truncate to 255 characters\n")
    runid = runid[0:255]

rowIdx = 1
for line in fileinput.input("-"):
    spl = line.strip('\r\n').split('\t')
    for gene in spl:
        string = str(runid) + "\t" + str(rowIdx) + "\t" + gene
        print string
    rowIdx = rowIdx + 1
