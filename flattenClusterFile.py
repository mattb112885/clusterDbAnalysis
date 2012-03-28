#!/usr/bin/python

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

parser = optparse.OptionParser()
parser.add_option("-I", "--inflation", help="Inflation parameter used for MCL clustering", action="store", type="float", dest="inflation", default=None)
parser.add_option("-c", "--cutoff", help="Scoring cutoff used for MCL clustering", action="store", type="float", dest="cutoff", default=None)
parser.add_option("-m", "--method", help="Scoring method used for MCL clustering", action="store", type="string", dest="method", default=None)
(options, args) = parser.parse_args()

if options.method == None or options.cutoff == None or options.inflation == None:
    sys.stderr.write("ERROR in flattenClusterFile: -I, -c, and -m are all required")
    exit(2)

runid = "I_%1.2f_cutoff_%1.2f_method_%s" %(options.inflation, options.cutoff, options.method)
if len(runid) > 64:
    print sys.stderr.write("WARNING: Resulting clusterID is too long - will truncate to 64 characters")
    runid = runid[0:64]

rowIdx = 1
for line in fileinput.input("-"):
    spl = line.strip().split('\t')
    for gene in spl:
        string = str(runid) + "\t" + str(rowIdx) + "\t" + gene
        print string
    rowIdx = rowIdx + 1
