#!/usr/bin/env python

# This is a pipe command
# Pipe in a tab-delimited file containing gene IDs
#
# Returns all (run ID, clusterID) pairs containing those gene IDs to stdout
#
# Option -g: Column containing gene ID (default = 1), starting from 1 in the first column
#

import sqlite3, fileinput, optparse
from FileLocator import *
from ClusterFuncs import *

header = [ "Run_id", "Cluster_id", "Provided_Gene_id" ]

usage = """%prog [options] < gene_ids > clusters_containing_genes

Output: """ + " ".join(header)
description = """Given a list of gene IDs, gets a list of clusters containing those genes (in all run IDs).
Only returns the genes you specified (not all genes in those clusters)."""
parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-g", "--gcolumn", help="Column number (start from 1) for gene ID", action="store", type="int", dest="genecolumn", default=1)
parser.add_option("-r", "--runid", help="Restrict results to the given run ID", action="store", type="str", dest="runid", default=None)
parser.add_option("--header", help="Specify to add header to the output file (useful if you want to take the results and put into Excel or similar programs)",
                  action="store_true", default=False)
(options, args) = parser.parse_args()
gc = options.genecolumn - 1 # Convert to Pythonic indexes                                                                                                                                                      
con = sqlite3.connect(locateDatabase())
cur = con.cursor()

desiredGenes = []
for line in fileinput.input("-"):
    spl = line.strip('\r\n').split("\t")
    # If the file does not include fig|, go ahead and add it automatically.
    if spl[gc].startswith("fig|"):
        desiredGenes.append(spl[gc])
    else:
        desiredGenes.append("fig|%s" %(spl[gc]))

res = getClustersContainingGenes(desiredGenes, cur, runid=options.runid)
if options.header:
    print "\t".join(header)
for tup in res:
    print "\t".join(tup)

con.close()
