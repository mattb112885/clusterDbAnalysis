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

usage = "%prog [options] < gene_ids > clusters_containing_genes"
description = "Given a list of gene IDs, gets a list of clusters containing those genes (in all run IDs)"
parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-g", "--gcolumn", help="Column number (start from 1) for gene ID", action="store", type="int", dest="genecolumn", default=1)
parser.add_option("-r", "--runid", help="Restrict results to the given run ID", action="store", type="str", dest="runid", default=None)
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
for tup in res:
    print("\t".join(tup))

con.close()
