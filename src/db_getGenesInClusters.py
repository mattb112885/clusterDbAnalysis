#!/usr/bin/env python

# This is a pipe command. Pipe in a 2-column table
# containing the run ID in the first column and the
# cluster ID in the second column (DEFAULT) OR
# 
# -r [runcolumn] -c [clustercolumn]
# (starting from 1 as the first column)
# (Note - type db_clusterToTblastResults --help for help details)
#
# Returns a list of the provided run IDs and cluster IDs, and the gene IDs for all provided clusters.
#
# Note - this is quite slow so I don't suggest using it.

import fileinput, sqlite3, optparse
from FileLocator import *
from ClusterFuncs import *

# Get input arguments
usage = "%prog [options] < runid_clusterid_table > gene_id_list"
description = """Given a list of run IDs and cluster IDs, 
returns a list of all genes present in those run ID \ cluster ID pairs"""
parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-r", "--rcolumn", help="Column number (start from 1) for run ID", action="store", type="int", dest="runcolumn", default=1)
parser.add_option("-c", "--ccolumn", help="Column number (start from 1) for cluster ID", action="store", type="int", dest="clustercolumn", default=2)
(options, args) = parser.parse_args()

rc = options.runcolumn - 1 # Convert to Pythonic indexes
cc = options.clustercolumn - 1

con = sqlite3.connect(locateDatabase())
cur = con.cursor()

run_cluster_ids = set()
for line in fileinput.input("-"):
    spl = line.strip('\r\n').split("\t")
    runid = spl[rc]
    clusterid = spl[cc]
    run_cluster_ids.add( (runid, clusterid) )

for run_cluster_id in run_cluster_ids:
    geneids = getGenesInCluster(run_cluster_id[0], run_cluster_id[1], cur)
    for gene in geneids:
        print "%s\t%s\t%s" %(run_cluster_id[0], run_cluster_id[1], gene)

con.close()
