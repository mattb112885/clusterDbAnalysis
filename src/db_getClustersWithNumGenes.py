#!/usr/bin/env python

# This is a pipe command.
#
# Pipe in a run ID or list of run IDs
#
# For each specified run ID, returns a 2-column table containing the run ID in the first column
# and the cluster ID's with the specified number of genes (regardless of organisms) in the
# second column.

import fileinput, sqlite3, optparse, sys
from FileLocator import *

header = [ "run_id", "cluster_id" ]
usage = """%prog -n numgenes [options] < run_ids > clusters_with_specified_num_genes

Output: """ + " ".join(header)
description = "Get all of the clusters with the specified number of genes in the specified cluster runs"
parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-r", "--rcolumn", help="Column number (start from 1) for run ID", action="store", type="int", dest="runcolumn", default=1)
parser.add_option("-n", "--numcluster", help="Desired number of genes in each cluster to extract", action="store", type="int", dest="numcluster", default=None)
parser.add_option("--header", help="Specify to add header to the output file (useful if you want to take the results and put into Excel or similar programs)",
                  action="store_true", default=False)
(options, args) = parser.parse_args()

rc = options.runcolumn - 1 # Convert to Pythonic indexes                                                                                                                                                      

if options.numcluster == None:
    sys.stderr.write("ERROR in db_getClustersWithNumGenes: Must specify number of elements desired in each cluster!\n")
    raise IOError
else:
    nc = options.numcluster

con = sqlite3.connect(locateDatabase())
cur = con.cursor()

if options.header:
    print "\t".join(header)
for line in fileinput.input("-"):
    spl = line.strip('\r\n').split("\t")
    cur.execute("""SELECT runid, clusterid FROM clusters WHERE clusters.runid = ?
                   GROUP BY clusterid HAVING count(*) = ?;""", (spl[rc], int(nc)))
    for l in cur:
        s = list(l)
        stri = "\t".join(str(t) for t in s)
        print stri

con.close()
