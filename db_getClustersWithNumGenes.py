#!/usr/bin/python

# This is a pipe command.
#
# Pipe in a run ID or list of run IDs
#
# For each specified run ID, returns a 2-column table containing the run ID in the first column
# and the cluster ID's with the specified number of genes (regardless of organisms) in the
# second column.

import fileinput, sqlite3, optparse, os

parser = optparse.OptionParser()
parser.add_option("-r", "--rcolumn", help="Column number (start from 1) for run ID", action="store", type="int", dest="runcolumn", default=1)
parser.add_option("-n", "--numcluster", help="Desired number of genes in each cluster to extract", action="store", type="int", dest="numcluster", default=None)
(options, args) = parser.parse_args()

rc = options.runcolumn - 1 # Convert to Pythonic indexes                                                                                                                                                      

if options.numcluster == None:
    print "ERROR in db_getClustersWithNumGenes: Must specify number of elements desired in each cluster!"
    os.system("./src/db_getClustersWithNumGenes.py -h");
    raise IOError
else:
    nc = options.numcluster

con = sqlite3.connect("db/methanosarcina")
cur = con.cursor()

for line in fileinput.input("-"):
    spl = line.strip().split("\t")
    cur.execute("""SELECT runid, clusterid FROM clusters WHERE clusters.runid = ?
                   GROUP BY clusterid HAVING count(*) = ?;""", (spl[rc], int(nc)))

    for l in cur:
        s = list(l)
        stri = "\t".join(str(t) for t in s)
        print stri

