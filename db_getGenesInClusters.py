#!/usr/bin/python

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

# Get input arguments
parser = optparse.OptionParser()
parser.add_option("-r", "--rcolumn", help="Column number (start from 1) for run ID", action="store", type="int", dest="runcolumn", default=1)
parser.add_option("-c", "--ccolumn", help="Column number (start from 1) for cluster ID", action="store", type="int", dest="clustercolumn", default=2)

(options, args) = parser.parse_args()
rc = options.runcolumn - 1 # Convert to Pythonic indexes
cc = options.clustercolumn - 1

con = sqlite3.connect("db/methanosarcina")
cur = con.cursor()

# Make a temporary table with all of the clusters that we want to actually extract

query = "CREATE TEMPORARY TABLE allclusters( runid VARCHAR(256), clusterid INT );"
cur.execute(query);
for line in fileinput.input("-"):
    spl = line.strip('\r\n').split("\t")
    query = "INSERT INTO allclusters VALUES (?, ?);"
    con.execute(query, (spl[rc], spl[cc]))
    
# Remove duplication... especially important since we split this thing up into multiple queries now.
cur.execute("CREATE TEMPORARY TABLE desiredclusters AS SELECT DISTINCT * FROM allclusters;")

# Get table of gene IDs
cur.execute("""SELECT clusters.* FROM clusters
               INNER JOIN desiredclusters ON desiredclusters.runid = clusters.runid
                      AND desiredclusters.clusterid = clusters.clusterid;""")

for l in cur:
    s = list(l)
    stri = "\t".join(str(t) for t in s)
    print stri

con.close()
