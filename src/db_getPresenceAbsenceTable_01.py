#!/usr/bin/python

# Generates a table that looks like this:
#
#[]	    []	    [Org 1] [Org 2]    ...  [Org n]
#[RunID]    [ClID]  [0/1]   [0/1]      ...  [0/1]
#
# from the "clusterorgs"
#
# Does not take any input arguments.

import sqlite3
import optparse
from locateDatabase import *

# This is only for self-documenting purposes - this function takes no arguments.
usage="%prog > presence_absence_01_table"
description="Generates a 0/1 presence/absence table for every cluster in every run in the database. This one is the one you should use as a starting point for a presence/absence fasta file if you want to make a tree out of that. Takes no input arguments and exports the table to stdout. NOTE: Any organisms not included in a cluster run will be given 0's for all clusters in that run - be aware of this!"
parser = optparse.OptionParser(usage=usage, description=description)
(options,args) = parser.parse_args()

con = sqlite3.connect(locateDatabase())
cur = con.cursor()

def getOneClusterOrganisms(runid, clusterid, cur):
    cur.execute("""SELECT organism, annotation FROM X
                   WHERE X.runid = ? AND X.clusterid = ?;""",
                  (runid, clusterid) )

    orgs = []
    annotes = []
    for l in cur:
        s = list(l)
        orgs.append(s[0])
        annotes.append(s[1])

    return set(orgs), max(set(annotes), key=annotes.count)

# Generate table list
cur.execute("SELECT organism FROM organisms;")
orgList = set()
for l in cur:
    s = list(l)
    orgList.add(s[0])
titleRow = "\t".join(orgList)
print "%s\t%s\t%s\t%s" %("RunID", "ClusterID", "SampleAnnote", titleRow)

# Iterate over clusters and Run IDs and get +/- for each...
cur.execute("SELECT DISTINCT runid, clusterid FROM clusters;")
execList = []
for l in cur:
    execList.append(list(l))

cur.execute("""CREATE TEMPORARY TABLE X AS
               SELECT clusters.*, processed.organism, processed.annotation  
               FROM processed                                                                                                                              
               INNER JOIN clusters ON clusters.geneid = processed.geneid;""")

cur.execute("CREATE INDEX tmpXRuns ON X(runid);")
cur.execute("CREATE INDEX tmpXclusters ON X(clusterid);")

for e in execList:
    oneRunOrglist, annote = getOneClusterOrganisms(e[0], e[1], cur)
    # Note- comma after the last argument means no newline
    print "%s\t%s\t%s" %(e[0], e[1], annote),
    
    for org in orgList:
        if org in oneRunOrglist:
            print "\t%d" %(1),
        else:
            print "\t%d" %(0),
    # Now we want a new line
    print ""
