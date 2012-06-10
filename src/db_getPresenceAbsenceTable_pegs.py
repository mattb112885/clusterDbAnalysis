#!/usr/bin/python

# Generates a table that looks like this:
#
#[]	    []	    [Org 1] [Org 2]    ...  [Org n]
#[RunID]    [ClID]  [PEGIDs][PEGIDs]   ...  [PEGIDs]
#
# from the "clusterorgs"
# Prints NONE if no PEG IDs exist in that cluster for that organism...

import sqlite3
import optparse

# This is only for self-documenting purposes - this function takes no arguments.                                                                                                                              
usage="%prog > presence_absence_pegid_table"
description="Generates a presence/absence table for every cluster in every run in the database and puts the peg IDs for any genes in the cluster in the appropriate row or NONE for absent genes. Takes no input arguments and exports the table to stdout. NOTE: Any organisms not included in a cluster run will be given NONE's for all clusters in that run - be aware of this!"
parser = optparse.OptionParser(usage=usage, description=description)
(options,args) = parser.parse_args()


con = sqlite3.connect("db/methanosarcina")
cur = con.cursor()

def getOneClusterOrganismsPegs(runid, clusterid, cur):
    cur.execute("""SELECT organism, annotation, geneid FROM X
                   WHERE X.runid = ? AND X.clusterid = ?;""",
                (runid, clusterid) )
    orgs = []
    annotes = []
    genes = []

    for l in cur:
        s = list(l)
        orgs.append(s[0])
        annotes.append(s[1])
        genes.append(s[2])

    # Build dictionary from organisms to a list of genes
    orgdict = {}
    for ii in range(len(orgs)):
        if orgs[ii] in orgdict:
            orgdict[orgs[ii]].append(genes[ii])
        else:
            orgdict[orgs[ii]] = [genes[ii]]

    return orgdict, max(set(annotes), key=annotes.count)


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
    orgdict, annote = getOneClusterOrganismsPegs(e[0], e[1], cur)
    # Note- comma after the last argument means no newline
    print "%s\t%s\t%s" %(e[0], e[1], annote),
    
    for org in orgList:
        if org in orgdict:
            print "\t%s" %(";".join(orgdict[org])),
        else:
            print "\t%s" %("NONE"),
    # Now we want a new line
    print ""
