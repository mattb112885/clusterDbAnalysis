#!/usr/bin/env python

# Generates a table that looks like this:
#
#[]	    []	    [Org 1] [Org 2]    ...  [Org n]
#[RunID]    [ClID]  [PEGIDs][PEGIDs]   ...  [PEGIDs]
#
# from the "clusterorgs"
# Prints NONE if no PEG IDs exist in that cluster for that organism...

import sqlite3
import optparse
from FileLocator import *
from sanitizeString import *
from ClusterFuncs import *

# This is only for self-documenting purposes - this function takes no arguments.                                                                                                                              
usage="%prog > presence_absence_pegid_table"
description="""Generates a presence/absence table for every cluster in every run 
in the database and puts the peg IDs for any genes in the cluster in the appropriate 
row or NONE for absent genes. Takes no input arguments and exports the table to 
stdout. NOTE: Any organisms not included in a cluster run will be given NONE's for 
all clusters in that run - be aware of this! Also... if user-specified data are
included, we include them in the presence-absence table by joining the clusters."""
parser = optparse.OptionParser(usage=usage, description=description)
(options,args) = parser.parse_args()

con = sqlite3.connect(locateDatabase())
cur = con.cursor()

cur.execute("SELECT organism FROM organisms;")
orgList = set()
for rec in cur:
    orgList.add(rec[0])

# [runid, clusterid] --> {organism --> genelist}
rc2go = {}
cur.execute("SELECT * from clusterorgs;")
for rec in cur:
    rectup = (rec[0], rec[1])
    if rectup in rc2go:
        orgdict = rc2go[rectup]
        if rec[3] in orgdict:
            orgdict[rec[3]].append(rec[2])
        else:
            orgdict[rec[3]] = [ rec[2] ]
    else:
        orgdict = {}
        orgdict[rec[3]] = [ rec[2] ]
        rc2go[rectup] = orgdict

# Attempt to add user-added genes to this list
# We can only do that if the cluster and run are defined for the organism
# in question.
try:
    cur.execute("SELECT runid, clusterid, organismid, user_geneid FROM user_genes WHERE runid IS NOT NULL and clusterid IS NOT NULL;")
    for res in cur:
        runid = res[0]
        clusterid = res[1]
        # We have to do this to avoid overwriting the cursor we are looping over.
        cur2 = con.cursor()
        organismname = organismIdToName(res[2], cur2, issanitized=False)
        cur2.close()
        geneid = res[3]
        rectup = (runid, clusterid)
        if rectup not in rc2go:
            sys.stderr.write("WARNING: User-specified cluster-run pair (%s, %s) not found in the database\n" %(runid, clusterid))
            continue
        if organismname in rc2go[rectup]:
            rc2go[rectup][organismname].append(geneid)
        else:
            rc2go[rectup][organismname] = [ geneid ]
except sqlite3.OperationalError:
    sys.stderr.write("No user-specified genes are loaded in the database (if you have user-specified genes use setup_step5.sh to load them\n")
    pass
except:
    sys.stderr.write(str(sys.exc_info()))
    raise

# Generate title row
titleRow = "\t".join(orgList)
tileRow =  "%s\t%s\t%s\t%s" %("RunID", "ClusterID", "SampleAnnote", titleRow)

mytable = []
# Generate table
for rectup in rc2go:
    myList = []
    myrunid = rectup[0]
    myclusterid = rectup[1]
    myorgdict = rc2go[rectup]
    # Get most common annotation
    myannote = findRepresentativeAnnotation(myrunid, myclusterid, cur)
    # Get the organism peg list
    myList.append(myrunid)
    myList.append(myclusterid)
    if '\t' in myannote:
        sys.stderr.write("WARNING: A gene in run ID %s with cluster ID %s had an annotation with tabs in it. This probably indicates a problem with parsing one of the input files.\n" %(myrunid, myclusterid))
        myannote = myannote.replace("\t", "")
    myList.append(myannote)
    for org in orgList:
        if org in myorgdict:
            myList.append(";".join(myorgdict[org]))
        else:
            myList.append("NONE")
    mytable.append(tuple(myList))

# Generate SQL table with this info in it.
cur.execute("DROP TABLE IF EXISTS presenceabsence;")

cmd = """CREATE TABLE presenceabsence (
"runid" VARCHAR(128),
"clusterid" INT,
"annote" VARCHAR(2048)"""

for org in orgList:
    cmd += """, %s VARCHAR(128)""" %(sanitizeString(org, False))

cmd += ");"

cur.execute(cmd)

for ln in mytable:
    cmd = "INSERT INTO presenceabsence VALUES ("
    for l in ln:
        cmd += "?,"
    cmd = cmd.rstrip(",")
    cmd += ");"
    cur.execute(cmd, ln)
# We have to commit because we added a table.
con.commit()

cur = con.cursor()
# These help us more quickly get subsets of a presence absence table.
cur.execute("CREATE INDEX presenceabsenceclusters ON presenceabsence(clusterid);")
cur.execute("CREATE INDEX presenceabsenceruns ON presenceabsence(runid);")
con.close()
