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
        organismname = organismIdToName(res[2], cur, issanitized=False)
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
    sys.stderr.write("No user-specified genes are loaded in the database (if you have user-specified genes use main5.sh to load them\n")
    pass

# Generate title row
titleRow = "\t".join(orgList)
tileRow =  "%s\t%s\t%s\t%s" %("RunID", "ClusterID", "SampleAnnote", titleRow)

mytable = []
# Generate table
for rectup in rc2go:
    myrunid = rectup[0]
    myclusterid = rectup[1]
    myorgdict = rc2go[rectup]
    # Get most common annotation
    myannote = findRepresentativeAnnotation(myrunid, myclusterid, cur)
    # Get the organism peg list
    myorgstr = ""
    for org in orgList:
        if org in myorgdict:
            myorgstr = "%s\t%s" %(myorgstr, ";".join(myorgdict[org]))
        else:
            myorgstr = "%s\t%s" %(myorgstr, "NONE")
    myorgstr = myorgstr.lstrip()
    myline = "%s\t%s\t%s\t%s" %(myrunid, myclusterid, myannote, myorgstr)
    mytable.append(myline)

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
    sp = ln.split("\t")
    cmd = "INSERT INTO presenceabsence VALUES ("
    for s in sp:
        cmd += "?,"
    cmd = cmd.rstrip(",")
    cmd += ");"
    cur.execute(cmd, tuple(sp))
# We have to commit because we added a table.
con.commit()
con.close()
    
