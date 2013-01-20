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

# This is only for self-documenting purposes - this function takes no arguments.                                                                                                                              
usage="%prog > presence_absence_pegid_table"
description="Generates a presence/absence table for every cluster in every run in the database and puts the peg IDs for any genes in the cluster in the appropriate row or NONE for absent genes. Takes no input arguments and exports the table to stdout. NOTE: Any organisms not included in a cluster run will be given NONE's for all clusters in that run - be aware of this!"
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
rc2genelist = {}
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
    if rectup in rc2genelist:
        rc2genelist[rectup].append(rec[2])
    else:
        rc2genelist[rectup] = [ rec[2] ]

# gene --> annotation
cur.execute("SELECT geneid, annotation FROM processed;")
gene2annote = {}
for rec in cur:
    gene2annote[rec[0]] = rec[1]

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
    annotelist = []
    for gene in rc2genelist[rectup]:
        annotelist.append(gene2annote[gene])
    myannote = max(set(annotelist), key=annotelist.count)
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

# Generate MySQL table with this info in it.
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
# Needed since we added a table...
con.commit()
con.close()
    
