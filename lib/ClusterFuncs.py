#!/usr/bin/python

'''This module contains functions for more general cluster data manipulation
such as cross-referencing with annotations, organisms, etc...'''

import sys
from FileLocator import *
from sanitizeString import *

def findRepresentativeAnnotation(runid, clusterid, cur):
    '''Identifies the most common annotation in a cluster/runID pair.
    cur a SQLite cursor'''

    q = """SELECT annotation FROM clusters
           INNER JOIN processed ON processed.geneid = clusters.geneid
           WHERE runid = ? AND clusterid = ?"""
    cur.execute(q, (runid, clusterid) )
    annotes = [ str(s[0]).lower() for s in cur ]
    maxcount = 0
    bestannote = None
    for annote in set(annotes):
        if annotes.count(annote) > maxcount:
            maxcount = annotes.count(annote)
            bestannote = annote
    return bestannote
    
def getGeneInfo(genelist, cur):
    '''Given a list of gene IDs, returns the "geneinfo" as a list of tuples
    in the same format as expected from output of db_getGeneInformation.py and
    db_getClusterGeneInformation.py...'''

    q = "SELECT processed.* from processed WHERE processed.geneid = ?;"
    res = []
    for gene in genelist:
        cur.execute(q, (gene,))
        for k in cur:
            res.append( [ str(s) for s in k ] )
    return res

def organismNameToId(orgname, cur, issanitized = False):
    '''Given an organism name, return the ID for that organism name.
    Use issanitized = True if the provided organism name has been sanitized
    with the sanitizeString() function'''
    q = "SELECT organism, organismid FROM organisms;"
    cur.execute(q)
    orgToId = {}
    for res in cur:
        if issanitized:
            orgToId[sanitizeString(res[0], False)] = res[1]
        else:
            orgToId[res[0]] = res[1]
    if orgname in orgToId:
        return orgToId[orgname]
    else:
        raise ValueError("ERROR: Organism name %s not found in database")

'''
# Test the organism ID finder
from FileLocator import *
import sqlite3
con = sqlite3.connect(locateDatabase())
cur = con.cursor()
print organismNameToId("Methanosarcina acetivorans C2A", cur)
print organismNameToId("Methanosarcina_acetivorans_C2A", cur, issanitized = True)
'''

'''
# Test
from FileLocator import *
import sqlite3
con = sqlite3.connect(locateDatabase())
cur = con.cursor()
print getGeneInfo(["fig|339860.1.peg.123", "fig|339860.1.peg.124"], cur)
'''
