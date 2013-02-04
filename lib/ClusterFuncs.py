#!/usr/bin/python

'''This module contains functions for more general cluster data manipulation
such as cross-referencing with annotations, organisms, etc...'''

import sys
from FileLocator import *

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

'''
# Test
from FileLocator import *
import sqlite3
con = sqlite3.connect(locateDatabase())
cur = con.cursor()
print getGeneInfo(["fig|339860.1.peg.123", "fig|339860.1.peg.124"], cur)
'''
