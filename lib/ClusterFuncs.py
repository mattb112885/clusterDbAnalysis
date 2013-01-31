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
    
