#!/usr/bin/python

'''This module contains functions for more general cluster data manipulation
such as cross-referencing with annotations, organisms, BLAST results, etc...'''

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

def getBlastResultsBetweenSpecificGenes(geneids, cur, blastn=False):
    '''Given a list of gene IDs, query the BLAST table to get a list of BLAST results
    containing the genes. The table is in -m9 format but with query and target self-bit scores
    added as the last two columns.

    blastn: TRUE if you want BLASTN results and FALSE if you want blastp results'''

    # FIXME - Can I get equivalent performance by passing in lots of queries at once instead of making a temporary table?
    # Expunging the temporary tables would allow us not to have to give "w" to anyone that wants to use the database.
    cur.execute("""CREATE TEMPORARY TABLE desiredgenes ("geneid" VARCHAR(128), FOREIGN KEY(geneid) REFERENCES rawdata(geneid));""")
    for geneid in geneids:
        cur.execute("INSERT INTO desiredgenes VALUES (?);", (geneid, ) )

    # Generate a list of blast results with query matching one of the desiredgenes
    if blastn:
        tbl = "blastnres_selfbit"
    else:
        tbl = "blastres_selfbit"   

    cmd = """SELECT %s.* FROM %s
         WHERE %s.targetgene IN (select geneid from desiredgenes)
         AND %s.querygene IN (select geneid from desiredgenes);""" %(tbl, tbl, tbl, tbl)
    cur.execute(cmd)

    resultTable = []
    for k in cur:
        resultTable.append( [ str(s) for s in k ] )

    # Since we're preserving cur we should clean up here.
    cur.execute("DROP TABLE desiredgenes;")
    return resultTable

def getGenesInCluster(runid, clusterid, cur):
    '''Get the genes in a cluster with ID clusterid from run with ID runid.
    cur is a SQLite cursor.  Returns a list of gene IDs'''

    q = """SELECT geneid FROM clusters
          WHERE runid = ? AND clusterid = ?"""
    cur.execute(q, ( runid, clusterid) )
    geneids = [ str(s[0]) for s in cur ]
    return geneids    
    
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

'''
# Test the cluster to BLAST results function
from FileLocator import *
import sqlite3
con = sqlite3.connect(locateDatabase())
cur = con.cursor()
genelist =  getGenesInCluster("all_I_1.7_c_0.4_m_maxbit", 1524, cur)
blastres = getBlastResultsBetweenSpecificGenes(genelist, cur, False)
for res in blastres:
    print "\t".join(res)
'''
