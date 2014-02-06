#!/usr/bin/env python

'''This module contains functions for more general cluster data manipulation
such as cross-referencing with annotations, organisms, BLAST results, etc...'''

import math
import sys
from FileLocator import *
from sanitizeString import *

def findRepresentativeAnnotation(runid, clusterid, cur):
    '''
    Identifies the most common annotation in a cluster/runID pair and returns a string
    containing that annotation.

    cur is a SQLite cursor
    '''

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
    '''
    Given a list of gene IDs, query the BLAST table to get a list of BLAST results
    containing the genes. The table is in -m9 format but with query and target self-bit scores
    added as the last two columns.

    blastn: TRUE if you want BLASTN results and FALSE if you want blastp results
    '''

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


def getValidBlastScoreMethods():
    '''
    List currently-supported BLAST scoring metrics
    '''
    symmetric_scores = ['maxbit', 'minbit', 'avgbit', 'normhsp']
    not_symmetric_scores = [ 'loge', 'selfbit', 'otherbit' ]
    return symmetric_scores, not_symmetric_scores

def calculateScoreFromBlastres(blastres, method, cutoff, include_zeros=False, needsymmetric = False):
    '''
    Standard function for computing BLAST scores from a list of BLAST results
    (from a library function or from parsing a BLAST results table from stdin).

    blastres is a list of lists where each internal list is the split representation
    of a BLAST resutls table (query, target, %ID, etc...)

    method must be symmetric if needsymmetric = True.

    cutoff is a floating point value. The score is treated as 0 if it is less than the cutoff.

    If include_zeros is FALSE then (query,target) pairs with a score of 0 (after applying the cutoff) are skipped. 
    Otherwise they are included with a score of 0

    Returns a list of tuples (query gene, target gene, score).
    '''

    symmetric, non_symmetric = getValidBlastScoreMethods()
    okMethods = symmetric + non_symmetric
    if method not in okMethods:
        raise ValueError("BLAST scoring method %s is not supported (valid methods: %s)" %(method, " ".join(okMethods)))
    if needsymmetric and method not in symmetric:
        raise ValueError("Calling function required a symmetric scoring method (query->target is the same as target->query) but method %s is not symmetric" %(method) )

    score_list = []
    for res in blastres:
        qgene = res[0]
        tgene = res[1]
        evalue = float(res[10])
        bitscore = float(res[11])
        qselfbit = float(res[12])
        tselfbit = float(res[13])
        hsplen = float(res[3])
        if method == "minbit":
            score = bitscore/min(qselfbit, tselfbit)
        elif method == "maxbit":
            score = float(bitscore)/max(qselfbit,tselfbit)
        elif method == "avgbit":
            score =  bitscore* 2 / (qselfbit+tselfbit)
        elif method == "normhsp":
            score = bitscore/hsplen
        elif method == "loge":
            score = -math.log10(evalue + 1E-200)
        elif method == "selfbit":
            score =  bitscore / (qselfbit)
        elif method == "otherbit":
            score =  bitscore / (tselfbit)
        if score < cutoff:
            score = 0
        if score == 0 and not include_zeros:
            continue
        score_list.append( ( qgene, tgene, score ) )

    return score_list

def getGeneNeighborhoods(geneid, clusterrunid, cur):
    '''
    Call the SQLITE database with cursor "cur" to obtain gene neighborhoods

    Returns a list of gene neighborhood arrays for gene ID "geneid" (as many as there are cached in
    the database - no cutoff for number of neighbors is applied)

    If no match is identified then we return an empty list.
    '''
    newgeneid = geneid
    cur.execute("""SELECT * from neighborhoods
                   WHERE neighborhoods.centergene=?;""", (newgeneid,))

    results = cur.fetchall()
    geneids = [l[1] for l in results]
    # Append the cluster ID to the table (needed for visualization purposes) - this will depend on the run ID used.
    # We want to do an IN query, but need to format w. correct number of ?s, so generate this string
    sql = "SELECT geneid, clusterid FROM clusters WHERE geneid IN ({seq}) AND runid = ?;".format(seq=','.join(['?']*len(geneids)))
    geneids.append(clusterrunid)
    cur.execute(sql, geneids)
    lookupcluster = dict(cur.fetchall())
    outdata = []
    for l in results:
        if l[1] in lookupcluster:
            outdata.append(l + (lookupcluster[l[1]], ) )
        else:
            sys.stderr.write("WARNING: Gene ID %s not found in cluster run %s\n" %(l[1], clusterrunid))
        pass
        
#    outdata = [l + (lookupcluster[l[1]],) for l in results]
    return outdata

def getGenesInRegion(contig_id, start, stop, cur, overhang=0):
    '''
    Call the SQLITE database with cursor "cur" to obtain genes in a particular region of DNA

    Returns a list of gene IDs for all genes between start and stop on the specified contig.
    The gene is allowed to "overhang" by at most "overhang" nucleotides (default = 0 , which means the
    entire gene must be present in the region)

    If start is less than stop we just reverse them (this is useful so that we can directly
    pass the numbers from TBLASTn results into this function...).
    '''

    interval_start = min(start, stop)
    interval_end = max(start, stop)

    # The MAX(processed.genestart, processed.geneend) >= start and MIN(processed.genestart, processed.geneend) <= stop
    # are needed to prevent us from getting genes that lie entirely outside the interval if the overhang is too large
    sql = """SELECT geneid FROM processed 
             WHERE processed.contig_mod = ?
             AND MIN(processed.genestart, processed.geneend) >= ?
             AND MAX(processed.genestart, processed.geneend) >= ?
             AND MIN(processed.genestart, processed.geneend) <= ?
             AND MAX(processed.genestart, processed.geneend) <= ?; """
    cur.execute(sql, (contig_id, interval_start - overhang, interval_start, interval_end, interval_end + overhang) )

    genelist = []
    for res in cur:
        genelist.append( res[0] )

    return genelist

def getGenesInCluster(runid, clusterid, cur):
    '''
    Get the genes in a cluster with ID clusterid from run with ID runid.
    cur is a SQLite cursor.  Returns a list of gene IDs
    '''

    q = """SELECT geneid FROM clusters
          WHERE runid = ? AND clusterid = ?"""
    cur.execute(q, ( runid, clusterid) )
    geneids = [ str(s[0]) for s in cur ]
    return geneids    
    
def getGeneInfo(genelist, cur):
    '''
    Given a list of gene IDs, returns the "geneinfo" as a list of tuples
    in the same format as expected from output of db_getGeneInformation.py and
    db_getClusterGeneInformation.py...
    '''

    q = "SELECT processed.* from processed WHERE processed.geneid = ?;"
    res = []
    for gene in genelist:
        cur.execute(q, (gene,))
        for k in cur:
            res.append( [ str(s) for s in k ] )
    return res

def getOrganismsInClusterRun(runid, cur):
    '''
    Get a list of organism names in a cluster run. Returns them as a list.
    '''
    organisms = []
    cur.execute("SELECT DISTINCT organism FROM clusterorgs WHERE runid=?", (runid, ))
    organisms = [ str(s[0]) for s in cur ]
    return organisms

def organismNameToId(orgname, cur, issanitized = False):
    '''
    Given an organism name, return the ID for that organism name.
    Use issanitized = True if the provided organism name has been sanitized
    with the sanitizeString() function
    '''
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
        raise ValueError("ERROR: Organism name %s not found in database" %(orgname))

def organismIdToName(orgid, cur, issanitized=False):
    '''
    Convert an organism ID (in format \d+\.\d+) into the name of the organism.

    If issanitized is True we expect the format \d+_\d+ instead.
    '''
    q = "SELECT organism, organismid FROM organisms;"
    cur.execute(q)
    idToOrg = {}
    for res in cur:
        if issanitized:
            idToOrg[sanitizeString(res[1], False)] = res[0]
        else:
            idToOrg[res[1]] = res[0]
    if orgid in idToOrg:
        return idToOrg[orgid]
    else:
        raise ValueError("ERROR: Organism ID %s not found in database")

def getContigIds(cur, orgid=None, orgname=None, issanitized=False):
    '''
    Obtain a list of contig IDs and return a list of them.

    By default, grabs ALL contigs.
    If orgid isn't None, grabs contigs only for organism with id "orgid".
    If the organism's name isn't none, grabs contigs only for the specified organism.
    '''

    if orgid is not None and orgname is not None:
        raise IOError("Cannot specify both an organism ID and an organism name")

    if orgname is not None:
        orgid = organismNameToId(orgname, cur, issanitized=issanitized)

    # Note - I didn't use the contigs table here because
    # the user might not have loaded up the contigs.
    q = "SELECT DISTINCT contig_mod FROM processed"

    if orgid is None and orgname is None:
        cur.execute(q)
    elif orgid is not None:
        cur.execute(q + " WHERE organismid = ?", (orgid,))

    contig_ids = []
    for res in cur:
        contig_ids.append(str(res[0]))
    return contig_ids

def getContigSequence(cur, contig_list):
    '''
    Accepts as input a list of contig IDs (in ITEP format).

    Returns dictionary from contig ID to sequence.
    '''
    q = "SELECT contig_mod, seq FROM contigs WHERE contig_mod = ?"
    contig_dict = {}
    for contig in contig_list:
        cur.execute(q, (contig, ))
        for res in cur:
            contig_dict[res[0]] = res[1]

    return contig_dict

'''
# Test the contig sequence ID and sequence getters
from FileLocator import *
import sqlite3
con = sqlite3.connect(locateDatabase())
cur = con.cursor()
contigids = getContigIds(cur, orgid="192952.1")
contigids_2 = getContigIds(cur, orgname="Methanosarcina acetivorans C2A")
print contigids_2
#print getContigSequence(cur, contigids)
'''

'''
# Test the organism name finder
from FileLocator import *
import sqlite3
con = sqlite3.connect(locateDatabase())
cur = con.cursor()
print organismIdToName("192952.1", cur)
print organismIdToName("192952_1", cur, issanitized=True)
print organismIdToName("7809798710412412.89713", cur)
'''
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

'''
from FileLocator import *
import sqlite3
con = sqlite3.connect(locateDatabase())
cur = con.cursor()
genes = getGenesInRegion("451756.88888.NZ_ABDX01000032.1", 30000, 40000, cur)
print genes
'''
