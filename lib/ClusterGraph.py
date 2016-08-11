#!/usr/bin/env python

'''This library contains functions for generating a GML file for a particular
cluster of interest.

It includes helper files for making pretty colors. The intention is for you to import
the resulting GML file into Cytoscape or similar and use their rendering capabilities
to view it.

WARNING - Cytoscape 3 beta's GML reader is broken; it does not read attributes correctly (this has
been fixed and targeted for 3.0.1). If you want to read the attributes in go back to cytoscape
2.8.3

'''

import matplotlib
import networkx
import operator
import sys

from sanitizeString import *
from ClusterFuncs import *

def getHexFromScore( value, minvalue, maxvalue ):
    '''Given a score (value) and a minimum and maximum possible value for that score
    (minvalue and maxvalue), calculates a ratio between 0 and 1 and then uses MatPlotlib
    to convert this into a Hex color for export into GraphML (Note - depending on what
    Cytoscape can handle I might need to change this into something else)

    If value is less than minvalue it is treated as equal to minvalue, and if it is greater
    than maxvalue it is treated as equal to maxvalue for this purpose...

    Also it is assumed that higher values are better (so 0 gets white and is essentially invisible)'''

    if minvalue >= maxvalue:
        raise ValueError("minvalue must be less than maxvalue")

    if value < minvalue:
        ratio = 0
    elif value > maxvalue:
        ratio = 1
    else:
        ratio = ( float(value) - float(minvalue) ) / ( float(maxvalue) - float(minvalue) )

    # We want 0 to be white, not black.
    ratio = 1 - ratio

    # MATPLOTLIB requires the ratio to be put into the function as a string
    ratio = str(ratio)
    CC = matplotlib.colors.ColorConverter()
    rgb = CC.to_rgb(ratio)
    hex_value = matplotlib.colors.rgb2hex(rgb)
    return hex_value

def makeNetworkObjectFromBlastResults( blastres, score_method, cutoff, cur ):
    '''blastres is a list of lists of BLAST results in the same
    order as in the blastres_selfbit table from the database.

    cur is the sqlite object for the database (needed to get metadata)
    '''

    G = networkx.Graph(layoutAlgorithm="Prefuse Force Directed Layout", )

    getqueries = operator.itemgetter(0)
    gettargets = operator.itemgetter(1)

    querygenes = list(map(getqueries, blastres))
    targetgenes = list(map(gettargets, blastres))
    
    genelist = list(set(querygenes + targetgenes))
    for gene in genelist:
        geneinfo = getGeneInfo( [ gene ], cur )
        print(gene, geneinfo)
        # Not sure if the string sanitizing will be necessary.
        G.add_node(gene, organism=geneinfo[0][1], annotation=geneinfo[0][9])

    # Note - we don't care here if the score is symmetric or not, we just want to visualize it.
    scores = calculateScoreFromBlastres( blastres, score_method, cutoff, include_zeros = False, needsymmetric = False )
    minscore = 0
    maxscore = max( list(map(operator.itemgetter(2), scores )) )
    for score in scores:
        # Omit self-hits
        if score[0] == score[1]:
            continue
        hex_color = getHexFromScore(score[2], minscore, maxscore)
        G.add_edge(score[0], score[1], weight=score[2], metric=score_method, cutoff = cutoff, graphics = { "fill" : hex_color })
        
    return G

def getGraphForCluster( runid, clusterid, score_method, cutoff, cur, blastn=False ):
    '''Get a networkx graph object for a cluster (clusterid) in a specific cluster run (runid)
    given a scoring metric and a cutoff.

    The graph is made of all the edges with BLAST results above the cutoff (by the specified metric).
    
    cur is a sqlite cursor object. If blastn is specified as TRUE, we make the graph based on
    BLASTN results instead of BLASTP results (BLASTP is the default)
    '''

    geneid_list = getGenesInCluster(runid, clusterid, cur)
    blastres = getBlastResultsBetweenSpecificGenes(geneid_list, cur, blastn=blastn)
    G = makeNetworkObjectFromBlastResults(blastres, score_method, cutoff, cur)

    return G

def exportGraphToGML(G, filename):
    '''Export a networkx graph object G to GML with filename "filename"
    Note - I use GML mostly because GraphML doesn't work with the graphics attribute.'''    
    networkx.write_gml(G, filename)

'''
# Example usage
from FileLocator import *
import sqlite3
con = sqlite3.connect(locateDatabase())
cur = con.cursor()
G = getGraphForCluster("all_I_1.7_c_0.4_m_maxbit", 1524, "maxbit", 0.4, cur)
exportGraphToGML(G, "all_1524_maxbit_0.4_gml.gml")

con.close()
'''
