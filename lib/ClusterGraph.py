#!/usr/bin/python

'''This library contains functions for generating a GraphML file for a particular
cluster of interest.

It includes helper files for making pretty colors. The intention is for you tu import
the resulting GraphML file into Cytoscape or similar and use their rendering capabilities
to view it.

'''

import matplotlib
import networkx
import sys
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

def makeNetworkObjectFromBlastResults( blastres, score_method, cutoff ):
    '''blastres is a list of lists of BLAST results in the same
    order as in the blastres_selfbit table from the database'''
    G = networkx.Graph(layoutAlgorithm="Prefuse Force Directed Layout", )
    for res in blastres:
        # We don't want self-hits.
        if res[0] == res[1]:
            continue

        # Query gene (FIXME - should add some metadata in here - organism, location, annotation, etc...)
        G.add_node(res[0])
        # Target gene
        G.add_node(res[1])
        # Scoring method for now is implemented as maxbit. FIXME: Need to make this more generic and add
        # a function for it so I don't just keep re-implementing it differently in different functions
        # Also - note that I want the minscore to be the cutoff - anything below that gets white...
        minscore = cutoff
        maxscore = 1

        bitscore = float(res[11])
        query_selfbit = float(res[12])
        target_selfbit = float(res[13])
        score = bitscore/max(query_selfbit, target_selfbit)

        if score < cutoff:
            continue

        hex_color = getHexFromScore(score, minscore, maxscore)

        G.add_edge(res[0], res[1], weight=score, metric=score_method, graphics = { "fill" : hex_color })
        
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
    G = makeNetworkObjectFromBlastResults(blastres, score_method, cutoff)

    return G

def exportGraphToGML(G, filename):
    '''Export a networkx graph object G to GML with filename "filename"
    Note - I use GML mostly because GraphML doesn't work with the graphics attribute.'''    
    networkx.write_gml(G, filename)

# Example usage
'''
from FileLocator import *
import sqlite3
con = sqlite3.connect(locateDatabase())
cur = con.cursor()
G = getGraphForCluster("all_I_1.7_c_0.4_m_maxbit", 1524, "maxbit", 0.4, cur)
exportGraphToGML(G, "all_1524_maxbit_0.4_gml.gml")

con.close()
'''
