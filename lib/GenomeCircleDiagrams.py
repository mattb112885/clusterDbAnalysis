#!/usr/bin/env python

from __future__ import print_function
import Bio.SeqFeature as SeqFeature
import Bio.Graphics.GenomeDiagram as GenomeDiagram
from reportlab.lib import colors
from ClusterFuncs import *
import operator

'''Contains scripts for generating circular genome objects from
the biopython GenomeDiagram library

Some of this is derived from the documentation
at http://www.bio-cloud.info/Biopython/en/ch15.html
'''

def createSeqFeature(start, stop, strand, idnum, featuretype="CDS"):
    '''Create a SeqFeature object (which needs to be added to a track
    to generate an image)'''
    location = SeqFeature.FeatureLocation(start, stop, strand)
    feature = SeqFeature.SeqFeature(location, type=featuretype, id=idnum)
    return feature


def makeGenomeDiagram(seqFeatureNameToList, genomeDiagramName, seqFeatureToColor = None):
    '''INPUT: seqFeatureNameToList: A dictionary
    ring_name --> [ seq_features_within_ring ]

    genomeDiagramName: A string defining the name of your object.

    optional: seqFeatureToColor: A dictionary
    ring_name --> color (a reportlab.lib color object)

    OUTPUT: A GenomeDiagram object storing all of those seqFeature objects
    within the appropriate rings.

    The order of rings will be in alphabetical order by name by default.
    Switch the order with gd_diagram.move_track(old_loc,new_loc)
    '''
    
    gd_diagram = GenomeDiagram.Diagram(genomeDiagramName)
    tracknum = 0
    for name in sorted(seqFeatureNameToList.keys()):
        tracknum += 1
        gd_track_for_features = gd_diagram.new_track(tracknum, name=name)
        gd_feature_set = gd_track_for_features.new_set()
        # Add the things we want to actually put ON the circle.
        for feature in seqFeatureNameToList[name]:
            if seqFeatureToColor is None:
                color = colors.blue
            else:
                color = seqFeatureToColor[name]
                gd_feature_set.add_feature(feature, color=color, label=False, sigil="ARROW", arrowshaft_height=0.5, name=feature.id )
#            gd_feature_set.add_feature(feature, color=color, label=True, label_size = 12, label_position = "middle", sigil="ARROW", arrowshaft_height=0.5, name=test_feature.id )
    return gd_diagram

def geneListToSeqFeatureList(geneid_list, cur):
    '''From a list of gene IDs, come up with a list of seqFeatures.
    Throws an error if not all of the genes are in the same organism...

    cur is a Sqlite3 cursor object'''

    # Get "geneinfo" which contains location and strand for all of those genes.
    geneinfo = getGeneInfo(geneid_list, cur)
    orglist = set(map(operator.itemgetter(1), geneinfo))
    if len(orglist) > 1:
        raise ValueError("ERROR: Specified gene ID list contains multiple organisms and doesnt make sense to add to a single genome diagram! Organisms: %s" %("\t".join(orglist)))

    seqFeatureList = []
    for info in geneinfo:
        myid = info[0]
        startloc = int(info[5])
        stoploc = int(info[6])
        # Need +/- 1
        strand = int(info[8])
        seq = createSeqFeature(startloc, stoploc, strand, myid)
        seqFeatureList.append(seq)

    return seqFeatureList

# Example usage
'''
import sqlite3
from FileLocator import *

con = sqlite3.connect(locateDatabase())
cur = con.cursor()

mazei_geneids = [ line.strip("\r\n") for line in open("mazei_core_genes") ]
methanosarcina_geneids = [ line.strip("\r\n") for line in open("methanosarcina_core_genes") ]
methanosarcinales_geneids = [ line.strip("\r\n") for line in open("methanosarcinales_core_genes") ]
methanosarcinales_methanocella_geneids = [ line.strip("\r\n") for line in open("methanosarcinales_methanocella_genelist") ]

mazei_seqFeatures = geneListToSeqFeatureList(mazei_geneids, cur)
methanosarcina_seqFeatures = geneListToSeqFeatureList(methanosarcina_geneids, cur)
methanosarcinales_seqFeatures = geneListToSeqFeatureList(methanosarcinales_geneids, cur)
methanosarcinales_methanocella_seqFeatures = geneListToSeqFeatureList(methanosarcinales_methanocella_geneids, cur)

nameToFeatureList = {}
nameToColor = {}
#nameToFeatureList["mazei"] = mazei_seqFeatures
#nameToColor["mazei"] = colors.blue
#nameToFeatureList["methanosarcina"] = methanosarcina_seqFeatures
#nameToColor["methanosarcina"] = colors.red
#nameToFeatureList["methanosarcinales"] = methanosarcinales_seqFeatures
#nameToColor["methanosarcinales"] = colors.green
nameToFeatureList["methanosarcinales_methanocella"] = methanosarcinales_methanocella_seqFeatures
nameToColor["methanosarcinales_methanocella"] = colors.blue

myDiagram = makeGenomeDiagram(nameToFeatureList, "mazei_conservation", nameToColor)
#myDiagram.move_track(3,0)
#myDiagram.move_track(1,3)
#myDiagram.move_track(0,1)
myDiagram.move_track(1, 2)
myDiagram.draw(format="circular", circular=True, orientation="landscape", pagesize='A4', start=0, end=4096345)
myDiagram.write("methanosarcinales_methanocella.svg", "SVG")

con.close()
'''
