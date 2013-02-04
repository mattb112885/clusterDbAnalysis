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


import sqlite3
from FileLocator import *

con = sqlite3.connect(locateDatabase())
cur = con.cursor()

all_geneids = []
for line in open("acetivorans_all_genes", "r"):
    all_geneids.append(line.strip("\r\n"))

core_geneids = []
for line in open("acetivorans_core_genes", "r"):
    core_geneids.append(line.strip("\r\n"))

all_seqFeatures = geneListToSeqFeatureList(all_geneids, cur)
core_seqFeatures = geneListToSeqFeatureList(core_geneids, cur)

nameToFeatureList = {}
nameToColor = {}
nameToFeatureList["all"] = all_seqFeatures
nameToColor["all"] = colors.blue
nameToFeatureList["core"] = core_seqFeatures
nameToColor["core"] = colors.red

myDiagram = makeGenomeDiagram(nameToFeatureList, "core_vs_all", nameToColor)
myDiagram.draw(format="circular", circular=True, orientation="landscape", pagesize='A4', start=0, end=5751491)
myDiagram.write("example.svg", "SVG")

con.close()
