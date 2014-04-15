#!/usr/bin/env python

'''
This library file contains functions for generating and manipulating Biopython
graphics objects.
'''

import colorsys
import itertools
import math
import numpy
import os
import tempfile

from Bio.Graphics import GenomeDiagram
from Bio.SeqFeature import SeqFeature, FeatureLocation
from reportlab.lib import colors as rcolors

from sanitizeString import *
from ClusterFuncs import *
from FileLocator import *
from TreeFuncs import splitTblastn

###########################################
# Functions for making SeqFeature objects #
###########################################

def makeSeqFeature(geneid, cur):
    '''
    Make a BioPython SeqFeature object for a gene with ITEP ID geneid
    '''

    geneinfo = getGeneInfo( [ geneid ], cur )
    geneinfo = geneinfo[0]
    start = int(geneinfo[5])
    stop = int(geneinfo[6])
    strand = int(geneinfo[8])
    feature = SeqFeature(FeatureLocation(start, stop), strand=strand, id=geneid)
    # This can be overwritten by other functions but we need a placeholder.
    feature.qualifiers["cluster_id"] = -1
    return feature

def makeSeqFeaturesForGeneNeighbors(genename, runid, cur):
    '''                                                                                                              
    Create seqFeature objects for a gene and its neighbors.

    genename is the ITEP ID for a gene.
    runid is a run ID

    The function returns a list of BioPython SeqFeature objects for the specified gene 
    and its neighbors.
    
    If the gene is not found it returns an empty list.
    '''
    outdata = getGeneNeighborhoods(genename, runid, cur)
    seqfeatures = []
    for neargene in outdata:
        feature = makeSeqFeature(neargene[1], cur)
        feature.qualifiers["cluster_id"] = int(neargene[8])
        seqfeatures.append(feature)
    return seqfeatures

def makeSeqObjectsForTblastnNeighbors(tblastn_id, clusterrunid, cur, N=200000):
    '''
    Given a tBLASTn ID and a dictionary from sanitized contig IDs (which is what will be 
    present in the TBLASTN id) to non-sanitized IDs (which are what is in the database),
    returns a list of seq objects INCLUDING the TBLASTN hit itself so that we can show that
    on the region drawing.

    We pick an N large enough to get at least one gene and then pick the closest one and get
    all of its neighbors with a call to makeSeqFeaturesForGeneNeighbors() and just tack the TBLASTN
    onto it. 
    '''
    # Lets first get the contig and start/stop locations (which tell us teh strand) out of    
    # the TBLASTN id. This returns a ValueError if it fails which the calling function can catch if needed. 
    sanitizedToNot = getSanitizedContigList(cur)

    # The tBLASTn ID holds information on where the hit was located.
    contig,start,stop = splitTblastn(tblastn_id)
    if contig in sanitizedToNot:
        contig = sanitizedToNot[contig]

    # Create a seq object for the TBLASTN hit itself
    start = int(start)
    stop = int(stop)
    if start < stop:
        strand = +1
    else:
        strand = -1
    tblastn_feature = SeqFeature(FeatureLocation(start, stop), strand=strand, id=tblastn_id)
    tblastn_feature.qualifiers["cluster_id"] = -1

    # Find the neighboring genes.
    neighboring_genes = getGenesInRegion(contig, start-N, stop+N, cur)
    if len(neighboring_genes) == 0:
        sys.stderr.write("WARNING: No neighboring genes found for TBLASTN hit %s within %d nucleotides in contig %s\n" %(tblastn_id, N, contig))
        return [ tblastn_feature ]
    else:
        neighboring_geneinfo = getGeneInfo(neighboring_genes, cur)

    # Find the closest gene to ours and get the clusters for those neighbors based on the specific clusterrunid                                                                                               
    minlen = N
    mingene = None
    minstrand = None
    for geneinfo in neighboring_geneinfo:
        genestart = int(geneinfo[5])
        geneend = int(geneinfo[6])
        distance = min( abs(genestart - start), abs(geneend - start), abs(genestart - stop), abs(geneend - stop))
        if distance < minlen:
            mingene = geneinfo[0]
            minlen = distance

    neighboring_features = makeSeqFeaturesForGeneNeighbors(mingene, clusterrunid, cur)

    # Add the TBLASTN itself and return it.
    neighboring_features.append(tblastn_feature)
    return neighboring_features

###############################
# Drawing functions           #
###############################

def regionlength(seqfeatures):
    ''' Find the beginning and end of nucleotides spanning a set of gene locations '''
    location = [(int(loc.location.start), int(loc.location.end)) for loc in seqfeatures]
    starts, ends = zip(*location)
    #have to compare both, as some are reversed
    start = max(max(starts),max(ends))
    end = min(min(starts),min(ends))
    return start, end

def make_region_drawing(seqfeatures, getcolor, centergenename, maxwidth, tempdir=None, imgfileloc = None, label=False, labeltype = 'clusterid' ):
    '''
    Makes a PNG figure for regions with a given color mapping, set of gene locations... 

    seqfeatures is a list of SeqFeature objects (with the cluster_id qualifier)
    getcolor is a map from cluster ID to the desired color
    centergenename is the ID (as in the seqFeature) for the gene you wish to have in the middle.
    maxwidth is the maximum width of the image (in pixels)

    if label is TRUE we add a label to each of the arrows.

    labeltype: 'clusterid' : Add numeric cluster ID to each feature
               'aliases'   : Add a underscore-delimited list of aliases to each feature (aliases file is located in $ITEP_ROOT/aliases/aliases)

    If an output file is unspecified the region drawing will be made temporary in tempdir (or if that isn't specified either, it will be made temporary
    in a new directory created in /tmp/ or the default temporary directory for Python)
    '''

    # The files are not automatically deleted
    # but at least this prevents collisions.
    # A user who wants to clean up should specify a temporary directory and delete it afterward.
    if imgfileloc is None:
        if tempdir is None:
            tempdir = tempfile.gettempdir()
        imghandle = tempfile.NamedTemporaryFile(delete=False, dir=tempdir)
        imgfileloc = imghandle.name

    # Set up an entry genome diagram object                                                                                                                                                                   
    gd_diagram = GenomeDiagram.Diagram("Genome Region")
    gd_track_for_features = gd_diagram.new_track(1, name="Annotated Features")
    gd_feature_set = gd_track_for_features.new_set()

    # Some basic properties of the figure itself
    arrowshaft_height = 0.3
    arrowhead_length = 0.3
    scale = 20     #AA per px for the diagram

    if label:
        # y-margins need to be bigger if we have labels
        # Since labels could be on the top or the bottom we need to account for either possibility
        # Labels are longer for other types than for cluster IDs (which are typically < 6 characters)
        if labeltype == 'clusterid':
            yt = 0.1
            yb = 0.1
            ht = 350
            default_fontsize = 25
        else:
            yt = 0.3
            yb = 0.3
            ht = 725
            default_fontsize = 16 # Font size for genome diagram labels
    else:
        yt = 0
        yb = 0
        ht = 250
        default_fontsize = 16 # Doesnt matter

    geneIdToAlias = {}
    if label:
        if labeltype == "aliases":
            if os.path.exists(locateAliasesFile()):
                for line in open(locateAliasesFile(), "r"):
                    spl = line.strip("\r\n").split("\t")
                    if spl[0] in geneIdToAlias:
                        geneIdToAlias[spl[0]].append(spl[1])
                    else:
                        geneIdToAlias[spl[0]] = [ spl[1] ]

    # Build arrow objects for all of our features.
    for feature in seqfeatures:
        bordercol=rcolors.white

        if feature.id == centergenename:
            bordercol=rcolors.red
            centerdstart, centerend = int(feature.location.start), int(feature.location.end)
            centerdstrand = feature.strand
        color = getcolor[feature.qualifiers["cluster_id"]]

        name = feature.id
        if label:
            if labeltype == "aliases":
                if feature.id in geneIdToAlias:
                    name = "_".join(geneIdToAlias[feature.id])                   
                    if len(name) > 30:
                        name = name[0:30]
            elif labeltype == "clusterid":
                name = str(feature.qualifiers["cluster_id"])
            else:
                raise IOError("Invalid labeltype")

        # 90 degrees to avoid the rightmost feature going off the screen
        gd_feature_set.add_feature(feature, name = name,
                                   color=color, border = bordercol,
                                   sigil="ARROW", arrowshaft_height=arrowshaft_height, arrowhead_length = arrowhead_length,
                                   label=label,  label_angle=90, label_size = default_fontsize, label_position = 'middle'
                                   )

    start, end = regionlength(seqfeatures)
    pagew_px = maxwidth / scale
    #offset so start of gene of interest lines up in all the figures
    midcentergene = abs(centerend - centerdstart)/2 + min(centerdstart, centerend)
    l2mid = abs(midcentergene - start)
    r2mid = abs(midcentergene - end)
    roffset = float((pagew_px/2) - (l2mid/scale))
    loffset = float((pagew_px/2) - (r2mid/scale))

    gd_diagram.draw(format="linear", start=start, end=end, fragments=1, pagesize=(ht, pagew_px), xl=(loffset/pagew_px), xr=(roffset/pagew_px), yt=yt, yb=yb )

    gd_diagram.write(imgfileloc, "PNG")

    #flip for reversed genes
    if centerdstrand == -1:
        os.system("convert -rotate 180 %s %s" % (imgfileloc, imgfileloc))
    return imgfileloc


def makeClusterColorMap(seqfeatures, greyout):
    '''
    seqfeatures is a list of SeqFeature objects.

    Generate a color map based on presence and absence of genes in clusters.

    If the number of genes in all the provided seqfeatures is less than the cutoff,
    the color is set to gray. If you have seqfeatures for multiple genes you want to test at once,
    concatenate them before calling this function.
    '''
    allclusters = []
    for feature in seqfeatures:
        allclusters.append(feature.qualifiers["cluster_id"])

    uniqueclusters = set(allclusters)

    # Get clusters that have enough members to bother trying to color them (as determined by
    # the greyout keyword)
    multipleclusters = [c for c in uniqueclusters if allclusters.count(c) >= greyout]

    # Don't die if nothing has enough clusters...
    if len(multipleclusters) > 0:
        getcolor = colormap(multipleclusters)
    else:
        getcolor = {}

    #also add in grey (0.5,0.5,0.5 in RGB) for all others                                                                                                                                                     
    singleclusters = [c for c in uniqueclusters if allclusters.count(c) < greyout]
    getcolor.update([(sc, (0.5,0.5,0.5)) for sc in singleclusters])
    return getcolor

##############################
# Putting it all together... #
##############################
def makeSingleGeneNeighborhoodDiagram(geneid, runid, cur, tempdir=None, imgfileloc = None, labeltype = 'clusterid'):
    '''
    Make a genome context diagram for a single gene with ITEP ID geneid.
    '''
    seqfeatures = makeSeqFeaturesForGeneNeighbors(geneid, runid, cur)
    getcolor = makeClusterColorMap(seqfeatures, 1)
    start, end = regionlength(seqfeatures)
    width = abs(end - start)
    imgfileloc = make_region_drawing(seqfeatures, getcolor, geneid, width, tempdir=tempdir, imgfileloc = imgfileloc, label=True, labeltype = labeltype)
    return imgfileloc

##########################
# Other utilities        #
##########################

def RGB_to_hex(RGBlist):
    '''
    Convert an RGB color into a HEX string (required for some display functions)
    '''
    n = lambda x: int(x*255)
    RGB256 = [(n(r),n(g),n(b)) for r,g,b in RGBlist]
    colors = ['#%02x%02x%02x' % (r, g, b) for r, g, b in RGB256]
    return colors

def colormap(valuelist):
    '''
    Generate a list of divergent colors for use with labeling SeqFeature objects
    '''
    values = numpy.unique(valuelist)
    N = len(values)
    #we will vary in 2 dimensions, so this is how many steps in each
    perm = int(math.ceil(math.sqrt(N)))
    #need offset, as humans can't tell colors that are unsaturated apart
    H = [(x*1.0/perm) for x in range(perm)]
    S = [(x*1.0/perm)+0.2 for x in range(perm)]
    #we will use this to truncate at correct length
    V = [0.7]*N
    # Create all combinations of our colors.                                                                                                                                                           
    HS = itertools.product(H, S)
    H, S = zip(*HS)
    HSV = zip(H,S,V)
    RGB = [colorsys.hsv_to_rgb(h,s,v) for h, s, v in HSV]
    colorlookup = dict(zip(values, RGB[:N]))
    return colorlookup
