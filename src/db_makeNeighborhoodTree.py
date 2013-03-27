#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This script plots regions around proteins using a newick tree of the proteins and the results of iTEP's db_getGeneNeighborhoods.py

Created on Wed Oct 17 12:40:42 2012

@author: jamesrh

NOTE: If you get an error complaining about missing fonts:

reportlab.graphics.renderPM.RenderPMError: Can't setFont(Times-Roman) missing the T1 files?

it means your installation of reportlab is missing the font files for Times New Roman
and you should install them. How to install them:

1: Download the .pfb files from the reportlab website. The current location is:
http://www.reportlab.com/ftp/fonts/pfbfer.zip

2: Extract them to one of the places reportlab looks for the files.
One of those locations is ${HOME}/fonts so follow these directions:

$ cd ~
$ mkdir fonts
$ cd fonts
$ wget http://www.reportlab.com/ftp/fonts/pfbfer.zip
$ unzip pfbfer.zip

When you run this again, it should work.

See the following for more details.

http://blog.gmane.org/gmane.comp.python.reportlab.user/month=20060701/page=5
"""

import colorsys
import math
import itertools
import optparse
import os
import sqlite3
import sys

from FileLocator import *
from TreeFuncs import *
from ClusterFuncs import *
from sanitizeString import *

import numpy
from reportlab.lib import colors as rcolors
from Bio.Graphics import GenomeDiagram
from Bio.SeqFeature import SeqFeature, FeatureLocation
from ete2 import Tree, faces, TreeStyle, TextFace, PhyloTree
from ete2 import Phyloxml, phyloxml

#from tempfile import NamedTemporaryFile
#from Bio import SeqIO

def RGB_to_hex(RGBlist):
    n = lambda x: int(x*255)
    RGB256 = [(n(r),n(g),n(b)) for r,g,b in RGBlist]
    colors = ['#%02x%02x%02x' % (r, g, b) for r, g, b in RGB256]
    return colors

def colormap(valuelist):
    #generate list of divergent colors
    clusters = numpy.unique(valuelist)
    N = len(clusters)
    #we will vary in 2 dimensions, so this is how many steps in each
    perm = int(math.ceil(math.sqrt(N)))
    #need offset, as human's can't tell colors that are unsaturated apart
    H = [(x*1.0/perm) for x in range(perm)]
    S = [(x*1.0/perm)+0.2 for x in range(perm)]
    #we will use this to truncate at correct length
    V = [0.7]*N
    #all combanations
    HS = itertools.product(H, S)
    H, S = zip(*HS)
    HSV = zip(H,S,V)
    RGB = [colorsys.hsv_to_rgb(h,s,v) for h, s, v in HSV]
    colorlookup = dict(zip(clusters, RGB[:N]))
    return colorlookup

def makeSeqFeaturesForGeneNeighbors(genename, clusterrunid, cur):
    outdata = getGeneNeighborhoods(genename, clusterrunid, cur)    
    genelocs = []
    for neargene in outdata: 
        neargeneid = neargene[1]
        start = int(neargene[4])
        stop = int(neargene[5])
        strandsign = neargene[6]
        if strandsign =='-': 
            strand = -1
        if strandsign =='+': 
            strand = +1
        feature = SeqFeature(FeatureLocation(start, stop), strand=strand, id = neargeneid)
        feature.qualifiers["cluster_id"] = int(neargene[8])
        genelocs.append(feature)
    return genelocs

def GetGeneToAlias():
    '''
    Get a dictionary from gene IDs to aliases (Note though that each gene can have more than one alias)
    '''
    geneToAlias ={}
    rootpath = os.path.split(os.path.split(locateDatabase())[0])[0]
    alias_file = os.path.join(rootpath, 'aliases', 'aliases')
    if not os.path.exists(alias_file):
        sys.stderr.write("WARNING: No aliases file found in expected location:%s\n" %(alias_file))
        return geneToAlias
    for line in open(alias_file, "r"):
        spl = line.strip('\r\n').split("\t")
        geneToAlias[spl[0]] = spl[1]
    return geneToAlias

def removeLeadingDashes(t):
    '''James could you explain this convention to me?'''
    tblastnadded = []
    for genename in t.get_leaf_names():
        #if it is from tblastn, we want to change it in the tree and have a record to indecate this in the plot
        if genename.startswith('-'): 
            genename = unsanitizeGeneId(genename)
            tblastn_leaf = t&genename
            genename = genename.lstrip('-')
            tblastnadded.append(genename)
            tblastn_leaf.name = genename
    return t, tblastnadded

def regionlength(genelocs):
    location = [(int(loc.location.start), int(loc.location.end)) for loc in genelocs]
    starts, ends = zip(*location)
    #have to compare both, as some are reversed
    start = max(max(starts),max(ends))
    end = min(min(starts),min(ends))
    return start, end

def make_region_drawing(genelocs, getcolor, centergenename, maxwidth):
    '''
    Makes a PNG figure for regions with a given color mapping, set of gene locations...
    
    TODO - Needs better documentation
    TODO make auto-del tempfiles, or pass svg as string
    '''

    geneToAlias = GetGeneToAlias()
    org, peg = splitrast(centergenename, removefigpeg = True)

    # Set up an entry genome diagram object
    gd_diagram = GenomeDiagram.Diagram("Genome Region")
    gd_track_for_features = gd_diagram.new_track(1, name="Annotated Features")
    gd_feature_set = gd_track_for_features.new_set()

    # Some basic properties of the figure itself
    arrowshaft_height = 0.3
    arrowhead_length = 0.3
    default_fontsize = 30 # Font size for genome diagram labels
    scale = 20     #AA per px for the diagram

    for feature in genelocs:
        bordercol=rcolors.white
        if feature.id == centergenename:
            bordercol=rcolors.red
            centerdstart, centerend = int(feature.location.start), int(feature.location.end)
            centerdstrand = feature.strand
        color = getcolor[feature.qualifiers["cluster_id"]]

        thisorg, thispeg = splitrast(feature.id, removefigpeg = True)
        gd_feature_set.add_feature(feature, name = thispeg,
                                   color=color, border = bordercol, 
                                   sigil="ARROW", arrowshaft_height=arrowshaft_height, arrowhead_length = arrowhead_length,
                                   label=False,  label_angle=20, label_size = default_fontsize
                                   )
    start, end = regionlength(genelocs)
    pagew_px = maxwidth / scale
    #offset so start of gene of interest lines up in all the figures
    midcentergene = abs(centerend - centerdstart)/2 + min(centerdstart, centerend)
    l2mid = abs(midcentergene - start)
    r2mid = abs(midcentergene - end)
    roffset = float((pagew_px/2) - (l2mid/scale))
    loffset = float((pagew_px/2) - (r2mid/scale))

    gd_diagram.draw(format="linear", start=start, end=end, fragments=1, pagesize=(225, pagew_px), xl=(loffset/pagew_px), xr=(roffset/pagew_px) )

    imgfileloc = "/tmp/%s%s_temp.png" %(org, peg)
    gd_diagram.write(imgfileloc, "PNG")
    #flip for reversed genes
    if centerdstrand == -1:
        os.system("convert -rotate 180 %s %s" % (imgfileloc, imgfileloc))
    return imgfileloc

def draw_tree_regions(clusterrunid, t, ts, cur, greyout=3):
    '''
    Draw the neighborhoods around each of the genes in a gene tree given the cluster and run IDs and the tree (t)

    clusterrunid is the run ID to use to identify homoloous clusters and ts is the treeStyle object associeted with the
    ETE tree t

    cur is a SQLite cursor object for the database

    The arrows are grayed out if less than "greyout" genes appear in a given cluster.
    '''
    t, tblastnadded = removeLeadingDashes(t)
    unsanitized = []
    for genename in t.get_leaf_names():
        unsanitized.append(unsanitizeGeneId(genename))

    # Create a list of SeqFeature objects for each gene in the tree
    regionindex={}
    for genename in unsanitized:
        sys.stderr.write("Getting gene neighborhoods for gene %s...\n" %(genename) )
        regionindex[genename]= makeSeqFeaturesForGeneNeighbors(genename, clusterrunid, cur)

    # Get a list of clusters containing these genes
    allclusters = []
    for gene in regionindex:
        for feature in regionindex[gene]:
            allclusters.append(feature.qualifiers["cluster_id"])

    uniqueclusters = set(allclusters)

    # Get clusters that have enough members to bother trying to color them (as determined by
    # the greyout keyword)
    multipleclusters = [c for c in uniqueclusters if allclusters.count(c) > greyout]
    getcolor = colormap(multipleclusters)

    #also add in grey (0.5,0.5,0.5 in RGB) for all others
    singleclusters = [c for c in uniqueclusters if allclusters.count(c) <= greyout]
    getcolor.update([(sc, (0.5,0.5,0.5)) for sc in singleclusters])

    #generate the region images for any leaf that has them, and map onto the tree
    #we will want to know the max width to make the figures
    widths = []
    for genelocs in regionindex.values():
        start, end = regionlength(genelocs)
        widths.append(abs(end - start))
    maxwidth = max(widths)
    for leaf in t.iter_leaves():
        newname = unsanitizeGeneId(leaf.name)
        # Not all genes necessarily are in the database and we don't want to crash if that happens.
        # Instead, Just don't print a neighborhood for them.
        try: 
            genelocs = regionindex[newname]
        except KeyError: 
            continue 
        sys.stderr.write("Making region drawing for gene ID %s...\n" %(newname))
        imgfileloc = make_region_drawing(genelocs, getcolor, newname, maxwidth)
        org, peg = splitrast(newname, removefigpeg = True)
        imageFace = faces.ImgFace(imgfileloc)
        leaf.add_face(imageFace, column=2, position = 'aligned')
        if newname in tblastnadded:
            leaf.add_face(TextFace("TBlastN added", fsize=30), column=3, position = 'aligned')
    #add legend for clusters
    ts = treelegend(ts, getcolor, greyout)
    return t, ts

def treelegendtext(cluster, color):
    text = TextFace(" %s " % cluster)
    text.hz_align = False
    text.fsize = 30
    text.fstyle = 'Bold'
    text.background.color = color
    return text

def treelegend(ts, getcolor, greyout):
    '''
    Add legend to the tree with cluster numbers corresponding to each color
    '''
    #needs hex, not 0 to 1 RGB, this function wants a list so unpack and pack  back up
    clusters, colors = zip(*getcolor.items())
    colors = RGB_to_hex(colors)
    colorlist = zip(clusters, colors)
    colorlist.sort()
    greynum = len([uc for uc, color in colorlist if color == '#7f7f7f'])
    # to make close to a box
    greycols =  int(math.ceil(math.sqrt(greynum))) 
    colornum = len(colorlist) - greynum
    # Because the color palette varies in 2 dimensions, this will make a box with the H and S indexed
    colorcols =  int(math.ceil(math.sqrt(colornum)))
    # Put legend with colors at bottom to display cluster IDs.
    ts.legend_position=1
    gnum = 0
    cnum = 0
    for cluster, color in colorlist:
        #drop grey ones
        if color == '#7f7f7f': 
            #offset the greys
            col = (gnum%greycols) + colorcols + 1 + 1 #offset from colors and grey def
            gnum += 1
        else:             
            col = cnum%colorcols
            cnum += 1
        text = treelegendtext(cluster, color)
        #placement of this legend
        ts.legend.add_face(text, column=col)
    ts.legend.add_face(treelegendtext("> %s occurrences        < or = %s occurrences " % (greyout, greyout),'#FFFFFF'), column=colorcols + 1)
    return ts

if __name__ == "__main__":
    '''
    Given a protein tree, make a figure contaiing that tree with neighborhoods overlayed onto the tree.
    '''

    usage="%prog -p protein_tree -r runid [options]"
    description="""Generates a tree with gene regions"""
    parser = optparse.OptionParser(usage=usage, description=description)
    parser.add_option("-r", "--runid", help="Run ID (required)", action="store", type="str", dest="clusterrunid", default=None)
    parser.add_option("-p", "--prottree", help="Protein tree (required)", action="store", type="str", dest="treeinfile", default=None)
    parser.add_option("-t", "--treetitle", help="Tree title (D: same as run ID)", action="store", type="str", dest="gene", default=None)
    parser.add_option("-o", "--outfile", help="Base name for output file (D: Same as input tree)", action="store", type="str", dest="outfile", default=None)
    parser.add_option("-d", "--display", help="Display the resulting tree (D: Don't display, just save)", action="store_true", dest="display", default=False)
    (options,args) = parser.parse_args()

    if options.treeinfile is None:
        sys.stderr.write("ERROR: -p (protein input tree) is required\n")
        exit(2)

    # James - how did you plan on using muliple clusters? Doesnt make a lot of sense to me.
    if options.clusterrunid is None:
        sys.stderr.write("ERROR - -r (runid) is a required argument\n")
        exit(2)

    if options.outfile is None:
        options.outfile = options.treeinfile

    clusterrunid = options.clusterrunid
    treeinfile = options.treeinfile 

    if options.gene is None:
        gene = options.clusterrunid
    else:
        gene = options.gene 

    # FIXME - if we declare this as a PhyloTree rather than just a Tree, ETE for some reason creates an
    # EXTRA COPY of the names (try it and see!) and therefore if you set show_leaf_names to True it shows
    # TWO copies of the leaf name at each leaf node (and if you set it to False it shows one, making it look
    # like it's ignoring you when you set it to False). The Tree class works as expected.
    #
    # Since we only need PhyloTree when we want to export to PhyloXML (not implemented yet), we should only create one
    # in that instance (due to the brokenness)
#    t = PhyloTree(treeinfile, sp_naming_function=parse_sp_name)
    t = Tree(treeinfile)
    ts= TreeStyle()
    
    con = sqlite3.connect(locateDatabase())
    cur = con.cursor()

    t, ts = draw_tree_regions(clusterrunid, t, ts, cur)

    # Now that we don't need to reference anything with the gene IDs any more, try to change them into
    # annotations
    for node in t.traverse():
        if node.is_leaf():
            unsanitized = unsanitizeGeneId(node.name)
            geneinfo = getGeneInfo( [ unsanitized ], cur)
            if len(geneinfo) > 0:
                organism = geneinfo[0][1]
                annotation = geneinfo[0][9]
            node.name = sanitizeString("%s_%s_%s" %(organism, annotation[0:63], unsanitized), False)

    t, ts = prettifyTree(t, title = gene + " cluster regions", show_bootstraps = False, ts=ts)

    os.system("rm test.svg 2> /dev/null")
    t.render("%s.svg" %(options.outfile), tree_style=ts)
    os.system("convert -trim -depth 32 -background transparent %s.svg %s.png" %(options.outfile, options.outfile))

    if options.display:
        t.show(tree_style=ts)

    con.close()
