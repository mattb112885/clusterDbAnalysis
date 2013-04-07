#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
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
'''

import math
import itertools
import optparse
import os
import sqlite3
import sys

from BioPythonGraphics import *
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

#######################################
## DEPRECIATED ########################
#######################################
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

###################################################
### NEIGHBORHOOD HANDLING for genes in database ###
###################################################

def makeSeqFeaturesForGeneNeighbors(genename, clusterrunid, cur):
    '''
    Given a gene ID, compute the neighbors of that gene and create a SeqFeature
    object for each.

    Returns a list of neigboring genes or an empty array [] if we couldn't do anything
    with the gene Id given (not found in database or similar issue)
    '''
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

####################################################
#### Functions for TBLASTN neighborhood support ####
####################################################
def getSanitizedContigList(cur):
    '''
    Get a list of sanitized Contig IDs. Returns a dictionary from sanitized to unsanitized contig IDs
    present in the database.
    '''
    q = "SELECT DISTINCT contig_mod FROM contigs;"
    cur.execute(q)
    sanitizedToNot = {}
    for res in cur:
        sanitizedToNot[sanitizeString(res[0], False)] = res[0]
    return sanitizedToNot

def makeSeqObjectsForTblastnNeighbors(tblastn_id, clusterrunid, cur, N=200000):
    '''
    Given a tBBLSATn ID and a dictionary from sanitized contig IDs (which is what will be
    present in the TBLASTN id) to non-sanitized IDs (which are what is in the database),
    returns a list of seq objects INCLUDING the TBLASTN hit itself (so that we can show that
    on the region drawing).

    We pick an N large enough to get at least one gene and then pick the closest one and get
    all of its neighbors with a call to makeSeqFeaturesForGeneNeighbors() and just tack the TBLASTN
    onto it.
    '''
    # Lets first get the contig and start/stop locations (which tell us teh strand) out of
    # the TBLASTN id. This returns a ValueError if it fails which the calling function can catch if needed.
    contig,start,stop = splitTblastn(tblastn_id)
    start = int(start)
    stop = int(stop)

    # Create a seq object for the TBLASTN hit itself
    if start < stop:
        strand = +1
    else:
        strand = -1
    tblastn_feature = SeqFeature(FeatureLocation(start, stop), strand=strand, id=tblastn_id)
    tblastn_feature.qualifiers["cluster_id"] = -1
    
    # Find the neighboring genes.
    neighboring_genes = getGenesInRegion(contig, start-N, stop+N, cur)
    if len(neighboring_genes) == 0:
        sys.stderr.write("WARNING: No neighboring genes found for TBLASTN hit %s within %d nucleotides" %(tblastn_id, N))
        return feature
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

###########################
#### Drawing functions ####
###########################

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

    imgfileloc = "/tmp/%s.png" %(sanitizeString(centergenename, False))
    
    # Set up an entry genome diagram object
    gd_diagram = GenomeDiagram.Diagram("Genome Region")
    gd_track_for_features = gd_diagram.new_track(1, name="Annotated Features")
    gd_feature_set = gd_track_for_features.new_set()

    # Some basic properties of the figure itself
    arrowshaft_height = 0.3
    arrowhead_length = 0.3
    default_fontsize = 30 # Font size for genome diagram labels
    scale = 20     #AA per px for the diagram

    # Build arrow objects for all of our features.
    for feature in genelocs:
        bordercol=rcolors.white
        if feature.id == centergenename:
            bordercol=rcolors.red
            centerdstart, centerend = int(feature.location.start), int(feature.location.end)
            centerdstrand = feature.strand
        color = getcolor[feature.qualifiers["cluster_id"]]

        
        gd_feature_set.add_feature(feature, name = feature.id,
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

    # DEPRECIATED
    t, tblastnadded = removeLeadingDashes(t)

    unsanitized = []
    for genename in t.get_leaf_names():
        unsanitized.append(unsanitizeGeneId(genename))


    # Create a list of SeqFeature objects for the neighbors of each gene in the tree
    # If passed a TBLASTN hit it will create seq objects for every gene surrounding the TBLASTN hit and
    # for the TBLASTN hit itself.
    #
    # Nothing is added if we can't find that ID in the database or the ID is badly formatted.
    seqfeatures={}
    for genename in unsanitized:
        sys.stderr.write("Getting gene neighborhoods for gene %s...\n" %(genename) )
        features_for_genename = makeSeqFeaturesForGeneNeighbors(genename, clusterrunid, cur)
        if len(features_for_genename) > 0:
            seqfeatures[genename] = features_for_genename
        else:
            # Try TBLASTN and if that doesn't work, just give up.
            try:
                features_for_tblastn = makeSeqObjectsForTblastnNeighbors(genename, clusterrunid, cur)
                seqfeatures[genename] = features_for_tblastn
            except ValueError:
                sys.stderr.write("WARNING: Unable to find entries for gene or TBLASTN hit %s in database\n" %(genename) )
                pass

    # Don't bother trying the rest if nothing matches at all.
    if len(seqfeatures.keys()) == 0:
        sys.stderr.write("WARNING: No genes in input tree had entries in the database so no neighborhoods will be drawn\n")
        return t, ts

    # Get a list of clusters containing these genes
    allclusters = []
    for gene in seqfeatures:
        for feature in seqfeatures[gene]:
            allclusters.append(feature.qualifiers["cluster_id"])

    uniqueclusters = set(allclusters)

    # Get clusters that have enough members to bother trying to color them (as determined by
    # the greyout keyword)
    multipleclusters = [c for c in uniqueclusters if allclusters.count(c) > greyout]

    # Don't die if nothing has enough clusters...
    if len(multipleclusters) > 0:
        getcolor = colormap(multipleclusters)
    else:
        getcolor = {}

    #also add in grey (0.5,0.5,0.5 in RGB) for all others
    singleclusters = [c for c in uniqueclusters if allclusters.count(c) <= greyout]
    getcolor.update([(sc, (0.5,0.5,0.5)) for sc in singleclusters])

    #generate the region images for any leaf that has them, and map onto the tree
    #we will want to know the max width to make the figures
    widths = []
    for genelocs in seqfeatures.values():
        start, end = regionlength(genelocs)
        widths.append(abs(end - start))
    maxwidth = max(widths)

    for leaf in t.iter_leaves():
        newname = unsanitizeGeneId(leaf.name)
        # Not all genes necessarily are in the database and we don't want to crash if that happens.
        # Instead, Just don't print a neighborhood for them.
        try: 
            genelocs = seqfeatures[newname]
        except KeyError: 
            continue 
        sys.stderr.write("Making region drawing for gene ID %s...\n" %(newname))
        imgfileloc = make_region_drawing(genelocs, getcolor, newname, maxwidth)
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
            else:
                # FIxME - Attempt to get the organism name from the contig for TBLASTN IDs. If that fails we just give up.
                annotation = ""
                try:
                    contig,start,stop = splitTblastn(unsanitized)
                    q = "SELECT organism FROM organisms INNER JOIN contigs ON contigs.organismid = organisms.organismid WHERE contigs.contig_mod=?;"
                    cur.execute(q, (contig,) )
                    for res in cur:
                        organism = res[0]
                except ValueError:
                    organism = ""
            node.name = sanitizeString("%s_%s_%s" %(organism, annotation[0:63], unsanitized), False)
    
    t, ts = prettifyTree(t, title = gene + " cluster regions", show_bootstraps = False, ts=ts)

    os.system("rm test.svg 2> /dev/null")
    t.render("%s.svg" %(options.outfile), tree_style=ts)
    os.system("convert -trim -depth 32 -background transparent %s.svg %s.png" %(options.outfile, options.outfile))

    if options.display:
        t.show(tree_style=ts)

    con.close()
