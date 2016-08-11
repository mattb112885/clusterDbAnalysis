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
import shutil
import sys
import tempfile

from BioPythonGraphics import *
from FileLocator import *
from TreeFuncs import *
from ClusterFuncs import *
from sanitizeString import *

import numpy
from ete2 import Tree, faces, TreeStyle, TextFace, PhyloTree
from ete2 import Phyloxml, phyloxml


def regionlength(genelocs):
    location = [(int(loc.location.start), int(loc.location.end)) for loc in genelocs]
    starts, ends = list(zip(*location))
    #have to compare both, as some are reversed
    start = max(max(starts),max(ends))
    end = min(min(starts),min(ends))
    return start, end

def draw_tree_regions(clusterrunid, t, ts, cur, greyout=3, tempdir=None, label=False):
    '''
    Draw the neighborhoods around each of the genes in a gene tree given the cluster and run IDs and the tree (t)

    clusterrunid is the run ID to use to identify homologous clusters and ts is the treeStyle object associated with the
    ETE tree t

    cur is a SQLite cursor object for the database

    The arrows are grayed out if less than "greyout" genes appear in a given cluster.

    tempdir is a temporary directory in which to store the results. The user is responsible for deleting this
    directory afterwards to clean up.
    '''

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
    if len(list(seqfeatures.keys())) == 0:
        sys.stderr.write("WARNING: No genes in input tree had entries in the database so no neighborhoods will be drawn\n")
        return t, ts

    allseqfeatures = []
    for gene in seqfeatures:
        allseqfeatures += seqfeatures[gene]
    getcolor = makeClusterColorMap(allseqfeatures, greyout)

    #generate the region images for any leaf that has them, and map onto the tree
    #we will want to know the max width to make the figures
    widths = []
    for genelocs in list(seqfeatures.values()):
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
        imgfileloc = make_region_drawing(genelocs, getcolor, newname, maxwidth, tempdir=tempdir, label=label)
        imageFace = faces.ImgFace(imgfileloc)
        leaf.add_face(imageFace, column=2, position = 'aligned')

    #add legend for clusters
    ts = treelegend(ts, getcolor, greyout, clusterrunid, cur)

    return t, ts

def treelegendtext(whattoprint, color):
    text = TextFace(" %s " % whattoprint)
    text.hz_align = False
    text.fsize = 30
    text.fstyle = 'Bold'
    text.background.color = color
    return text

def treelegend(ts, getcolor, greyout, clusterrunid, cur):
    '''
    Add legend to the tree with cluster numbers corresponding to each color
    '''
    #needs hex, not 0 to 1 RGB, this function wants a list so unpack and pack  back up
    clusters, colors = list(zip(*list(getcolor.items())))
    colors = RGB_to_hex(colors)
    colorlist = list(zip(clusters, colors))
    colorlist.sort()
    greynum = len([uc for uc, color in colorlist if color == '#7f7f7f'])
    # to make close to a box
    greycols =  int(math.ceil(math.sqrt(greynum))) 
    colornum = len(colorlist) - greynum
    # Because the color palette varies in 2 dimensions, this will make a box with the H and S indexed
#    colorcols =  int(math.ceil(math.sqrt(colornum)))
    # Adding long annotations means we don't want this to be huge horizontally.
    colorcols = 1
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
            # For gray boxes we don't care what the function is.
            text = treelegendtext(cluster, color)
        else:             
            col = cnum%colorcols
            cnum += 1
            samplefunc = findRepresentativeAnnotation(clusterrunid, str(cluster), cur)
            text = treelegendtext(str(cluster) + ": " + samplefunc[0:63], color)
        #placement of this legend
        ts.legend.add_face(text, column=col)
    ts.legend.add_face(treelegendtext("> %s occurrences        < or = %s occurrences " % (greyout, greyout),'#FFFFFF'), column=colorcols + 1)
    return ts

if __name__ == "__main__":
    '''
    Given a protein tree, make a figure containing that tree with neighborhoods overlayed onto the tree.
    '''

    usage="%prog -p protein_tree -r runid [options]"
    description="""Generates a tree with gene regions"""
    parser = optparse.OptionParser(usage=usage, description=description)
    parser.add_option("-r", "--runid", help="Run ID (required)", action="store", type="str", dest="clusterrunid", default=None)
    parser.add_option("-p", "--prottree", help="Protein tree (required)", action="store", type="str", dest="treeinfile", default=None)
    parser.add_option("-t", "--treetitle", help="Tree title (D: same as run ID)", action="store", type="str", dest="gene", default=None)
    parser.add_option("-o", "--outfile", help="Base name for output file (D: Same as input tree)", action="store", type="str", dest="outfile", default=None)
    parser.add_option("-d", "--display", help="Display the resulting tree (D: Don't display, just save)", action="store_true", dest="display", default=False)
    parser.add_option("-c", "--cutoff", help="Number of members of a cluster below which a gene is greyed out (D: 3 - 2 or less are greyed out)", action="store", type="int", dest="cutoff", default=3)
    parser.add_option("-l", "--label", help="Add labels to the genes. The labels are the cluster IDs for the clusters in which the genes are found in the specified cluster run (D: Dont because its very messy)", action="store_true", dest="label", default=False)
    parser.add_option("--png", help="Save high-quality PNG and SVG images (with base name specified by -o or by default, with the same name as input file)", action="store_true", 
                      dest="savepng", default=False)
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

    tempdir = tempfile.mkdtemp()

    t, ts = draw_tree_regions(clusterrunid, t, ts, cur, greyout=options.cutoff, tempdir=tempdir, label=options.label)

    # Now that we don't need to reference anything with the gene IDs any more, try to change them into
    # annotations
    sanitizedToNot = getSanitizedContigList(cur)
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
                organism = ""
                try:
                    contig,start,stop = splitTblastn(unsanitized)
                    if contig in sanitizedToNot:
                        contig = sanitizedToNot[contig]
                    q = "SELECT organism FROM organisms INNER JOIN contigs ON contigs.organismid = organisms.organismid WHERE contigs.contig_mod=?;"
                    cur.execute(q, (contig,) )
                    for res in cur:
                        organism = res[0]
                except ValueError:
                    # Not a tblastn ID.
                    pass
            node.name = sanitizeString("%s_%s_%s" %(organism, annotation[0:63], unsanitized), False)
    
    t, ts = prettifyTree(t, title = gene + " cluster regions", show_bootstraps = False, ts=ts)

    if options.savepng:
        os.system("rm test.svg 2> /dev/null")
        t.render("%s.svg" %(options.outfile), tree_style=ts)
        os.system("convert -trim -depth 32 -background transparent %s.svg %s.png" %(options.outfile, options.outfile))

    if options.display:
        t.show(tree_style=ts)

    # Clean up
    con.close()
    try:
        sys.stderr.write("Deleting temporary directory %s...\n" %(tempdir))
        shutil.rmtree(tempdir)
    except OSError as exc:
        raise IOError("Unable to delete temporary directory %s" %(tempdir) )
