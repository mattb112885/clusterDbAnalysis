#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import optparse
from locateDatabase import *
from sanitizeString import *

####################################
# Lets read input arguments first.
####################################

usage = "%prog [options] Newick_file < runid or %prog -i runid [options] Newick_file . Default activity is to do nothing - one of -s, -p, or -d must be specified..."
description = "Draws gene context with specified number of genes as arrows with the appropraite direction (not to scale) alongside a phylogenetic tree. Only can draw arrows if the IDs in the tree agree with what is available in the database but does not crash when outgroups are also present (just nothing is drawn)."

parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-i", "--runid", help="Run id (default: read from stdin)", action="store", type="str", dest="runid", default=None)
parser.add_option("-n", "--neighborhood", help="Max number of genes away from target to display (D=3)", action="store", type="int", dest="MAXK", default=3)
parser.add_option("-a", "--annotation", help="Show annotations for neighboring genes underneath them (D: Show annotations in a legend underneath the tree)", action="store_true", dest="annotation", default=False)
parser.add_option("-d", "--display", help="Display result", action="store_true", dest="display", default=False)
parser.add_option("-s", "--savesvg", help="Convert the file to svg (requires -b)", action="store_true", dest="savesvg", default=False)
parser.add_option("-p", "--savepng", help="Convert the file to png (requires -b, implies -s)", action="store_true", dest="savepng", default=False)
parser.add_option("-b", "--basename", help="Base name for file outputs (ignored without -s or -p)", action="store", type="str", dest="basename", default=None)
parser.add_option("-r", "--rootgene", help="Root on this gene (default = keep same root as nwk file).", action="store", type="str", dest="rootgene", default=None)
(options, args) = parser.parse_args()

# Must specify a newick file
if len(args) < 1:
    sys.stderr.write("ERROR: Newick file must be specified. Use -h for help details\n")
    exit(2)

# Make sure the MAXK specified is possible given what we actually calculated to put into the database.
MAXK = options.MAXK
if MAXK > 5:
    sys.stderr.write("ERROR: Number of genes away must be <= 5 (this is the max that was pre-calculated)\n")
    exit(2)
if MAXK <= 0:
    sys.stderr.write("ERROR: Invalid number of genes away (must be more than 0)")
    exit(2)

# At least one of -p, -s, and -d must be specified
if not ( options.savesvg or options.savepng or options.display):
    sys.stderr.write("ERROR: At least one of -p, -s, or -d must be specified. Use -h for help details\n")
    exit(2)

# Must specify a base file name if you want to save svg or png
if ( options.savesvg or options.savepng ) and options.basename == None:
    sys.stderr.write("ERROR: Calling with -p (savepng) or -s (savesvg) requires specification of a base filename to save to\n")
    exit(2)

# Since we need SVG to get PNG the savepng option overrides the savesvg option
if options.savepng:
    options.savesvg = True

####################################
# Import lots of stuff for drawing...
####################################

from PyQt4 import QtCore
from PyQt4.QtCore import QPointF
from PyQt4.QtGui import QGraphicsRectItem, QGraphicsSimpleTextItem, \
    QGraphicsPolygonItem, QColor, QPen, QBrush, QPolygonF
from ete2 import Tree, faces, TreeStyle, NodeStyle, AttrFace
# And other stuff we need...
from operator import itemgetter
import fileinput
import sqlite3
import os

##################################################
# Create a QT object that is shaped like an arrow
# given its dimensions, directions, and color.
##################################################
def makeArrowNode(node, *args, **kargs):
    # node: An ETE Node object
    # Args: (arrowlength, arrowwidth, rectanglelength, rectanglewidth, direction (+ or -), annotation (gene ID), color, brushstyle)
    arrowLength = args[0][0]
    arrowWidth = args[0][1]
    rectLength = args[0][2]
    rectWidth = args[0][3]
    direction = args[0][4]
    annotation = args[0][5]
    hexcolor = args[0][6]

    ## Creates a main master Item that will contain all other elements
    masterItem = QGraphicsRectItem(0, 0, 40+rectLength, 40+2*arrowWidth)
    # Keep a link within the item to access node info
    masterItem.node = node 
    # I dont want a border around the masterItem
    masterItem.setPen(QPen(QtCore.Qt.NoPen))

    # Draw arrow in this order
    #              x         ^
    #              3\   arrowWidth
    #   x    2     x 4       |         ^
    #   1             x      V     rectWidth
    #   x   7      x 5                 V
    #              6/
    #              x
    #              <-> arrowLength
    #   <-rectLength->
    
    # Calculate points for "+" orientation first
    Sx0 = 20
    Sy0 = 20 + arrowWidth/2
    
    Sx1 = Sx0
    Sy1 = Sy0 + rectWidth
    
    Sx2 = Sx1 + (rectLength - arrowLength)
    Sy2 = Sy1

    Sx3 = Sx2
    Sy3 = Sy2 + (arrowWidth - rectWidth/2 )
    
    Sx4 = Sx3 + arrowLength
    Sy4 = Sy3 - arrowWidth
    
    Sx5 = Sx4 - arrowLength
    Sy5 = Sy4 - arrowWidth

    Sx6 = Sx5
    Sy6 = Sy5 + (arrowWidth - rectWidth/2)

    # Should be back at the start now.
    Sx7 = Sx6 - (rectLength - arrowLength)
    Sy7 = Sy6 

    # Deal with left-facing arrows. This is just a reflection...
    maxl = 40 + rectLength
    if direction == "-":
        Sx0 = maxl - Sx0
        Sx1 = maxl - Sx1
        Sx2 = maxl - Sx2
        Sx3 = maxl - Sx3
        Sx4 = maxl - Sx4
        Sx5 = maxl - Sx5
        Sx6 = maxl - Sx6
        Sx7 = maxl - Sx7

    # Actually draw the arrow
    points = [ QPointF(Sx0, Sy0),  
               QPointF(Sx1, Sy1),
               QPointF(Sx2, Sy2),
               QPointF(Sx3, Sy3),
               QPointF(Sx4, Sy4),
               QPointF(Sx5, Sy5),
               QPointF(Sx6, Sy6),
               QPointF(Sx7, Sy7) ]
    poly = QPolygonF(points)

    # Add the arrow to the master object.
    arrow = QGraphicsPolygonItem(poly, masterItem)
    pen = QPen(QtCore.Qt.black, 2, QtCore.Qt.SolidLine)
    arrow.setPen(pen)
    
    # Change arrow's color
    arrow.setBrush(QBrush(QColor(hexcolor)))

    # Add annotation (gene ID) to the arrow
    text = QGraphicsSimpleTextItem(annotation)
    text.setParentItem(arrow)
    text.setPen(QPen(QPen(QColor("black"), 2, QtCore.Qt.SolidLine)))

    # Center text according to masterItem size
    tw = text.boundingRect().width()
    th = text.boundingRect().height()
    center = masterItem.boundingRect().center()
    text.setPos(center.x()-tw/2, center.y()-th/2)
    
    return masterItem

# Main portion of function
colorTable = [ 
    "#000000", # Black
    "#808080", # Dark gray
    "#808000", # Olive
    "#B8860B", # DarkGoldenRod
    "#008080", # Teal
    "#008000", # Green
    "#FF00FF", # Magenta
    "#7FFFD4", # Aquamarine
    "#FF0000", # Red
    "#0000FF", # Blue
    "#800080", # Purple
    "#800000", # Maroon
    "#000080", # Navy blue
    "#A52A2A", # Brown
    "#9932CC", # Dark orchid
    "#FF4500", # Orange-red
    "#FF1493", # Deep pink
    "#00FFFF", # Cyan
    "#E9967A", # Dark salmon
    "#FFD700", # Gold
    "#FFFFFF"  # White (for nothing there)
    ]

# Read Newick file
sys.stderr.write("Reading tree file...\n")
t = Tree(args[0])

# If outgroup is specified, re-root now before doing anything else.
# This will just fail if the specified protein isn't present in the tree.
if not options.rootgene == None:
    done = False
    for node in t.traverse():
        if node.name == options.rootgene:
            t.set_outgroup(node)
            done = True
            break
    if not done:
        sys.stderr.write("ERROR: Specified outgroup %s not found in tree\n" %(options.rootgene))
        exit(2)

##############################################
# Get various gene / organism / annotation / neighborhood / cluster info out of the database
##############################################

# Annotation and organism
geneToAnnote = {}
geneToOrganism = {}
sys.stderr.write("Reading gene annotations and organisms from database...\n")

con = sqlite3.connect(locateDatabase())
cur = con.cursor()
cur.execute("SELECT * FROM processed;")

for l in cur:
    spl = [ str(s) for s in list(l) ]
    geneToAnnote[spl[0]] = sanitizeString(spl[9], False)
    geneToOrganism[spl[0]] = sanitizeString(spl[1], False)

# Neighborhoods
sys.stderr.write("Pulling gene neighborhoods out of the database...\n")

cur.execute("SELECT * FROM neighborhoods;")
geneToNeighbors = {}
for l in cur:
    spl = [ str(s) for s in list(l) ]
    spl[2] = int(spl[2])
    if abs(spl[2]) > MAXK:
        continue
    # A neighborPair is ( neighboring gene, number of genes away, strand )
    neighborPair = [ spl[1], spl[2], spl[5] ]
    if spl[0] in geneToNeighbors:
        geneToNeighbors[spl[0]].append(neighborPair)
    else:
        geneToNeighbors[spl[0]] = [ neighborPair ]

# We will need to make all of the genes have the right number of neighbors and
# add fillers for those that don't.
sys.stderr.write("Filling in neighbors for genes on the ends of contigs...\n")
for gene in geneToNeighbors:
    allPairs = geneToNeighbors[gene]
    for j in range(-MAXK, MAXK+1):
        ok = False
        for pair in allPairs:
            if pair[1] == j:
                ok = True
                break
        if not ok:
            neighborPair = [ "NONE", j, "+" ]
            geneToNeighbors[gene].append( neighborPair )

# If the center gene is on the "-" strand, reverse the direction of everything
# so that we are comparing apples to apples.
for gene in geneToNeighbors:
    allPairs = geneToNeighbors[gene]
    selfPair = []
    for pair in allPairs:
        if pair[1] == 0:
            selfPair = pair
            break
    if selfPair[2] == "-":
        newPairList = []
        for pair in allPairs:
            pair[1] = -pair[1]
            if pair[2] == "-":
                pair[2] = "+"
            else:
                pair[2] = "-"
            newPairList.append(pair)
        geneToNeighbors[gene] = newPairList

# Sort all of the neighborPair lists by k (number of genes away)
# This ensures that the genes are actually listed in the same order!
for gene in geneToNeighbors:
    sortTuples = sorted(tuple(geneToNeighbors[gene]), key=itemgetter(1))
    geneToNeighbors[gene] = sortTuples

################
# Get the run ID to search for
# and all the cluster info for that run ID
################

sys.stderr.write("Gathering cluster information for specified run ID...\n")
done = False
runid = options.runid
if runid == None:
    for line in fileinput.input("-"):
        if done == True:
            sys.stderr.write("ERROR: Must only pipe in ONE run ID\n")
            exit(2)
        done = True
        runid = line.strip('\r\n')

cur.execute("SELECT * FROM clusters WHERE clusters.runid=?;", (runid, ) )
atLeastOne = False
geneToCluster = {}
for l in cur:
    atLeastOne = True
    geneToCluster[str(l[2])] = str(l[1])

if not atLeastOne:
    sys.stderr.write("ERROR: Specified run ID %s not found in database\n" %(runid))
    exit(2)

############
# Set up cluster list to have the most-prevalent clusters colored first
############

# Count instances of each cluster and obtain a representative annotation for each one
clusterToNumber = {}
clusterToAnnote = {}
for node in t.traverse():
    if node.is_leaf():
        # We don't want to just crash if we have genes that don't have data in our database
        # (e.g. if we include an outgroup in our tree)
        if not node.name in geneToNeighbors:
            sys.stderr.write("WARNING: Gene %s in the tree was not present in the database - it will be ignored and not have neighbors specified! \n" %(node.name) )
            continue
        pairList = geneToNeighbors[node.name]
        for pair in pairList:
            if not pair[0] in geneToCluster:
                # Don't bother printing warnings about fake genes...
                if not pair[0] == "NONE":
                    sys.stderr.write("WARNING: The tree contains gene %s that is not present in the clustering run specified. These will not be colored!\n" %(pair[0]) )
                continue
            cluster = geneToCluster[pair[0]]
            clusterToAnnote[cluster] = geneToAnnote[pair[0]]
            if cluster in clusterToNumber:
                clusterToNumber[cluster] += 1
            else:
                clusterToNumber[cluster] = 1

# Assign color to each cluster
# according to the number of genes in that cluster within the neighbors
# of genes in the tree
clusterTuples = sorted(clusterToNumber.iteritems(), key=itemgetter(1), reverse=True)
clusterToColor = {}
currentIdx = 0
for tup in clusterTuples:
    clusterToColor[tup[0]] = colorTable[currentIdx]
    if currentIdx < len(colorTable) - 1:
        currentIdx += 1

clusterToColor["NONE"] = "#FFFFFF"
clusterToAnnote["NONE"] = "NONE or OTHER"
geneToCluster["NONE"] = "NONE"
geneToAnnote["NONE"] = "NONE"

#############
# Add actual arrow objects to the tree
#############
sys.stderr.write("Adding arrow objects to leaves of the tree...\n")
for node in t.traverse():
    if node.is_leaf():
        # Dont' crash because of e.g. outgroups put in. We already warned about this so don't need to do it again.
        if not ( node.name in geneToNeighbors and node.name in geneToOrganism and node.name in geneToAnnote ):
            # We still want to make the text face though so that it is visible!
            F = faces.TextFace(node.name, ftype="Times", fsize=32)
            node.add_face(F, 0, position="aligned")
            continue

        # Add an annotation text with larger font to replace the crappy size-10 ish font that comes by default...
        newname = "_".join( [ node.name, geneToOrganism[node.name], geneToAnnote[node.name] ] )
        F = faces.TextFace(newname, ftype="Times", fsize=32)
        node.add_face(F, 0, position="aligned")

        genename = node.name
        geneNeighbors = geneToNeighbors[genename]
        for neighborPair in geneNeighbors:
            label = neighborPair[0]
            color = clusterToColor[geneToCluster[neighborPair[0]]]
            direction = neighborPair[2]
            # Create arrow object and add to tree
            F = faces.DynamicItemFace(makeArrowNode,30,30,300,30, direction, label, color)
            node.add_face(F, neighborPair[1] + MAXK+1, position="aligned")
            # If desired, add annotations below the arrows.
            if options.annotation:
                anno = geneToAnnote[label]
                trimmedAnno = anno[0:40] + "\n" + anno[40:80]
                F = faces.TextFace(trimmedAnno, ftype="Times", fsize=20)
                node.add_face(F, neighborPair[1] + MAXK + 1, position="aligned")
    else:
        # Make the branch support bigger
        F = faces.TextFace(node._support, ftype="Times", fsize=32)
        node.add_face(F, 0, position="branch-top")


ts = TreeStyle()

#####################
# Make legend for arrow colors
#####################
currentcol = 0
for cluster in clusterToAnnote:
    # Limit the number of annotations in a row but let it be more than 1 since otherwise it's a huge waste of space.
    if currentcol > 7:
        currentcol = 0
    # Color box
    fc1 = faces.TextFace("")
    fc1.background.color = clusterToColor[cluster]
    boxsize = 20
    fc1.margin_left = boxsize
    fc1.margin_right = boxsize
    ts.legend.add_face(fc1, column=currentcol)
    # Annotation box
    fc2 = faces.TextFace(clusterToAnnote[cluster], ftype="Times", fsize=20)
    ts.legend.add_face(fc2, column=currentcol + 1)
    currentcol += 2

######################
# Other setup \ tree prettying
######################

# We already made big bootstrap labels so don't bother showing the crappy small ones.
# Same with the leaf labels.
ts.show_branch_support = False
ts.show_leaf_name = False

# Estimate a reasonable tree width according to the maximum distance on the tree
maxdist = 0
root = t.get_tree_root()
for node in t.traverse():
    if node.is_leaf():
        dist = t.get_distance(node, root)
        if dist > maxdist:
            maxdist = dist
# Old versions of ETE don't support the tree_width specifier
try:
    ts.tree_width = maxdist * 20
except ValueError:
    sys.stderr.write("WARNING: Your version of ETE does not support the tree_width specification - consider upgrading to get nicer looking trees!\n")

# If the trees are identical and identically rooted, we want them to always appear the same.
# Start by sorting by number of branches (ladderize)
#
# I then break ties according to the names of leaves descended from a given node.
# Essentially this amounts to sorting first by number of branches and then by alphabetical order
# This SHOULD always give a consistent result for identical trees (regardless of how they appear in the file).
t.ladderize(direction=0)

for node in t.traverse(strategy="levelorder"):
    if not node.is_leaf():
        children = node.get_children()
        if not len(children) == 2:
            sys.stderr.write("INTERNAL ERROR: Should always have 2 children per node?\n")
            continue;
        nl0 = len(children[0].get_leaves())
        nl1 = len(children[1].get_leaves())
        if nl0 == nl1:
            names0 = "".join(sorted(children[0].get_leaf_names()))
            names1 = "".join(sorted(children[1].get_leaf_names()))
            if names0 > names1:
                node.swap_children()

################
# Outputs
################

if options.savesvg:
    # Some versions of ETE create (and do not delete) a "test.svg" file automatically
    # and others do not.
    # To avoid confusion, lets remove that one and make our own (using the specified tree style)...
    os.system("rm test.svg 2> /dev/null")
    t.render("%s.svg" %(options.basename), tree_style=ts)

# Convert the svg file into a high-quality (300 dpi) PNG file...
# The PNG converter in ETE gives a terrible quality image
#
# Use convert to make something better and then trim the edges off.
if options.savepng:
    os.system("convert -trim -depth 16 -background transparent %s.svg %s.png" %(options.basename, options.basename))

if options.display:
    t.show(tree_style=ts)

con.close()
