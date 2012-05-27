#!/usr/bin/python
# -*- coding: utf-8 -*-

# Lets read input arguments first.
import sys

if not len(sys.argv) == 3 and not len(sys.argv) == 4:
    print "Usage: ./TutorialArrow.py [Newick File] [Output basename] [number of genes away (optional)]"
    print "Pipe in a run ID [only one!] to use for coloring."
    print "Output file is just a file name"
    print "Usage of this function requlres you to have a UI (i.e. you must have an X server up and running)"
    print "Number of genes away must be less than or equal to 5"
    exit(2)

if len(sys.argv) == 4:
    MAXK = int(sys.argv[3])
    print MAXK
    if MAXK > 5:
        print "ERROR: Number of genes away must be <= 5"
        exit(2)
    if MAXK <= 0:
        print "ERROR: Invalid maximum number of genes away"

# Import lots of stuff for drawing...
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

# Create a QT object that is shaped like an arrow
# given its dimensions, directions, and color.
def makeArrowNode(node, *args, **kargs):
    arrowLength = args[0][0]
    arrowWidth = args[0][1]
    rectLength = args[0][2]
    rectWidth = args[0][3]
    # If direction is "-" we just do rectLength - all of the coordinates...
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

    points = [ QPointF(Sx0, Sy0),  
               QPointF(Sx1, Sy1),
               QPointF(Sx2, Sy2),
               QPointF(Sx3, Sy3),
               QPointF(Sx4, Sy4),
               QPointF(Sx5, Sy5),
               QPointF(Sx6, Sy6),
               QPointF(Sx7, Sy7) ]
    poly = QPolygonF(points)

    # Set up a drawing object
    arrow = QGraphicsPolygonItem(poly, masterItem)
    pen = QPen(QtCore.Qt.black, 2, QtCore.Qt.SolidLine)
    arrow.setPen(pen)
    
    # Change arrow's color
    arrow.setBrush(QBrush(QColor(hexcolor)))

    # Add annotation to the arrow
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
t = Tree(sys.argv[1])

# Read geneinfo file (should include ALL of the genes across all the organisms)
# We need it to include all the genes because we need annotations not only for the genes
# on the tree but also for their neighbors.
geneToAnnote = {}
geneToOrganism = {}
sys.stderr.write("Reading gene annotations and organisms from database...\n")

con = sqlite3.connect("db/methanosarcina")
cur = con.cursor()
cur.execute("SELECT * FROM processed;")

for l in cur:
    spl = [ str(s) for s in list(l) ]
    geneToAnnote[spl[0]] = spl[9]
    geneToOrganism[spl[0]] = spl[1]

# Read neighborhood file - NOTE that we want BOTH strands
# rather than just neighbors on the same strand for this
# analysis. Make sure you call the correct function to generate
# the neighborhood file!
#
# Also this should be a concatinated file.
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
runid = ""
for line in fileinput.input("-"):
    if done == True:
        sys.stderr.write("ERROR: Must only pipe in ONE run ID\n")
        exit(2)
    done = True
    runid = line.strip()

cur.execute("SELECT * FROM clusters WHERE clusters.runid=?;", (runid, ) )
atLeastOne = False
geneToCluster = {}
for l in cur:
    atLeastOne = True
    geneToCluster[str(l[2])] = str(l[1])

if not atLeastOne:
    sys.stderr.write("ERROR: piped in run ID not found in database\n")
    exit(2)

############
# Set up cluster list to have the most-used ones first
############

# Count instances of each cluster
clusterToNumber = {}
for node in t.traverse():
    if node.is_leaf():
        pairList = geneToNeighbors[node.name]
        for pair in pairList:
            if not pair[0] in geneToCluster:
                # Don't bother printing warnings about fake genes...
                if not pair[0] == "NONE":
                    sys.stderr.write("WARNING: The tree contains gene %s that is not present in the clustering run specified. These will not be colored!\n" %(pair[0]) )
                continue
            cluster = geneToCluster[pair[0]]
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
    print tup
    clusterToColor[tup[0]] = colorTable[currentIdx]
    if currentIdx < len(colorTable) - 1:
        currentIdx += 1

clusterToColor["NONE"] = "#FFFFFF"
geneToCluster["NONE"] = "NONE"
geneToAnnote["NONE"] = "NONE"

#############
# Add actual arrow objects to the tree
#############
sys.stderr.write("Adding arrow objects to leaves of the tree...\n")
for node in t.traverse():
    if node.is_leaf():
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
            # For now I make all the arrows the same length.
            F = faces.DynamicItemFace(makeArrowNode,30,30,300,30, direction, label, color)
            node.add_face(F, neighborPair[1] + MAXK+1, position="aligned")
            anno = geneToAnnote[label]
            trimmedAnno = anno[0:30] + "\n" + anno[30:60]
            F = faces.TextFace(trimmedAnno, ftype="Times", fsize=20)
            node.add_face(F, neighborPair[1] + MAXK + 1, position="aligned")

# Note - PDF and PNG export are both horrid-quality - I will try to fix this but for now I'll just export to SVG...
# My suggestion is to find something that doesn't crash and that isn't horrible quality to convert this to another format.
# Is there such a thing though??? Sigh...
ts = TreeStyle()
# We DO want to show bootstraps
ts.show_branch_support = True
# We'll be putting thse in separately
ts.show_leaf_name = False
# The default width of the tree is too squished.
ts.tree_width = 1000
t.render("%s.svg" %(sys.argv[2]), tree_style=ts)

# Convert the svg file into a high-quality (300 dpi) PNG file...
# The PNG converter in ETE gives a terrible quality image
# as does the "convert" function (which is probably what ETE uses)
# so this is the best I could come up with...
os.system("inkscape -e %s_temp.png -d 300 %s.svg" %(sys.argv[2], sys.argv[2]) )
# Then trim off the edges
os.system("convert -trim %s_temp.png %s.png" %(sys.argv[2], sys.argv[2]))
os.system("rm %s_temp.png" %(sys.argv[2]))

# This is an interactive part. Uncomment this if you want interaction with the tree.
# As of right now the interaction doesn't DO anything. I will consider making it a hover-over thing
t.show(tree_style=ts)

con.close()
