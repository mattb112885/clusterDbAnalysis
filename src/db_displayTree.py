#!/usr/bin/python
# -*- coding: utf-8 -*-

# Display a tree or reroot.
#
# Comes with a decent set of options.

import sys
import optparse
import os
from locateDatabase import *
from sanitizeString import *

####################################
# Lets read input arguments first.
####################################
usage = "%prog (-s -b basename|-p -b basename|-n -b basename|-d) [options] Newick_file"
description="""Display a tree with annotations and specified root and formats.
There is no default acticity; one of -s, -p, -n, or -d must be specified.."""
parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-d", "--display", help="Display result", action="store_true", dest="display", default=False)
parser.add_option("-s", "--savesvg", help="Convert the file to svg (requires -b)", action="store_true", dest="savesvg", default=False)
parser.add_option("-p", "--savepng", help="Convert the file to png (requires -b, implies -s)", action="store_true", dest="savepng", default=False)
parser.add_option("-n", "--savenewick", help="Save re-rooted tree as a newick file (requires -b)", action="store_true", dest="savenewick", default=False)
parser.add_option("-b", "--basename", help="Base name for file outputs (ignored without -s or -p)", action="store", type="str", dest="basename", default=None)
parser.add_option("-r", "--rootgene", help="Root on this gene (default = keep same root as nwk file).", action="store", type="str", dest="rootgene", default=None)
parser.add_option("-o", "--rootorg", help="Root on this organism ID (e.g. 83333.1) (default = keep same root as nwk file)", action="store", type="str", dest="rootorg", default=None)
(options, args) = parser.parse_args()

# Must specify a newick file
if len(args) < 1:
    sys.stderr.write("ERROR: Newick file must be specified. Use -h for help details\n")
    exit(2)

# At least one of -p, -s, and -d must be specified
if not ( options.savesvg or options.savepng or options.display or options.savenewick):
    sys.stderr.write("ERROR: At least one of -p, -s, -n or -d must be specified. Use -h for help details\n")
    exit(2)

# Must specify a base file name if you want to save a file
if ( options.savesvg or options.savepng or options.savenewick ) and options.basename == None:
    sys.stderr.write("ERROR: Calling with -p (savepng) or -s (savesvg) requires specification of a base filename to save to\n")
    exit(2)

# Can only specify one of -o and -r
if ( not options.rootgene == None ) and ( not options.rootorg == None ):
    sys.stderr.write("ERROR: Can only specify one of -o and -r (rooting options)\n")
    exit(2)

# Since we need SVG to get PNG the savepng option overrides the savesvg option
if options.savepng:
    options.savesvg = True

####################################
# Import lots of stuff for drawing...
####################################

from ete2 import Tree, faces, TreeStyle, NodeStyle, AttrFace
import fileinput
import sqlite3
from TreeFuncs import *

# Read Newick file
sys.stderr.write("Reading tree file...\n")
t = Tree(args[0])

# If outgroup is specified, re-root now before doing anything else.
# This will just fail if the specified protein isn't present in the tree.
rerootEteTree(t, root_leaf = options.rootgene, root_leaf_part = options.rootorg)

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
    # The SVG parser whines with some special characters (i.e. ' )
    geneToAnnote[spl[0]] = sanitizeString(spl[9], False)
    geneToOrganism[spl[0]] = sanitizeString(spl[1], False)

######################
# Add annotations and
# larger bootstrap values to tree
######################

for node in t.traverse():
    if node.is_leaf():
        # Dont' crash because of e.g. outgroups put in. We already warned about this so don't need to do it again.
        if node.name in geneToOrganism and node.name in geneToAnnote:
            newname = "_".join( [ node.name, geneToOrganism[node.name], geneToAnnote[node.name] ] )
            node.name = newname

# Standardize font sizes and tree width
t, ts = prettifyTree(t)
# Standardize leaf order in equivalent trees (with same root)
t = standardizeTreeOrdering(t)

# Label face columns [This will be useful for labeling tables next to the tree!]
F = faces.TextFace("Annotation", ftype="Times", fsize=20)
ts.aligned_header.add_face(F, 0)

if options.savenewick:
    t.write(outfile="%s.nwk" %(options.basename), format=0)

if options.savesvg:
    # Some versions of ETE create a "test.svg" and others do not.
    # To avoid confusion (and in case TreeStyle isn't enforced)
    # I just create a new one.
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
