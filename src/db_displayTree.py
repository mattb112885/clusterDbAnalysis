#!/usr/bin/python
# -*- coding: utf-8 -*-

# Display a tree or reroot.
#
# Comes with a decent set of options.

import sys
import optparse
import os
from locateDatabase import *

####################################
# Lets read input arguments first.
####################################
usage = "%prog [options] Newick_file . Default activity is to do nothing - one of -s, -p, -n, or -d must be specified..."
description="Display a tree with annotations and specified root and formats"
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

# Look for the rootorg in the gene names.
if not options.rootorg == None:
    done = False
    for node in t.traverse():
        if options.rootorg in node.name and done:
            sys.stderr.write("WARNING: More than one node present in cluster that matches specified organism %s. Will just use the first one...\n" %(options.rootorg) )
            continue
        if options.rootorg in node.name:
            t.set_outgroup(node)
            done = True
    if not done:
        sys.stderr.write("WARNING No representatives of specified organism %s found - are you using the organism name instead of the ID?\n" %(options.rootorg) )
#        exit(2)

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
    geneToAnnote[spl[0]] = spl[9]
    geneToOrganism[spl[0]] = spl[1]

######################
# Add annotations and
# larger bootstrap values to tree
######################

for node in t.traverse():
    if node.is_leaf():
        # Dont' crash because of e.g. outgroups put in. We already warned about this so don't need to do it again.
        if not (  node.name in geneToOrganism and node.name in geneToAnnote ):
            # We still want to make the text face though so that it is visible!
            F = faces.TextFace(node.name, ftype="Times", fsize=32)
            node.add_face(F, 0, position="aligned")
            continue

        # Add an annotation text with larger font to replace the crappy size-10 ish font that comes by default...
        newname = "_".join( [ node.name, geneToOrganism[node.name], geneToAnnote[node.name] ] )
        F = faces.TextFace(newname, ftype="Times", fsize=32)
        node.add_face(F, 0, position="aligned")
    else:
        # Make the branch support bigger
        F = faces.TextFace(node._support, ftype="Times", fsize=20)
        node.add_face(F, 0, position="branch-top")

# Note - PDF and PNG export are both horrid-quality - I will try to fix this but for now I'll just export to SVG...
# My suggestion is to find something that doesn't crash and that isn't horrible quality to convert this to another format.
# Is there such a thing though??? Sigh...
ts = TreeStyle()
# I already made bigger numbers so no sense in keeping the tiny ones too.
ts.show_branch_support = False
# We'll be putting these in separately
ts.show_leaf_name = False

# Label face columns [This will be useful for labeling tables next to the tree!]
F = faces.TextFace("Annotation", ftype="Times", fsize=20)
ts.aligned_header.add_face(F, 0)

# The default width of the tree is too squished.
# Lets determine a width based on the maximum distance to the root.
maxdist = 0
root = t.get_tree_root()
for node in t.traverse():
    if node.is_leaf():
        dist = t.get_distance(node, root)
        if dist > maxdist:
            maxdist = dist

# Heuristically determine tree width to prevent the tree from getting too squished ...
ts.tree_width = maxdist * 20

# Make ordering the closer for identical trees.
# Direction = 0 - make the outgroup on top.
t.ladderize(direction=0)

# Ladderize doesn't always break ties the same way. Lets fix that, shall we?
# I break ties according to the names of leaves descended from a given node.
# Essentially this amounts to sorting first by number of branches and then by alphabetical order
for node in t.traverse(strategy="levelorder"):
    if not node.is_leaf():
        children = node.get_children()
        print children
        if not len(children) == 2:
            sys.stderr.write("WARNING: Node found with more than two children... Should always have 2 children per node?\n")
            continue
#            exit(2)
        nl0 = len(children[0].get_leaves())
        nl1 = len(children[1].get_leaves())
        if nl0 == nl1:
            names0 = "".join(sorted(children[0].get_leaf_names()))
            names1 = "".join(sorted(children[1].get_leaf_names()))
            if names0 > names1:
                node.swap_children()

if options.savenewick:
    t.write(outfile="%s.nwk" %(options.basename), format=0)

if options.savesvg:
    t.render("%s.svg" %(options.basename), tree_style=ts)

if options.savepng:
    # Convert the svg file into a high-quality (300 dpi) PNG file...
    # The PNG converter in ETE gives a terrible quality image
    # as does the "convert" function (which is probably what ETE uses)
    # so this is the best I could come up with...
    os.system("inkscape -e %s_temp.png -d 300 %s.svg" %(options.basename, options.basename) )
    # Then trim off the edges
    os.system("convert -trim %s_temp.png %s.png" %(options.basename, options.basename))
    os.system("rm %s_temp.png" %(options.basename))

if options.display:
    t.show(tree_style=ts)

con.close()
