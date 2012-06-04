#!/usr/bin/python
# -*- coding: utf-8 -*-

# Display a tree or reroot.
#
# Comes with a decent set of options.

import sys
import optparse

####################################
# Lets read input arguments first.
####################################
usage = "%prog [options] Newick_file . Default activity is to do nothing - one of -s, -p, -n, or -d must be specified..."
parser = optparse.OptionParser(usage=usage)
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
        sys.stderr.write("ERROR: No representatives of specified organism %s found - are you using the organism name instead of the ID?\n" %(options.rootorg) )
        exit(2)

##############################################
# Get various gene / organism / annotation / neighborhood / cluster info out of the database
##############################################

# Annotation and organism
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

#############
# Add actual arrow objects to the tree
#############
sys.stderr.write("Adding arrow objects to leaves of the tree...\n")
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
# The default width of the tree is too squished.
ts.tree_width = 1000

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
