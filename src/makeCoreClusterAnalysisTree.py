#!/usr/bin/python

from ete2 import Tree, TextFace, CircleFace
import optparse, os, sys
from CoreGeneFunctions import *
from TreeFuncs import *

usage="%prog [options] newick_file runid"
description="""Generate an ETE tree with internal node labels corresponding to the number of
clusters conserved in the nodes beneath it (conservation being defined by a variery of options
below). The input MUST be a Newick file with organism names replaced."""

parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-d", "--display", help="Display tree", action="store_true", dest="display", default=False)
parser.add_option("-v", "--savesvg", help="Save tree as SVG (requires -b)", action="store_true", dest="savesvg", default=False)
parser.add_option("-p", "--savepng", help="save tree as PNG (implies -v, requires -b)", action="store_true", dest="savepng", default=False)
parser.add_option("-r", "--reroot", help="Reroot tree to specified leaf before doing calculation (D: Use tree as is)",
                  action="store", type="str", dest="reroot_org", default=None)
parser.add_option("-a", "--all", 
                  help="""Count the clusters containing all of the descendent leaves of each internal node. Only specifying this
option, the program ignores every leaf that is not a descendent of the internal node.""",
                  action="store_true", dest="all", default=False)
parser.add_option("-u", "--uniq",
                  help = """Count only those clusters with EXACTLY ONE member
ofrom each descendent leaf (non-descendant leaves are ignoerd).""",
                  action= "store_true", dest="uniq", default=False)
parser.add_option("-s", "--only",
                  help = """Count clusters only if every gene is a member of an organism in a descendant node of the tree. (can be combined with -a or -u)""",
                  action = "store_true", dest="only", default=False)
parser.add_option("-n", "--none",
                  help = """Count clusters only if there are NO representatives in any descendent nodes from a given internal node (it doesn't matter which other organisms
are represented in the cluster)""",
                  action = "store_true", dest="none", default=False)
parser.add_option("-y", "--any",
                  help = "Count clusters if they have representatives in ANY descendent node from each internal node.""",
                  action = "store_true", dest="any", default = False)
parser.add_option("-b", "--basename", help="Output file base name (D: automatically generated)", action="store", dest="basename", default=None)

(options, args) = parser.parse_args()

if len(args) < 2:
    sys.stderr.write("ERROR: Newick file and runID are required\n")
    exit(2)

if not (options.display or options.savesvg or options.savepng):
    sys.stderr.write("ERROR: At least One of -v, -p, or -d is required\n")
    exit(2)

if options.savepng:
    options.savesvg = True

if options.basename is None:
    truthvalue_str = []
    if options.all:
        truthvalue_str.append("all")
    if options.any:
        truthvalue_str.append("any")
    if options.only:
        truthvalue_str.append("only")
    if options.none:
        truthvalue_str.append("none")
    if options.uniq:
        truthvalue_str.append("uniq")
    options.basename = "%s_%s" %(args[1], "_".join(truthvalue_str))


t = Tree(args[0])
runid = args[1]

if options.reroot_org is not None:
    t = rerootEteTree(t, root_leaf = options.reroot_org)

# Make something pretty out of it.
# We don't want bootstraps here since they just make the tree more cluttered.
t, ts = prettifyTree(t, show_bootstraps = False)

# The strategy really doesn't matter, it's just for aesthetics...
nodenum = 0
fid = open(options.basename, "w")
for node in t.traverse("postorder"):
    nodenum += 1
    # This function from ETE gives you all of the descendent leaf names in an array.
    leafnames = node.get_leaf_names()
    clusters = findGenesByOrganismList(leafnames, runid, sanitized=True, all_org = options.all, any_org = options.any,
                                       only_org = options.only, none_org = options.none, uniq_org = options.uniq)
    numclusters = len(clusters)
    # This is mostly so that I can keep track of progress.
    sys.stderr.write("%d (%d)\n" %(numclusters, nodenum))
#    numFace = TextFace("\n%d" %(numclusters), ftype="Times", fsize=24)
    numFace = TextFace("\%d (%d)" %(numclusters, nodenum), ftype="Times", fsize=24)
    # Numclusters will be as high as a few thousand for the bacteria (with -all) - with -any it can be way more than this so
    # use with caution...
#    cFace = CircleFace(radius=int(numclusters/50), color="Green", style="sphere")
#    node.add_face(cFace, 0, position="float")  
    node.add_face(numFace, 2, position="branch-bottom")
    for c in clusters:
        fid.write("%s\t%s\t%d\n" %(c[0], c[1], nodenum))

fid.close()

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
