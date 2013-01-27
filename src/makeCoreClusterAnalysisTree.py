#!/usr/bin/python

from ete2 import Tree, TextFace, CircleFace
import optparse, sys
from CoreGeneFunctions import *
from TreeFuncs import *

usage="%prog [options] newick_file runid"
description="""Generate an ETE tree with internal node labels corresponding to the number of
clusters conserved in the nodes beneath it (conservation being defined by a variery of options
below). The input MUST be a Newick file with organism names replaced."""

parser = optparse.OptionParser(usage=usage, description=description)
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

(options, args) = parser.parse_args()

if len(args) < 2:
    sys.stderr.write("ERROR: Newick file and runID are required\n")
    exit(2)

t = Tree(args[0])
runid = args[1]

if options.reroot_org is not None:
    t = rerootEteTree(t, root_leaf = options.reroot_org)

# Make something pretty out of it.
t, ts = prettifyTree(t)

# The strategy really doesn't matter, it's just for aesthetics...
for node in t.traverse("postorder"):
    # This function from ETE gives you all of the descendent leaf names in an array.
    leafnames = node.get_leaf_names()
    clusters = findGenesByOrganismList(leafnames, runid, sanitized=True, all_org = options.all, any_org = options.any,
                                       only_org = options.only, none_org = options.none, uniq_org = options.uniq)
    numclusters = len(clusters)
    print "%d clusters for organisms:" %(numclusters)
    print leafnames
    print ""
    numFace = TextFace(str(numclusters), ftype="Times", fsize=24)
    # Numclusters will be as high as a few thousand for the bacteria (with -all) - with -any it can be way more than this so
    # use with caution...
    cFace = CircleFace(radius=int(numclusters/100)+10, color="Green", style="sphere")
    node.add_face(cFace, 0, position="float")
    node.add_face(numFace, 0, position="branch-right")
    

t.show(tree_style = ts)
