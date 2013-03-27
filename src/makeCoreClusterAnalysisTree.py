#!/usr/bin/python

import optparse
import os
import sqlite3
import sys

from ete2 import Tree, TextFace, CircleFace

from CoreGeneFunctions import *
from ClusterFuncs import *
from FileLocator import *
from TreeFuncs import *

usage="%prog [options] newick_file runid"
description="""
Generate an ETE tree with internal node labels corresponding to the number of
clusters conserved in the nodes beneath it (conservation being defined by a variery of options
below). The input MUST be a Newick file with organism IDs REPLACED with their names.

The function alternatively (or in addition) exports a XLS file with sheet names equal to the node number printed
on the tree containing the cluster, runid pair and a representative annotation from each cluster
identified with these properties...
"""

parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-d", "--display", help="Display tree", action="store_true", dest="display", default=False)
parser.add_option("-v", "--savesvg", help="Save tree as SVG (requires -b)", action="store_true", dest="savesvg", default=False)
parser.add_option("-p", "--savepng", help="save tree as PNG (implies -v, requires -b)", action="store_true", dest="savepng", default=False)
parser.add_option("-x", "--savexls", help="Output lists of clusters for each tree to XLS (D: Dont). Requires xlwt.", action="store_true", dest="savexls", default=False)
parser.add_option("-b", "--basename", help="Output file base name (D: automatically generated)", action="store", dest="basename", default=None)

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

if not (options.display or options.savesvg or options.savepng or options.savexls):
    sys.stderr.write("ERROR: At least One of -v, -p, -x or -d is required\n")
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

con = sqlite3.connect(locateDatabase())
cur = con.cursor()

t = Tree(args[0])
runid = args[1]

if options.reroot_org is not None:
    t = rerootEteTree(t, root_leaf_part = options.reroot_org)

# Make something pretty out of it.
# We don't want bootstraps here since they just make the tree more cluttered.
t, ts = prettifyTree(t, show_bootstraps = False)

t, data = addCoreDataToTree(t, runid, 
                            all_org = options.all, any_org = options.any, only_org = options.only, 
                            none_org = options.none, uniq_org = options.uniq)

if options.savexls:
    import xlwt
    wb = xlwt.Workbook()
    for nodenum in data:
        sheet = wb.add_sheet(str(nodenum))
        rownum = 0
        for clusterrun in data[nodenum]:
            sheet.write(rownum, 0, clusterrun[0])
            sheet.write(rownum, 1, clusterrun[1])
            sheet.write(rownum, 2, findRepresentativeAnnotation(clusterrun[0], clusterrun[1], cur))
            rownum += 1
    wb.save("%s.xls" %(options.basename) )        

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
