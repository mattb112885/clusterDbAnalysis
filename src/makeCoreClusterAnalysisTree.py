#!/usr/bin/env python

from __future__ import print_function
import optparse
import os
import sqlite3
import sys
import csv

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
parser.add_option("-d", "--display", 
                  help="Display tree", 
                  action="store_true", dest="display", default=False)
parser.add_option("-v", "--savesvg", 
                  help="Save tree as SVG (requires -b)", 
                  action="store_true", dest="savesvg", default=False)
parser.add_option("-p", "--savepng", 
                  help="save tree as PNG (implies -v, requires -b)", 
                  action="store_true", dest="savepng", default=False)
parser.add_option("-x", "--savexls", 
                  help="Output lists of clusters for each tree to XLS (D: Dont). Requires xlwt.", 
                  action="store_true", dest="savexls", default=False)
parser.add_option("-m", "--save-multi-txt", 
                  help="Output lists of clusters for each tree to MANY tab delimited files, names by the node number (D: Don't). ", 
                  action="store_true", dest="save_multi_txt", default=False)
parser.add_option("-t", "--savetxt", 
                  help="Output a lists of nodes and clusters for the tree to ONE tab delimited file (D: Don't). ", 
                  action="store_true", dest="savetxt", default=False)
parser.add_option("-b", "--basename", 
                  help="Output file base name (D: automatically generated)", 
                  action="store", dest="basename", default=None)
parser.add_option("-e", "--no_annotation", 
                  help="Do NOT look up Representative Annotations (D: Look them up)",
                  action="store_true", dest="no_annotations", default=False)
parser.add_option("-r", "--reroot", 
                  help="Reroot tree to specified leaf before doing calculation (D: Use tree as is)",
                  action="store", type="str", dest="reroot_org", default=None)
parser.add_option("-a", "--all", 
                  help="""Count the clusters containing all of the descendent leaves of each internal node. Only specifying this option, the program ignores every leaf that is not a descendent of the internal node.""",
                  action="store_true", dest="all", default=False)
parser.add_option("-u", "--uniq", 
                  help = """Count only those clusters with EXACTLY ONE member from each descendent leaf (non-descendant leaves are ignoerd).""",
                  action= "store_true", dest="uniq", default=False)
parser.add_option("-s", "--only",
                  help = """Count clusters only if every gene is a member of an organism in a descendant node of the tree. (can be combined with -a or -u)""",
                  action = "store_true", dest="only", default=False)
parser.add_option("-n", "--none",
                  help = """Count clusters only if there are NO representatives in any descendent nodes from a given internal node (it doesn't matter which other organisms are represented in the cluster)""",
                  action = "store_true", dest="none", default=False)
parser.add_option("-y", "--any",
                  help = "Count clusters if they have representatives in ANY descendent node from each internal node.""",
                  action = "store_true", dest="any", default = False)
parser.add_option("-c", "--clades",
                  help = "Instead of comparing these statistics to the whole cluster run as an outgroup, restrict the outgroup to the sister clade(s) of the tree.",
                  action = "store_true", dest="clades", default=False)

(options, args) = parser.parse_args()

if len(args) < 2:
    sys.stderr.write("ERROR: Newick file and runID are required\n")
    exit(2)

if not (options.display or options.savesvg or options.savepng or options.savexls or options.savetxt or options.save_multi_txt):
    sys.stderr.write("ERROR: At least one of -v, -p, -x, -d, -t or -m is required\n")
    exit(2)

#more logic will be checked later, but let's make sure that the user is asking for something, as the latter library error messages might be criptic
if not (options.all or options.uniq or options.only or options.none or options.any):
    sys.stderr.write("ERROR: At least one of -a, -u, -s, -n, -y is required\n")
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
                            none_org = options.none, uniq_org = options.uniq, compare_to_adj_clade = options.clades)

#we will cache in a dictionary, as the lookup is slow
RepresentativeAnnotation={}
def get_annotation(clusterrun):
    if options.no_annotations: 
        annotation = ""
    else:
        try: 
            annotation = RepresentativeAnnotation[clusterrun[1]]
        except KeyError: 
            annotation = findRepresentativeAnnotation(clusterrun[0], clusterrun[1], cur)
            RepresentativeAnnotation[clusterrun[1]]=annotation 
    return annotation

if options.savexls:
    import xlwt
    wb = xlwt.Workbook()
    for nodenum in data:
        sheet = wb.add_sheet(str(nodenum))
        rownum = 0
        for clusterrun in data[nodenum]:
            sheet.write(rownum, 0, clusterrun[0])
            sheet.write(rownum, 1, clusterrun[1])
            sheet.write(rownum, 2, get_annotation(clusterrun))
            rownum += 1
    xlsname = "%s.xls" %(options.basename) 
    wb.save(xlsname) 
    print("Output to %s with tabs of nodes, and collumns:\n\tclusterrunID\tclusterID\tRepresentativeAnnotation\n" % xlsname)
    
if options.savetxt:
    tsvfilename="%s_all_nodes.tsv" % (options.basename)
    with open(tsvfilename,"wb") as tsvfile:
        tsv=csv.writer(tsvfile, delimiter='\t', quoting=False)
        for nodenum in data:
            for clusterrun in data[nodenum]:
                tsv.writerow((nodenum,) + clusterrun + (get_annotation(clusterrun),))
        tsvfile.close()
    print("Output to %s with collumns:\n\tnode\tclusterrunID\tclusterID\tRepresentativeAnnotation\n" % tsvfilename)

if options.save_multi_txt:
    for nodenum in data:
        tsvfilename="%s_node_%s.tsv" % (options.basename, nodenum)
        with open(tsvfilename,"wb") as tsvfile:
            tsv=csv.writer(tsvfile, delimiter='\t', quoting=False)
            for clusterrun in data[nodenum]:
                tsv.writerow(clusterrun + (get_annotation(clusterrun),))
        tsvfile.close()
    print("Output to files like %s with collumns:\n\tnode\tclusterrunID\tclusterID\tRepresentativeAnnotation\n" % tsvfilename)

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
