#!/usr/bin/env python

from __future__ import print_function
import fileinput
import optparse
import sys

from ete2 import Tree

from TreeFuncs import *

usage = """ %prog -l root_leaf_name < Newick_tree > Newick_tree_rerooted
%prog -p root_leaf_part < Newick_tree > Newick_tree_reroorted """
description = """ Given a Newick file and either the name or part of the
name of a leaf on which to reroot the tree, reroots the tree and returns a
rerooted Newick file (does no further processing, unlike other scripts like
db_displayTree). If given a part of the name only, the match must be
unique. If given the entire name it must be present exactly as written in the
tree."""

parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-l", "--leafname", help="Leaf name on which to reroot (must exactly match the text in a leaf name)", action="store", type="str", dest="leafname", default=None)
parser.add_option("-p", "--leafpart", help="PART of a leaf name on which to reroot. The string must be part of exactly one of the leaf names.", action="store", type="str", dest="leafpart", default=None)
(options, args) = parser.parse_args()

if options.leafname is None and options.leafpart is None:
    sys.stderr.write("ERROR: Either leafname (-l) or part of a leafname (-p) on which to reroot is required.\n")
    exit(2)

if (options.leafname is not None and options.leafpart is not None):
    sys.stderr.write("ERROR: Can only specify exactly one of -l or -p (leafname or leaf part)\n")
    exit(2)

newick_str = "".join( [ line.strip("\r\n") for line in fileinput.input("-") ] )

ete_tree = Tree(newick_str)

root_leaf = options.leafname
root_leaf_part = options.leafpart
ete_tree = rerootEteTree(ete_tree, root_leaf = root_leaf, root_leaf_part = root_leaf_part)

print(ete_tree.write())
