#!/usr/bin/python

import fileinput
import optparse
import sys

from ete2 import Tree

usage = "%prog < tree.nwk > list_of_leaf_names"
description = '''Given a Newick tree, 
returns a list of leaf names in the tree
'''

parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-f", "--file", help="Read single Newick tree from specified file (default: -, stdin)", action="store", type="str", dest="infile", default="-")
(options, args) = parser.parse_args()

ii = 0
for line in fileinput.input(options.infile):
    if ii > 0:
        sys.stderr.write("ERROR: Should only specify a single Newick tree in the input file.\n")
        exit(2)
    t = Tree(line)
    for node in t.traverse(strategy="postorder"):
        if node.is_leaf():
            print node.name
    ii += 1
