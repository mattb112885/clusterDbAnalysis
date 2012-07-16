#!/usr/bin/python

# Given a tree, generates a list of all leaves beneath each node of the tree (i.e. each combination of things
# for which we might want to generate a core list)

from ete2 import Tree
import optparse, os, sys

usage="%prog [options] newick_file > leaflist"
description="""Given a newick file as input, generate a list of all bunches of leaves
in specified format.
Default format is tab-delimited - each bunch gets its own row, and each leaf has a column.
Can also specify to make a separate file for each bunch with one row for each leaf."""
parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-f", "--files", help="Create a separate file for each group in the specified folder (D: False)", action="store", type="str", dest="filedir", default=None)
(options, args) = parser.parse_args()

if len(args) < 1:
    sys.stderr.write("ERROR: Newick file is required\n")
    exit(2)

t = Tree(args[0])

n = 0
# The strategy really doesn't matter, it's just for aesthetics...
for node in t.traverse("postorder"):
    n = n + 1
    if options.filedir is not None:
        p = os.path.join(options.filedir, "%d" %(n) )
        fid = open(p, "w")
        fid.write( "\n".join(node.get_leaf_names()))
        fid.write("\n")
        fid.close()
    else:
        leaflist = "\t".join(node.get_leaf_names())
        print leaflist
    
        
