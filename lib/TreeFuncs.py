#!/usr/bin/python

'''Convenient functions for tree-drawing and manipulation with ETE
(common operations used by our scripts)'''

from ete2 import Tree, faces, TreeStyle, NodeStyle, AttrFace
import re
import sys

def rerootEteTree(ete_tree, root_leaf = None, root_leaf_part = None):
    '''Given an ETE tree, re-root by either the whole name of a leaf (i.e. a gene ID) or a
    part of that name (i.e. an organism ID).

    You are not allowed to re-root by both.
    
    '''

    if root_leaf is None and root_leaf_part is None:
        return ete_tree

    if root_leaf is not None and root_leaf_part is not None:
        sys.stderr.write("ERROR: Only one of rootgene and rootorg can be specified\n")
    
    if root_leaf is not None:
        done = False
        for node in ete_tree.traverse():
            if root_leaf == node.name:
                ete_tree.set_outgroup(node)
                done = True
                break
        if not done:
            sys.stderr.write("ERROR: Specified outgroup %s not found in tree\n" %(root_leaf))
            raise ValueError

    if root_leaf_part is not None:
        done = False
        for node in ete_tree.traverse():
            if root_leaf_part in node.name and done:
                sys.stderr.write("WARNING: More than one node present in cluster that matches query string %s. Will just use the first one...\n" %(root_leaf_part) )
                continue
            if root_leaf_part in node.name:
                ete_tree.set_outgroup(node)
                done = True
        if not done:
            sys.stderr.write("ERROR: No leaves found that match specified query string %s found in your tree\n" %(root_leaf_part) )
            raise ValueError

    return ete_tree
