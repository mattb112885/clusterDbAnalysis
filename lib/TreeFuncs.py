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


def prettifyTree(ete_tree, leaf_font_size = 32, branch_support_size = 20):
    ''' Perform standardized functions to make the ETE trees easier to read:
    - Make the branch support bigger
    - Make the leaf fonts bigger
    - Change both to standard font (Times)
    - Standardize the tree's width (calculate based on the maximum length from the root to a tip)
    '''

    for node in ete_tree.traverse():
        if node.is_leaf():
            # Make text faces with name = the existing node name but with big font.
            # A side effect of this is that we also get the annotations lined up
            F = faces.TextFace(node.name, ftype="Times", fsize=leaf_font_size)
            node.add_face(F, 0, position="aligned")
        else:
            # Make branch support bigger
            F = faces.TextFace(node._support, ftype="Times", fsize=branch_support_size)
            node.add_face(F, 0, position="branch-top")
    
    ts = TreeStyle()
    ts.show_branch_support = False
    ts.show_leaf_name = False

    # The default width of the tree is too squished.
    # Lets determine a width based on the maximum distance to the root.
    maxdist = 0
    root = ete_tree.get_tree_root()
    for node in ete_tree.traverse():
        if node.is_leaf():
            dist = ete_tree.get_distance(node, root)
            if dist > maxdist:
                maxdist = dist

    # I haven't figured out why but this doesn't work on some servers.
    # For now I just catch the error and move on... 
    try:
        ts.tree_width = maxdist * 20
    except ValueError:
        sys.stderr.write("WARNING: Your ETE setup does not appear to support tree_width - the resulting tree might look strange!\n")

    return ete_tree, ts

def standardizeTreeOrdering(ete_tree):
    '''Make ordering of the leaves the same for identical trees.'''

    # Direction = 0 - make the outgroup on top.
    ete_tree.ladderize(direction=0)

    # Ladderize doesn't always break ties the same way. Lets fix that, shall we?
    # I break ties according to the names of leaves descended from a given node.
    # Essentially this amounts to sorting first by number of branches and then by alphabetical order 
    for node in ete_tree.traverse(strategy="levelorder"):
        if not node.is_leaf():
            children = node.get_children()
            if not len(children) == 2:
                sys.stderr.write("WARNING: Node found with more than two children... Should always have 2 children per node?\n")
                continue
            nl0 = len(children[0].get_leaves())
            nl1 = len(children[1].get_leaves())
            if nl0 == nl1:
                names0 = "".join(sorted(children[0].get_leaf_names()))
                names1 = "".join(sorted(children[1].get_leaf_names()))
                if names0 > names1:
                    node.swap_children()
    
    return ete_tree
