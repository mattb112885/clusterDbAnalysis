#!/usr/bin/python

'''Convenient functions for tree-drawing and manipulation with ETE
(common operations used by our scripts)'''

from ete2 import Tree, faces, TreeStyle, NodeStyle, AttrFace, TextFace
import re
import sys
import warnings

def rerootEteTree(ete_tree, root_leaf = None, root_leaf_part = None):
    '''Given an ETE tree, re-root by either the whole name of a leaf (i.e. a gene ID) or a
    part of that name (i.e. an organism ID).

    You are not allowed to re-root by both.
    
    '''

    if root_leaf is None and root_leaf_part is None:
        return ete_tree

    if root_leaf is not None and root_leaf_part is not None:
        raise ValueError("ERROR: Only one of rootgene and rootorg can be specified\n")
    
    if root_leaf is not None:
        done = False
        for node in ete_tree.traverse():
            if root_leaf == node.name:
                ete_tree.set_outgroup(node)
                done = True
                break
        if not done:
            raise ValueError("ERROR: Specified outgroup %s not found in tree\n" %(root_leaf))

    if root_leaf_part is not None:
        done = False
        for node in ete_tree.traverse():
            if root_leaf_part in node.name and done:
                warnings.warn("WARNING: More than one node present in cluster that matches query string %s. Will just use the first one...\n" %(root_leaf_part) )
                continue
            if root_leaf_part in node.name:
                ete_tree.set_outgroup(node)
                done = True
        if not done:
            raise ValueError("ERROR: No leaves found that match specified query string %s found in your tree\n" %(root_leaf_part) )

    return ete_tree

def prettifyTree(ete_tree, leaf_font_size = 32, branch_support_size = 20, title=None, ts = None):
    ''' Perform standardized functions to make the ETE trees easier to read:
    - Make the branch support bigger
    - Make the leaf fonts bigger
    - Change both to standard font (Times)
    - Standardize the tree's width (calculate based on the maximum length from the root to a tip)
    - (optional) add title to tree
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

    #correct the long root node bug (fixed in next release)
    ete_tree.dist=0

    # Optionally create a new TreeStyle if we are passing in an old one.
    if ts is None:
        ts = TreeStyle()

    ts.show_branch_support = False
    ts.show_leaf_name = False

    if title is not None:
        ts.title.clear()
        title = TextFace(title)
        title.hz_align = True
        title.fsize = 52
        ts.title.add_face(title, 0)

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
                warnings.warn("WARNING: Node found with more than two children... Should always have 2 children per node?\n")
                continue
            nl0 = len(children[0].get_leaves())
            nl1 = len(children[1].get_leaves())
            if nl0 == nl1:
                names0 = "".join(sorted(children[0].get_leaf_names()))
                names1 = "".join(sorted(children[1].get_leaf_names()))
                if names0 > names1:
                    node.swap_children()
    
    return ete_tree

# These functions are PhyloTree-specific - PhyloTrees have species attached
# to them while normal ETE trees don't.
def rerootPhyloTree(phylo_tree, reroot_species = None):
    '''Unlike a normal ETE tree object a phyloTree object has a specific "species"
    attached to it. This function identifies the node required by a specific species
    and then calls the more generic reroot command with that node name.

    Note - getting an AttributeError here means you probably passed in the wrong type
    of tree.'''

    if reroot_species is None:
        raise ValueError("ERROR: must specify a species on which to root\n")

    done = False
    for node in phylo_tree.traverse():
        if node.is_leaf():
            if node.species == reroot_species:
                if done:
                    warnings.warn("WARNING: Multiple leaves found matching species %s so just rerooted on the first one we found\n" %(reroot_species) )
                    continue
                phylo_tree = rerootEteTree(phylo_tree, root_leaf = node.name)
                done = True

    if not done:
        sys.stderr.write("ERROR: Specified species %s not found in tree\n" %(reroot_species) )
        raise ValueError

    return phylo_tree

def unsanitizeGeneId(sanitized_geneid):
    '''Turns fig_\d+_\d+_peg_\d+ into fig|\d+.\d+.peg.\d+ for db recognition purposes'''
    return re.sub(r"fig_(\d+)_(\d+)_peg_(\d+)", r"fig|\1.\2.peg.\3", sanitized_geneid)

# TODO - this should probably go into a different library. I put it here for now since James's
# tree function is the only one using it at the moment.
def splitrast(geneid, removefigpeg = False):
    '''Takes an UN-SANITIZED geneid and splits off the organiosm and gene, optionally removing 
    the "fig" and "peg" parts'''

    if ".peg." not in geneid:
        raise ValueError("ERROR: The expected string .peg. not found in gene ID - did you pass in a sanitized ID instead of an unsanitized one?\n")

    unsanitized = geneid
    fig, peg = unsanitized.split('.peg.')
    if removefigpeg:
        fig=fig.lstrip('fig|')
    else:
        peg = 'peg.'+peg

    return fig, peg


def parse_sp_name(node_name):
    '''Parse a node name into an organism ID.
    It is assumed that the node name is in the format "fig|\d+\.\d+\.peg\.\d+".

    We intend for PhyloTrees based on the gene IDs in our database to be constructed using

    pt = phyloTree(sp_naming_function = parse_sp_name)

    The ETE package will then automatically apply this to get species names and add the species name
    to the ".species" field of every node.
    '''
    node_name = str(node_name)
    if node_name=='NoName':
        pass

    unsanitized = unsanitizeGeneId(node_name)
    try:
        orgname, peg = splitrast(unsanitized, removefigpeg = False)
    except ValueError:
        orgname = node_name

    return orgname

    
