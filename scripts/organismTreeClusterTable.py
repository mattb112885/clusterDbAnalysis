#!/usr/bin/env python

'''
Draw presence / absence trees from:
    db_getPresenceAbsenceTable_pegs.py > PresenceAbsenceTable_pegs.tsv 
    a list of genes of interest
    *_PresenceAbsenceTable_pegs.tsv with the pegs for each gene
    and maps these onto a organismal tree
'''
#to use Matt's iTEP scripts
import sys, os
sys.path.append("/home/jamesrh/.gvfs/SFTP on biome.life.uiuc.edu/data/Cluster_Files/src") 
from locateDatabase import *

#load Matt's imports from db_displayTree.py
#from ete2 import Tree, faces, TreeStyle, TextFace, NodeStyle, AttrFace
#import fileinput
#import sqlite3
#import sys
#import optparse
#import os
#from locateDatabase import *
#from sanitizeString import *

import itertools, colorsys, math
from sanitizeString import *

import numpy as np
import matplotlib.mlab as mlab
#import matplotlib.pyplot as plt


#aptitude install python-setuptools python-numpy python-qt4 python-scipy python-mysqldb python-lxml
#easy_install -U ete2
from ete2 import Tree, faces, TreeStyle, TextFace, NodeStyle

#lookup by organism abbrv
def distlookup(node1, node2):
    return TreeNode.get_distance(t&node1,t&node2)

def colormap(valuelist):
    #generate list of divergent colors
    clusters = np.unique(valuelist)
    N = len(clusters)
    perm = int(math.ceil(math.sqrt(N)))
    H = [x*1.0/perm for x in range(perm)]
    S = [x*1.0/perm for x in range(perm)]
    HSV = list(itertools.product(H, S, [0.5]))
    RGB = [colorsys.hsv_to_rgb(h,s,v) for h, s, v in HSV]
    n = lambda x: int(x*255)
    RGB256 = [(n(r),n(g),n(b)) for r,g,b in RGB]
    colors = ['#%02x%02x%02x' % (r, g, b) for r, g, b in RGB256]
    colorlookup = dict(zip(clusters, colors[:N]))
    return colorlookup


def drawtree(gene, table):
    #get the runs in this table
    RunIDs = np.unique(table['runid']).tolist()
    #made by: cat /data/Cluster_Files/results/CORE_all_methanogens_09_25_2012/cat_alignment_tree.nwk | replaceOrgWithAbbrev.py -a -f /data/Cluster_Files/organisms > cat_alignment_tree_abbrv.nwk
    #load tree (org abbrv)
    t = Tree("cat_alignment_tree_rep.nwk")
    #root species tree
    #use Methanopyrus kandleri AV19
    for node in t.traverse():
        if node.name == "Methanopyrus_kandleri_AV19":
            t.set_outgroup(node)
            #correct the long root node bug (finxed in next release)
            t.dist=0

    #from Matt's iTol db_displayTree.py
    for node in t.traverse():
        if node.is_leaf():
            # Add an annotation text with larger font to replace the crappy size-10 ish font that comes by default...
            #newname = "_".join( [ node.name, geneToOrganism[node.name], geneToAnnote[node.name] ] )
            F = faces.TextFace(node.name, ftype="Times", fsize=18)
            node.add_face(F, 0, position="aligned")
        else:
            # Make the branch support bigger
            F = faces.TextFace(node._support, ftype="Times", fsize=20, fgcolor = 'red')
            node.add_face(F, 0, position="branch-top")

    t.ladderize(direction=0)

    # Ladderize doesn't always break ties the same way. Lets fix that, shall we?
    # I break ties according to the names of leaves descended from a given node.
    # Essentially this amounts to sorting first by number of branches and then by alphabetical order
    for node in t.traverse(strategy="levelorder"):
        if not node.is_leaf():
            children = node.get_children()
            print children
            if not len(children) == 2:
                sys.stderr.write("WARNING: Node found with more than two children... Should always have 2 children per node?\n")
                continue
    #            exit(2)
            nl0 = len(children[0].get_leaves())
            nl1 = len(children[1].get_leaves())
            if nl0 == nl1:
                names0 = "".join(sorted(children[0].get_leaf_names()))
                names1 = "".join(sorted(children[1].get_leaf_names()))
                if names0 > names1:
                    node.swap_children()

    ts = TreeStyle()
    #some formatting
    margins = 20
    ts.branch_vertical_margin=margins
    ts.show_branch_length = True
    ts.show_branch_support=True
    ts.margin_left = margins
    ts.margin_right = margins
    ts.margin_top = margins
    ts.margin_bottom = margins
    ts.show_branch_support = False
    ts.show_leaf_name = False
    #overall title
    title = TextFace(gene)
    title.hz_align = True
    title.fsize = 52
    ts.title.add_face(TextFace(' '), 0)
    ts.title.add_face(TextFace(' '), 1)
    ts.title.add_face(title, 2)

    coltitle = TextFace("all pegs")
    coltitle.hz_align = True
    coltitle.fsize = 18
    coltitle.fstyle = 'Bold'
    ts.aligned_header.add_face(coltitle, column=len(RunIDs)*2+2)#move over to far left

    # The default width of the tree is too squished.
    # Lets determine a width based on the maximum distance to the root.
    maxdist = 0
    root = t.get_tree_root()
    for node in t.traverse():
        if node.is_leaf():
            dist = t.get_distance(node, root)
            if dist > maxdist:
                maxdist = dist

    # Heuristically determine tree width to prevent the tree from getting too squished ...
    # I haven't figured out why but this doesn't work on some servers.
    # For now I just catch the error and move on...
    try:
        ts.tree_width = maxdist * 20
    except ValueError:
        sys.stderr.write("WARNING: Your ETE setup does not appear to support tree_width - the resulting tree might look strange!\n")

    #see info in: http://packages.python.org/ete2/reference/reference_treeview.html?highlight=treestyle#ete2.TreeStyle
    #put header columns to label runs
    for n, RunID in enumerate(RunIDs):
        ClusterIDs = table[table['runid']==RunID]['clusterid']#can be more than one
        coltitle = TextFace("      " + RunID + "  ")
        coltitle.hz_align = True
        coltitle.fsize = 18
        coltitle.fstyle = 'Bold'
        ts.aligned_header.add_face(coltitle, column=n*2+1)#move over for organism name in col 0

        #put legend with colors at bottom to display cluster IDs.  if multiple lines are from the same run, there will be more than one.  
        for Cluster in ClusterIDs:
            color = getcolor[Cluster]
            text = TextFace("Cluster ID " + str(Cluster))
            text.hz_align = True
            text.fsize = 18
            text.fstyle = 'Bold'
            text.background.color = color
            ts.aligned_foot.add_face(text, column=n*2+1) #move over for organism name in col 0

    #make row for each leaf in tree
    for leaf in t.iter_leaves():
        pegsseen = []
        leaf.set_style(NodeStyle())
        name = leaf.name
        #numclust = 0
        #numpeg = 0
        for n, pegs in enumerate(table[name]):
            if pegs == 'NONE': continue
            #numclust =+1
            #ClusterID = table['ClusterID'][n]
            #colname = rename_RunID(RunID)
            #index which collumn on the table each run will be (as they are ordred)
            colnum = RunIDs.index(table['runid'][n])*2 +1 #to make space for the organism name
            peglist = pegs.split(";")
            peglist.sort()
            for m, peg in enumerate(peglist):
                #based on info from: http://packages.python.org/ete2/tutorial/tutorial_drawting.html#node-faces
                color = getcolor[table['clusterid'][n]]
                text = TextFace(peg)
                text.background.color = color
                leaf.add_face(text, colnum, "aligned") #to put pegs in table
                #numpeg += 1
                pegsseen.append(peg)
            #TODO: REVERT put numper of pegs in this run in this organism
            #text = TextFace(str(len(peglist)) + "   ")
            text = TextFace("   ")
            #text.background.color = color
            text.hz_align = True
            leaf.add_face(text, colnum + 1, "aligned") 
            #numpeg += 1
        numpegs = len(set(pegsseen))
        text = TextFace(numpegs)
        leaf.add_face(text, len(RunIDs)*2+2, "aligned") #all pegs seen
                

    #t.show(tree_style=ts)
    os.system("rm test.svg 2> /dev/null")
    t.render("%s_cluster_tree.svg" %gene, tree_style=ts)
    os.system("convert -trim -depth 32 -background transparent %s_cluster_tree.svg %s_cluster_tree.png" %(gene, gene))


if __name__=="__main__":
    #open tsv, should be tabbed
    #get headers from overall file
    #all = np.genfromtxt("PresenceAbsenceTable_pegs.tsv", names=True)
    headers =  open("PresenceAbsenceTable_pegs.tsv","r").next().strip().split("\t")
    #RunID, ClusterID, SampleAnnote = headers[0:3]
    organisms = headers[3:]
    #gene genes we are interested in 
    geneslist = mlab.csv2rec("genes_of_interest.txt", names="genes")['genes']
    files = [gene + "_PresenceAbsenceTable_pegs.tsv" for gene in geneslist]
    PresenceAbsenceTableFiles={}
    PresenceAbsenceTableFiles.update(zip(geneslist, files))
    for gene, tablefile in PresenceAbsenceTableFiles.iteritems():
        try: table = mlab.csv2rec(tablefile, delimiter="\t", names = headers)#, missing='NONE ')
        except ValueError: next
        table.dtype.names = [sanitizeString(sp, False) for sp in table.dtype.names]
        clusters = np.unique(table['clusterid'])
        getcolor = colormap(clusters)
        drawtree(gene, table)

