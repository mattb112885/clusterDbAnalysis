#!/usr/bin/env python

import re
import os
import sys
import sqlite3
import optparse
from FileLocator import *

usage="%prog [options] > presence_absence_table"
description="""Generates a presence - absence table (or slices thereof) based on
the one automatically loaded as part of setup_step2.sh. 
Default activity is to dump the database as is (pegs)."""
parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-n", "--number", help="Rather than printing PEGs, print the number of representatives in each organism (D: Prints pegs)", 
                  action="store_true", dest="number", default=False)
parser.add_option("-b", "--binary", help="Rather than printing PEGs, print 0 if there are no representatives and 1 if there are representatives (D: prints pegs)",
                  action = "store_true", dest="binary", default=False)
parser.add_option("-r", "--runid", help="Only print results for the specified run ID (D: Print the results for all of them)", action="store", type="str", dest="runid", default=None)
parser.add_option("-c", "--clusterid", help="Only print results for the specified cluster ID (D: Prints the table for all of them)", action="store", type="int", dest="clusterid", default=None)
parser.add_option("-t", "--treeorder", help="Given a newick file with the SAME organism names as the presence\absence table, orders the columns to conform with the tree (D: no ordering)",
                  action = "store", type="str", dest="treeorder", default=None)
parser.add_option("-u", "--useronly", help="""WARNING: HACKY and Requires setup_step5 to have been run. 
Only return those genes that were provided by the user (it is assumed that only ITEP genes match fig|\d+\.\d+\.peg\.\d+).""",
                  action = "store_true", dest="useronly", default=False)
parser.add_option("-i", "--iteponly", help="""WARNING: HACKY and Requires setup_step5 to have been run. 
Only return genes originating from ITEP (it is assumed that only ITEP genes match fig|\d+\.\d+\.peg\.\d+).""",
                  action = "store_true", dest="iteponly", default=False)
(options,args) = parser.parse_args()

def treeorder(treefile):
    from ete2 import Tree, faces, TreeStyle, NodeStyle, AttrFace
    t = Tree(treefile)
    rt = t.get_tree_root()
    nameorder = []
    for desc in rt.iter_descendants("preorder"):
        if not desc.is_leaf():
            continue
        nameorder.append(desc.name)
    return nameorder

if options.number and options.binary:
    sys.stderr.write("ERROR: Cannot specify both -n and -b (can either print 0\1 or number of represenatitives, not both)\n")
    exit(2)

if options.number and options.binary:
    sys.stderr.write("ERROR: Cannot specify both -r and -c (supply rin and clusterid to -c )\n")
    exit(2)

if options.iteponly and options.useronly:
    sys.stderr.write("ERROR: Cannot ask for both only itep and only user-specified genes\n")
    exit(2)

con = sqlite3.connect(locateDatabase())
cur = con.cursor()

if options.runid is None and options.clusterid is None:
    cur.execute("SELECT * FROM presenceabsence;")
else:
    if options.clusterid is None:
        cur.execute("SELECT * FROM presenceabsence WHERE runid = ?", (options.runid,))
    else: 
        runid = options.runid
        clustid = options.clusterid
        cur.execute("SELECT * FROM presenceabsence WHERE runid = ? AND clusterid = ?", (runid, clustid))

nameorder = []
if options.treeorder is not None:
    nameorder = treeorder(options.treeorder)

# Get a list of organism names. If we want to sort them in tree-order,
# We need to check that all of the names are consistent with what is in the
# tree.
collist = [tup[0] for tup in cur.description]
newcol2dbcol = {}
if options.treeorder is not None:
    for ii in range(len(collist)):
        # Mapping from the column order we want to the column order we have.
        if collist[ii] not in nameorder:
            # Don't bother giving me a warning if it's just a label for the run ID, cluster ID or annotation (first three columns)
            if ii >= 3:
                sys.stderr.write("WARNING: Organism name %s in the database was not found in the provided tree. It will be deleted!!\n" %(collist[ii]))
        else:
            idx = nameorder.index(collist[ii])
            newcol2dbcol[idx] = ii

if options.treeorder is None:
    print "\t".join(collist)
else:
    newcollist = collist[0:3]
    for ii in range(len(nameorder)):
        if ii in newcol2dbcol:
            newcollist.append(collist[newcol2dbcol[ii]])
    print "\t".join(newcollist)

istitle = True
titles = []
itep_finder = re.compile("fig\|\d+\.\d+\.peg\.\d+")
for rec in cur:
    lst = [ str(s) for s in rec ]
    # HACK ALERT!
    # We pull things out of the list that do not match our expected ITEP format
    if options.useronly:
        for ii in range(len(lst)):
            # Skip info columns and don't replace NONEs
            if ii < 3 or lst[ii] == "NONE":
                continue
            spl = lst[ii].split(";")
            keep = []
            for s in spl:
                if itep_finder.match(s) is None:
                    keep.append(s)
                pass
            if len(keep) == 0:
                lst[ii] = "NONE"
            else:
                lst[ii] = ";".join(keep)
                pass
            pass
        pass
    # Might as well add this too... pull things out of the list that DO match our ITEP format
    elif options.iteponly:
        for ii in range(len(lst)):
            if ii < 3 or lst[ii] == "NONE":
                continue
            spl = lst[ii].split(";")
            keep = []
            for s in spl:
                if itep_finder.match(s) is not None:
                    keep.append(s)
                pass
            if len(keep) == 0:
                lst[ii] = "NONE"
            else:
                lst[ii] = ";".join(keep)
                pass
            pass
        pass

    for ii in range(len(lst)):
        if ii < 3:
            continue
        if (options.number or options.binary) and lst[ii] == "NONE":
            lst[ii] = "0"
        elif options.binary:
            lst[ii] = "1"
        elif options.number:
            lst[ii] = str(lst[ii].count(";") + 1)

    # If we want to re-order, do it now.
    if options.treeorder is not None:
        newlst = lst[0:3]
        for ii in range(len(nameorder)):
            if ii in newcol2dbcol:
                newlst.append(lst[newcol2dbcol[ii]])
        lst = newlst
    print "\t".join(lst)

con.close()
