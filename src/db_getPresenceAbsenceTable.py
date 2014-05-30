#!/usr/bin/env python

import re
import os
import sys
import sqlite3
import optparse
from ClusterFuncs import *
from FileLocator import *
from sanitizeString import *

usage="""%prog [options] > presence_absence_table

Output: Run_id Cluster_id Representative_annotation (organism1_PA) (organism2_PA) ... """
description="""Generates a presence - absence table (or slices thereof) based on pre-computed clusters.
Default activity is to print all available presence-absence from all clusters and all cluster runs."""
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

# If a run ID is specified, we want to only return columns in that run.
orgsToInclude = None
if options.runid is not None:
    orgsToInclude = set()
    orgs = getOrganismsInClusterRun(options.runid, cur)
    if len(orgs) == 0:
        raise IOError("ERROR: Specified run ID %s does not exist in the database." %(options.runid) )        
    for org in orgs:
        orgsToInclude.add(sanitizeString(org, False))

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

# Get a list of organism names. 
collist = [tup[0] for tup in cur.description]

# Mapping from the column order we want to the column order we have.
newcol2dbcol = {}
for ii in range(len(collist)):
    if orgsToInclude is not None:
        if collist[ii] not in orgsToInclude:
            if ii >= 3:
                # Don't bother giving me a warning if it's just a label for the run ID, cluster ID or annotation (first three columns)
#                sys.stderr.write("WARNING: Organism name %s in the database was not found in the requested cluster run. It will be deleted!!\n" %(collist[ii]))
                continue
    if options.treeorder is None:
        # Keep the same ordering that exists
        newcol2dbcol[ii] = ii
    else:
        if collist[ii] not in nameorder:
            # Don't bother giving me a warning if it's just a label for the run ID, cluster ID or annotation (first three columns)
            if ii >= 3:
                sys.stderr.write("WARNING: Organism name %s in the database was not found in the provided tree. It will be deleted!!\n" %(collist[ii]))
                continue
        else:
            idx = nameorder.index(collist[ii]) + 3
            newcol2dbcol[idx] = ii

# Reorder columns in header
newcollist = collist[0:3]
for ii in range(len(collist)):
    if ii >= 3 and ii in newcol2dbcol:
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

    # Convert format if requested
    for ii in range(len(lst)):
        if ii < 3:
            continue
        if (options.number or options.binary) and lst[ii] == "NONE":
            lst[ii] = "0"
        elif options.binary:
            lst[ii] = "1"
        elif options.number:
            lst[ii] = str(lst[ii].count(";") + 1)

    # Reorder rows (and throw out rows that don't match organism lists we need - either from tree or from 
    # run ID)
    newlst = lst[0:3]
    for ii in range(len(lst)):
        if ii >= 3 and ii in newcol2dbcol:
            newlst.append(lst[newcol2dbcol[ii]])
    lst = newlst
    print "\t".join(lst)

con.close()
