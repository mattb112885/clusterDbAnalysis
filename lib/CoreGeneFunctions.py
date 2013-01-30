#!/usr/bin/python

# Given a list of organisms to stdin and a cluster run,
# identify clusters within that run that only have representatives
# within the specified list of organisms

import fileinput, operator, optparse, sqlite3, sys, os
from FileLocator import *
from sanitizeString import *
from ete2 import TextFace

def addCoreDataToTree(ete_tree, runid, sanitized = False, any_org = False, all_org = False, only_org = False, none_org = False, uniq_org = False, color = "Black"):
    '''A function to add data related to gene and organism distribution across clusters
    to a core gene tree.

    See http://packages.python.org/ete2/reference/reference_treeview.html#color-names for a list
    of valid color names.'''

    cl = getClusterOrgsByRun(runid)

    nodenum = 0
    clusterrunlist = []
    # The strategy really doesn't matter, it's just for aesthetics... and to make sure its always the same.
    for node in ete_tree.traverse("postorder"):
        nodenum += 1
        leafnames = node.get_leaf_names()
        clusters = findGenesByOrganismList(leafnames, runid, cl = cl, sanitized=True, 
                                           all_org = all_org, any_org = any_org, only_org = only_org, none_org = none_org, uniq_org = uniq_org)
        numclusters = len(clusters)
        # This is mostly so that I can keep track of progress.
        sys.stderr.write("%d (N%d)\n" %(numclusters, nodenum))
        numFace = TextFace("%d (N%d)" %(numclusters, nodenum), ftype="Times", fsize=24, fgcolor=color)
        node.add_face(numFace, 0, position="branch-bottom")
        for c in clusters:
            clusterrunlist.append( ( c[0], c[1], nodenum ) )

    return ete_tree, clusterrunlist

def getClusterOrgsByRun(runid):
    '''I separated this call from the findGenesByOrganismList below because we often need to call the latter
    many times and this one is the same for all of them (if the run ID is the same).

    The return object is a list of (runid, clusterid, organism) tuples sorted by run ID then by cluster ID.'''
    # From the sqlite database, download the list of clusterorgs

    con = sqlite3.connect(locateDatabase())
    cur = con.cursor()

    cl = []
    cur.execute("SELECT runid, clusterid, organism FROM clusterorgs WHERE clusterorgs.runid=?", (runid, ) )
    for res in cur:
        ls = [ str(s) for s in res ]
        cl.append(ls)

    # This is far faster than using ORDER BY in sqlite (at least unless I try to hack the SQL command to make it
    # get the order of operations right...)
    cl = sorted(cl, key=operator.itemgetter(1,2) )

    con.close()

    return cl

def findGenesByOrganismList(orglist, runid, cl = None, sanitized = False, any_org = False, all_org = False, only_org = False, none_org = False, uniq_org = False):
    '''Identify clusters that have a specific set of properties with respect to a given set of
    organisms. The valid properties are ANY, ALL, ONLY, and NONE.

    Specifiy sanitized=TRUE if the organism names passed here are sanitized (spaces, periods, etc. replaced by
    underscores - see sanitizeString.py for the standard way to sanitize names).

    If the list of runid, clusterid, organismid tuples has already been computed, pass it in via the "cl"
    argument to avoid computing it again. Otherwise, it will be (re)computed within this function.

    You can also use the "cl" argument to restrict analysis to a specific set of (run ID, cluster ID) pairs
    by just passing that subset to the function. If no "cl" is passed then it is assumed you want to compare against
    ALL clusters in a run.

    The organisms in "orglist" are considered the "ingroup" and any organisms in the given cluster run but
    NOT in the orglist are considered the "outgroup". Clusters are pulled out according to the following table
    where the number in the entry corresponds to the number of represented ORGANISMS (NOT GENES) IN THE INGROUP
    (other combinations are possible - this is just a representative set of examples):

      Property  | Ingroup |  Outgroup
    +-----------+---------+-----------
      ALL       |  == N   |    >= 0
    +-----------+---------+-----------
      ANY       |  >= 1   |    >= 0
    +-----------+---------+-----------
      ONLY      |  >= 1   |    == 0
    +-----------+---------+-----------
      NONE      |  == 0   |    >= 1*
    +-----------+---------+-----------
     ALL + ONLY |  == N   |    == 0    - Genes that are only found in the ingroup and that are found in all members of the ingroup
    +-----------+---------+-----------
     ANY + ONLY |  >= 1   |    == 0    - Genes that are found only in the ingroup (but not necessarily in all of its members)
    +-----------+---------+-----------
     ALL + NONE |
     ANY + NONE | Contradictions (raise errors).
     ONLY + NONE|     
    +-----------+---------+-----------    

    *: No clusters have 0 representatives

    N is the number of organisms in the ingroup and O is the number in the outgroup.

    UNIQ specifies that in addition to any other flags, genes in every organism in the ingroup
    must be uniquely represented in the cluster. Some groups definitions of "core genes" are
    satisfied by using AND and UNIQ as constraints.

    The function returns a list of (runid, clusterid) pairs that adhere to the user-specified criteria.

    (TODO - I need to check if it enforces it for
    only the ingroup or for both the ingroup AND the outgroup. We probably want it to only care
    about the ingroup I think).'''


    if all_org and none_org:
        raise ValueError("ERROR: all_org and none_org options are contradictory\n")
    if any_org and none_org:
        raise ValueError("ERROR: any_org and none_org options are contradictory\n")
    if only_org and none_org:
        raise ValueError("ERROR: only_org and none_org options are contradictory\n")
    if not (only_org or all_org or any_org or none_org):
        raise ValueError("ERROR: At least one of any_org, all_org, none_org, or only_org must be specified.\n")

    # Change sanitized gene names to un-sanitized gene names using the organisms file.
    if sanitized:
        allOrgsDict = {}
        p = locateOrganismFile()
        orgfile = open(p, "r")
        for line in orgfile:
            spl = line.strip("\r\n").split("\t")
            allOrgsDict[sanitizeString(spl[0], False)] = spl[0]

        for ii in range(len(orglist)):
            orglist[ii] = allOrgsDict[orglist[ii]]

    if cl is None:
        cl = getClusterOrgsByRun(runid)

    previd = -1
    orgset = set(orglist)
    currentorgs = set()
    goodClusters = []
    for l in cl:
        # Basically we slurp up all cluster,org pairs corresponding to a specific
        # cluster and then once we have all of them we check if they are unique, have all of the organisms of interest,
        # etc...
        if l[1] != previd:
            if previd != -1:
                (anyok, allok, noneok, onlyok) = False,False,False,False
                # Check ANY
                intersection = orgset & currentorgs
                if len(intersection) > 0:
                    anyok = True
                else:
                    noneok = True
                # Check ALL
                if len(intersection) == len(orgset):
                    allok = True
                # Check ONLY
                diff = currentorgs - orgset
                if len(diff) == 0:
                    onlyok = True

                # Our criteria: we can't have any of the options be TRUE and not have the corresponding condition also be true
                if not ( ( any_org and not anyok) or ( all_org and not allok) or (none_org and not noneok) or (only_org and not onlyok) or (uniq_org and not uniqok) ):       
                    goodClusters.append( (prevrun, previd) )

            # Reset
            uniqok = True
            currentorgs.clear()
            previd = l[1]
            prevrun = l[0]

        if l[2] in currentorgs:
            uniqok = False
        currentorgs.add(l[2])

    return goodClusters
