#!/usr/bin/python

# Given a list of organisms to stdin and a cluster run,
# identify clusters within that run that only have representatives
# within the specified list of organisms

import fileinput, operator, optparse, sqlite3, sys, os
from FileLocator import *
from sanitizeString import *

def findGenesByOrganismList(orglist, runid, sanitized = False, any_org = False, all_org = False, only_org = False, none_org = False, uniq_org = False):
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

    con = sqlite3.connect(locateDatabase())
    cur = con.cursor()

    # From the sqlite database, download the list of clusterorgs
    cl = []
    #cur.execute("SELECT runid, clusterid, organism FROM clusterorgs WHERE clusterorgs.runid=? ORDER BY runid,clusterid", (runid, ) )
    cur.execute("SELECT runid, clusterid, organism FROM clusterorgs WHERE clusterorgs.runid=?", (runid, ) )
    for res in cur:
        ls = [ str(s) for s in res ]
        cl.append(ls)

    # Lets see if this is faster than using ORDER BY in the SQL statement...
    cl = sorted(cl, key=operator.itemgetter(1,2) )

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

    con.close()
    return goodClusters
