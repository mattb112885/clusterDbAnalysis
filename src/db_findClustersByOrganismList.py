#!/usr/bin/python

# Given a list of organisms to stdin and a cluster run,
# identify clusters within that run that only have representatives
# within the specified list of organisms

import fileinput, optparse, sqlite3, sys, os
from locateDatabase import *
from sanitizeString import *

usage="%prog [options] run_id < organism_list > cluster_run_id_list"
description="""Find clusters with a paritcular quality relative to the list of organisms you specified.
Note: To find core gene clusters for a particular group, use both -a and -u
To find core genes only in a parituclar group (to the exclusion of all the others in that cluster run), use -a, -u, and -s
By some group's definitions core genes can be duplicates in some genomes. In such a case exclude the -u.
To find clusers that exclude all the specified organisms use -n
Using only -u will result in an error.
"""
parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-a", "--all", help="Only include clusters that have at least one representative in ALL of the specified organisms", action="store_true", dest="all", default=False)
parser.add_option("-y", "--any", help="Only include clusters that have at least one representative in AT LEAST ONE of the specified organisms", action="store_true", dest="any", default=False)
parser.add_option("-s", "--only", help="Only include clusters that ONLY has matches in the specified organisms", action="store_true", dest="only", default=False)
parser.add_option("-n", "--none", help="Only include clusters that have DOES NOT have a representative in the specified organisms", action="store_true", dest="none", default=False)
parser.add_option("-u", "--uniq", help="Only include clusters that contain exactly ONE representative in any matching organisms (D: Any number)", action="store_true", dest="uniq", default=False)
parser.add_option("-o", "--orgcol", help="Column number for organism starting from 1 (D=1)", action="store", type="int", dest="oc", default=1)
parser.add_option("-r", "--sanitized", help="If specified, the names in the input file have been sanitized (with sanitizeString.py) (D: False)", action="store_true", dest="sanitized", default=False)
#parser.add_option("-m", "--minorgs", help="Minimum number of organisms in clusters to be included (D=no minimum)", action="store", type="int", dest="minorg", default=None)
(options, args) = parser.parse_args()

if not len(args) == 1:
    sys.stderr.write("ERROR: Run ID must be provided\n")
    exit(2)

if options.all and options.none:
    sys.stderr.write("ERROR: ALL and NONE options are contradictory\n")
    exit(2)

if options.any and options.none:
    sys.stderr.write("ERROR: ANY and NONE options are contradictory\n")
    exit(2)

if options.only and options.none:
    sys.stderr.write("ERROR: ONLY and NONE options are contradictory\n")
    exit(2)

if not (options.only or options.all or options.any or options.none):
    sys.stderr.write("ERROR: At least one of -a, -y, -s, -n must be specified\n")
    exit(2)

oc = options.oc - 1

orglist = []
for line in fileinput.input("-"):
    spl = line.strip("\r\n").split("\t")
    orglist.append(spl[oc])

if options.sanitized:
    allOrgsDict = {}
    p = os.path.join(os.path.dirname(locateDatabase()), "..", "organisms")
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
cur.execute("SELECT runid, clusterid, organism FROM clusterorgs WHERE clusterorgs.runid=? ORDER BY runid,clusterid", (args[0], ) )
for res in cur:
    ls = [ str(s) for s in res ]
    cl.append(ls)

previd = -1
orgset = set(orglist)
currentorgs = set()
for l in cl:
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
            if not ( ( options.any and not anyok) or (options.all and not allok) or (options.none and not noneok) or (options.only and not onlyok) or (options.uniq and not uniqok) ):       
                print "%s\t%s" %(prevrun, previd)

        # Reset
        uniqok = True
        currentorgs.clear()
        previd = l[1]
        prevrun = l[0]

    if l[2] in currentorgs:
        uniqok = False
    currentorgs.add(l[2])


con.close()
