#!/usr/bin/env python

import fileinput, optparse, sys, sqlite3
from FileLocator import *

usage="""%prog [options] < runid_clusterid_pair > organism_list"""
description="Get a list of organisms included in the specified run / cluster ID pair"
parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-r", "--runcol", help="Column number for run ID starting from 1 (D=1)", action="store", type="int", dest="rc", default=1)
parser.add_option("-c", "--clustercol", help="Column number for cluster ID starting from 1 (D=2)", action="store", type="int", dest="cc", default=2)
(options, args) = parser.parse_args()

rc = options.rc - 1
cc = options.cc - 1

con = sqlite3.connect(locateDatabase())
cur = con.cursor()

numlines = 0
mystr = ""
for line in fileinput.input("-"):
    numlines = numlines+1
    if numlines > 1:
        sys.stderr.write("ERROR: Can only specify ONE cluster/runID pair\n")
        exit(2)
    spl = line.strip("\r\n").split("\t")
    run = spl[rc]
    cluster = spl[cc]
    cur.execute("SELECT DISTINCT organism FROM clusterorgs WHERE runid=? AND clusterid = ?", (run, cluster))
    for res in cur:
        # I store this in a string instead of printing it so that 
        # if we throw an error we don't print anything insidiously.
        # print "".join(str(s) for s in res)
        mystr = "%s\n%s" %("".join(str(s) for s in res), mystr)

print mystr.strip()
con.close()
