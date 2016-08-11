#!/usr/bin/env python

import fileinput, optparse, sys, sqlite3
from FileLocator import *
from ClusterFuncs import *

usage="%prog [options] < runid > organism_list"
description="""Get a list of organisms included in each piped-in cluster run 
(Note - the results are most useful if you only provide ONE)"""
parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-r", "--runcol", help="Column number for run ID starting from 1 (D=1)", action="store", type="int", dest="rc", default=1)
(options, args) = parser.parse_args()

rc = options.rc - 1

con = sqlite3.connect(locateDatabase())
cur = con.cursor()

for line in fileinput.input("-"):
    spl = line.strip("\r\n").split("\t")
    runid = spl[rc]
    orglist = getOrganismsInClusterRun(runid, cur)
    for org in orglist:
        print("\t".join( [ runid, org ] ))

con.close()
