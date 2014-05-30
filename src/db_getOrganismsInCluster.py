#!/usr/bin/env python

import fileinput, optparse, sys, sqlite3
from FileLocator import *
from ClusterFuncs import *

header = [ "organism_name" ]

usage="""%prog [options] < runid_clusterid_pair > organism_list

Output: """ + " ".join(header)
description="Get a list of organisms included in the specified run / cluster ID pair"
parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-r", "--runcol", help="Column number for run ID starting from 1 (D=1)", action="store", type="int", dest="rc", default=1)
parser.add_option("-c", "--clustercol", help="Column number for cluster ID starting from 1 (D=2)", action="store", type="int", dest="cc", default=2)
parser.add_option("--header", help="Specify to add header to the output file (useful if you want to take the results and put into Excel or similar programs)",
                  action="store_true", default=False)
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
    runid = spl[rc]
    clusterid = spl[cc]
    orglist = getOrganismsInCluster(runid, clusterid, cur)
    stri = "\n".join(orglist)

if options.header:
    print "\t".join(header)
print stri

con.close()
