#!/usr/bin/env python

import fileinput
import optparse
import sqlite3
import sys

from FileLocator import *

header = [ 'CDD_ID', 'Domain_ID', 'Domain_name', 'Description', 'Alignment_length' ]
usage = """%prog [options] < external_clusterids > descriptions

Output: """ + " ".join(header)
description = "Get descriptions of external clusters described by the specified external clusterIDs"
parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-c", "--column", help="Column number (starting from 1) for the external cluster ID (D=1)", 
                  action="store", type="int", dest="idc", default=1)
parser.add_option("--header", help="Specify to add header to the output file (useful if you want to take the results and put into Excel or similar programs)",
                  action="store_true", default=False)
(options, args) = parser.parse_args()

idc = options.idc - 1

con = sqlite3.connect(locateDatabase())
cur = con.cursor()
query = """SELECT * FROM external_clusters WHERE external_clusterid = ?"""

if options.header:
    print "\t".join(header)
for line in fileinput.input("-"):
    spl = line.strip("\r\n").split("\t")
    eid = spl[idc]
    cur.execute(query, (eid, ))
    for res in cur:
        print "\t".join([str(s) for s in res])
