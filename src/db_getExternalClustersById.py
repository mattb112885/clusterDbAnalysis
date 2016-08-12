#!/usr/bin/env python

from __future__ import print_function
import fileinput
import optparse
import sqlite3
import sys

from FileLocator import *

usage = "%prog [options] < external_clusterids > descriptions"
description = "Get descriptions of external clusters described by the specified external clusterIDs"
parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-c", "--column", help="Column number (starting from 1) for the external cluster ID (D=1)", 
                  action="store", type="int", dest="idc", default=1)
(options, args) = parser.parse_args()

idc = options.idc - 1

con = sqlite3.connect(locateDatabase())
cur = con.cursor()
query = """SELECT * FROM external_clusters WHERE external_clusterid = ?"""

for line in fileinput.input("-"):
    spl = line.strip("\r\n").split("\t")
    eid = spl[idc]
    cur.execute(query, (eid, ))
    for res in cur:
        print("\t".join([str(s) for s in res]))
