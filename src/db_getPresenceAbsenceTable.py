#!/usr/bin/python

import sqlite3
import sys
import optparse
from locateDatabase import *

usage="%prog [options] > presence_absence_table"
description="""Generates a presence\absence table (or slices thereof) based on
the one automatically loaded as part of main2.sh. 
Default activity is to dump the database as is (pegs)."""
parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-n", "--number", help="Rather than printing PEGs, print the number of representatives in each organism (D: Prints pegs)", 
                  action="store_true", dest="number", default=False)
parser.add_option("-b", "--binary", help="Rather than printing PEGs, print 0 if there are no representatives and 1 if there are representatives (D: prints pegs)",
                  action = "store_true", dest="binary", default=False)
parser.add_option("-r", "--runid", help="Only print results for the specified run ID (D: Prints the table for all of them)", action="store", type="str", dest="runid", default=None)
(options,args) = parser.parse_args()

if options.number and options.binary:
    sys.stderr.write("ERROR: Cannot specify both -n and -b (can either print 0\1 or number of represenatitives, not both)\n")
    exit(2)

con = sqlite3.connect(locateDatabase())
cur = con.cursor()

if options.runid is None:
    cur.execute("SELECT * FROM presenceabsence;")
else:
    cur.execute("SELECT * FROM presenceabsence WHERE runid = ?", (options.runid,))

collist = [tup[0] for tup in cur.description]
print "\t".join(collist)
for rec in cur:
    lst = [ str(s) for s in rec ]
    for ii in range(len(lst)):
        if ii < 3:
            continue
        if (options.number or options.binary) and lst[ii] == "NONE":
            lst[ii] = "0"
        elif options.binary:
            lst[ii] = "1"
        elif options.number:
            lst[ii] = str(lst[ii].count(";") + 1)
    print "\t".join(lst)
