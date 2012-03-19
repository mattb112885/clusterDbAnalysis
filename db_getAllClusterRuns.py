#!/usr/bin/python

# Takes no inputs.
#
# Output: A list of all run IDs from the database (if one does not desire to filter by run ID, and instead
# wishes to generate lists for ALL of the run IDs and process them later,
# one should call this function)
#
# Must be run from the root directory

import sqlite3, optparse

usage = "%prog [options]"
description = "Return list of all run IDs from the database"
parser = optparse.OptionParser(usage=usage, description=description)
(options, args) = parser.parse_args()

con = sqlite3.connect("db/methanosarcina")
cur = con.cursor()

cur.execute("SELECT DISTINCT runid FROM clusters;")

for l in cur:
    s = list(l)
    stri = "\t".join(str(t) for t in s)
    print stri

con.close()
