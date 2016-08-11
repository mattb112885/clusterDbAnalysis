#!/usr/bin/env python

# Returns a list of all BLAST results
# available in the database.
# 

import sqlite3
import optparse
from FileLocator import *

usage="%prog > all_blast_results"
description="Print all blast results available in the database (without further filtering)"
parser = optparse.OptionParser(usage=usage, description=description)
(options, args)=parser.parse_args()

con = sqlite3.connect(locateDatabase())
cur = con.cursor()

# Generate a list of blast results with query matching one of the desiredgenes
cur.execute("""SELECT * FROM blastres_selfbit;""")

for l in cur:
    s = list(l)
    stri = "\t".join(str(t) for t in s)
    print(stri)

con.close()
