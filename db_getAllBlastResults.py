#!/usr/bin/python

# Returns a list of all BLAST results
# available in the database.
# 

import sqlite3

con = sqlite3.connect("db/methanosarcina")
cur = con.cursor()

# Generate a list of blast results with query matching one of the desiredgenes
cur.execute("""SELECT * FROM blastres_selfbit;""")

for l in cur:
    s = list(l)
    stri = "\t".join(str(t) for t in s)
    print stri

con.close()
