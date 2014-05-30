#!/usr/bin/env python

# Returns a list of all BLAST results
# available in the database.
# 

import sqlite3
import optparse
from FileLocator import *

header = [ 'Query_gene', 'Target_gene', 'Percent_identity', 'HSP_length', 'Percent_mismatch',
           'Gap_opens', 'Query_start', 'Query_end', 'Target_start', 'Target_end', 'E_value',
           'Bit_score', 'Query_selfbit', 'Target_selfbit' ]

usage="""%prog > all_blast_results

Output: """ + " ".join(header)
description="Print all blast results available in the database (without further filtering)"
parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("--header", help="Specify to add header to the output file (useful if you want to take the results and put into Excel or similar programs)",
                  action="store_true", default=False)
(options, args)=parser.parse_args()

con = sqlite3.connect(locateDatabase())
cur = con.cursor()

# Generate a list of blast results with query matching one of the desiredgenes
cur.execute("""SELECT * FROM blastres_selfbit;""")

if options.header:
    print "\t".join(header)
for l in cur:
    s = list(l)
    stri = "\t".join(str(t) for t in s)
    print stri

con.close()
