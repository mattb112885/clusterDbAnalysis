#!/usr/bin/python

# This is a pipe command.
# Specify the column containing a gene ID
# and we will add the organism's name corresponding to that gene ID as an
# extra column on the output.

import sqlite3
import fileinput
import optparse

usage="%prog [options] < gene_ids > gene_ids_with_organism"
description = "Add the organism name to a tab-delimited file containing gene IDs"
parser = optparse.OptionParser(description=description, usage=usage)
parser.add_option("-g", "--genecol", help="Column number for gene IDs starting from 1 (D=1)", action="store", type="int", dest="genecol", default=1)
(options, args) = parser.parse_args()

gc = options.genecol - 1

con = sqlite3.connect("db/methanosarcina")
cur = con.cursor()

for line in fileinput.input("-"):
    spl = line.strip('\r\n').split("\t")
    cur.execute("SELECT organism FROM processed WHERE processed.geneid=?;", (spl[gc],))
    for k in cur:
        org = str(k[0])
        spl.append(org)
        print "\t".join(spl)

con.close()
