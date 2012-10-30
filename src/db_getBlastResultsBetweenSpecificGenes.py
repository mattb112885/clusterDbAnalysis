#!/usr/bin/python

# Provide a list of organisms to match [can match any portion of the organism
# so if you give it just "mazei" it will return to you a list of Methanosarcina mazei]
#
# Returns a list of blast results specific to those organisms to stdout
# (which can subsequently be used to do clustering...)
# 

import sqlite3, optparse, fileinput
from locateDatabase import *

usage = "%prog [options] < gene_ids > blast_results"
description = "Given list of genes to match, returns a list of BLAST results between genes in the list only"
parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-g", "--gcolumn", help="Column number (start from 1) for gene ID", action="store", type="int", dest="genecolumn", default=1)
parser.add_option("-n", "--blastn", help="Get BLASTN results (D: BLASTP results)", action="store_true", dest="blastn", default=False)
(options, args) = parser.parse_args()
gc = options.genecolumn - 1

con = sqlite3.connect(locateDatabase())
cur = con.cursor()

# Generate a table of BLAST results
cur.execute("""CREATE TEMPORARY TABLE desiredgenes ("geneid" VARCHAR(128), FOREIGN KEY(geneid) REFERENCES rawdata(geneid));""")
for line in fileinput.input("-"):
    spl = line.strip('\r\n').split("\t")
    cur.execute("INSERT INTO desiredgenes VALUES (?);", (spl[gc], ) )


# Generate a list of blast results with query matching one of the desiredgenes
if options.blastn:
    tbl = "blastnres_selfbit"
else:
    tbl = "blastres_selfbit"

cmd = """SELECT %s.* FROM %s
         WHERE %s.targetgene IN (select geneid from desiredgenes)
         AND %s.querygene IN (select geneid from desiredgenes);""" %(tbl, tbl, tbl, tbl)

cur.execute(cmd)

for l in cur:
    s = list(l)
    stri = "\t".join(str(t) for t in s)
    print stri

con.close()
