#!/usr/bin/python

# This is a pipe command. Pipe in a list of 
# genes between which to find TBLASTN results
# Default: column 1
# Can also specify column with
# -g [genecolumn]
# (starting from 1 as the first column)

import fileinput, sqlite3, optparse

# Get input arguments
description = "Given a list of genes (from stdin), finds all TBLASTN results between those genes only"
parser = optparse.OptionParser()
parser.add_option("-g", "--gcolumn", help="Column number (start from 1) for gene ID", action="store", type="int", dest="genecolumn", default=1)
(options, args) = parser.parse_args()
gc = options.runcolumn - 1 # Convert to Pythonic indexes

con = sqlite3.connect("db/methanosarcina")
cur = con.cursor()

# Get list of gene IDs
geneids = []
for line in fileinput.input("-"):
    spl = line.strip('\r\n').split("\t")
    geneids.append(spl[gc])

# Unique IDs
geneids = set(geneids)

# Generate a table of gene IDs to join
cur.execute("""CREATE TEMPORARY TABLE desiredgenes (
                 "geneid" VARCHAR(128) 
                  FOREIGN KEY(geneid) REFERENCES rawdata(geneid));""")

for g in geneids:
    cur.execute("""INSERT INTO desiredgenes VALUES (?);""", (g, ))

# Join them up
cur.execute("""SELECT 
                   tblast_converted.querygene, processed.genestart AS querystart, processed.geneend AS queryend, 
                        (processed.geneend - processed.genestart) AS querydiff, processed.annotation, processed.organism,
                   
                   tblast_converted.targetgene, tblast_converted.targetorganism, tblast_converted.targetannotation,
                   tblast_converted.tblast_targetstart, tblast_converted.tblast_targetend, 
                        (tblast_converted.tblast_targetend - tblast_converted.tblast_targetstart) AS tblastn_targetdiff,
                   tblast_converted.actual_targetstart, tblast_converted.actual_targetend, 
                        (tblast_converted.actual_targetend - tblast_converted.actual_targetstart) AS actual_targetdiff,

                   tblast_converted.pctid, tblast_converted.evalue

               FROM tblast_converted
               INNER JOIN processed ON tblast_converted.querygene = processed.geneid
               INNER JOIN desiredgenes ON desiredgenes.geneid = tblast_converted.querygene
               WHERE tblast_converted.targetgene IN (SELECT geneid FROM desiredgenes);""")

for l in cur:
    s = list(l)
    # Two X columns allow us to pipe the results from this into filterTblastResults.py
    stri = "X\tX\t" + "\t".join(str(t) for t in s)
    print stri
