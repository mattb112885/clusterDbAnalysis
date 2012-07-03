#!/usr/bin/python

# This is a pipe command.
#
# Unlike db_getBlastResultsBetweenSpecificGenes (which only gets blast results containing both
# query and target inside them),
# this command will get blast results if either the query OR the target
# is contained within them.
#
# The code is identical except the sql command has an OR instead of an AND...

import fileinput, optparse, sqlite3
from locateDatabase import *

usage = "%prog [options] < gene_ids > blast_results"
description = "Given list of genes to match, returns a list of BLAST results containing any gene ID in your list as either a query or a target (for blast results only BETWEEN the query genes, see db_getBlastResultsBetweenSpecificGenes.py)"
parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-g", "--gcolumn", help="Column number (start from 1) for gene ID", action="store", type="int", dest="genecolumn", default=1)
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
cur.execute("""SELECT blastres_selfbit.* FROM blastres_selfbit
               WHERE blastres_selfbit.targetgene IN (select geneid from desiredgenes)
               OR blastres_selfbit.querygene IN (select geneid from desiredgenes);""");

for l in cur:
    s = list(l)
    stri = "\t".join(str(t) for t in s)
    print stri

con.close()
