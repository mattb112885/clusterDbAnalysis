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
parser.add_option("-c", "--cutoff", help="E-value cutoff (D: Show all results in database)", action="store", type="float", dest="cutoff", default=10)
parser.add_option("-n", "--blastn", help="Base the results on BLASTN instead of BLASTP (D: BLASTP)", action="store_true", dest="blastn", default=False)
(options, args) = parser.parse_args()

gc = options.genecolumn - 1

con = sqlite3.connect(locateDatabase())
cur = con.cursor()

# Generate a table of BLAST results                                                                                                                                                                            
cur.execute("""CREATE TEMPORARY TABLE desiredgenes ("geneid" VARCHAR(128), FOREIGN KEY(geneid) REFERENCES rawdata(geneid));""")
for line in fileinput.input("-"):
    spl = line.strip('\r\n').split("\t")
    gn = spl[gc]
    # fig| is optional...
    if not gn.startswith("fig|"):
        gn = "fig|%s" %(gn)
    cur.execute("INSERT INTO desiredgenes VALUES (?);", (gn, ) )

# Generate a list of blast results with query matching one of the desiredgenes

if options.blastn:
    tbl = "blastnres_selfbit"
else:
    tbl = "blastres_selfbit"

cmd = """SELECT %s.* FROM %s
         WHERE (%s.evalue < ?) AND 
             ( %s.targetgene   IN (select geneid from desiredgenes)
               OR %s.querygene IN (select geneid from desiredgenes) ) ;""" %(tbl, tbl, tbl, tbl, tbl)

cur.execute(cmd, (options.cutoff,));

for l in cur:
    s = list(l)
    stri = "\t".join(str(t) for t in s)
    print stri

con.close()
