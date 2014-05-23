#!/usr/bin/env python

import fileinput
import optparse
import sqlite3
import sys

from FileLocator import *

header = " ".join([ "querygene", "cdd_id", "pctid", "alnlen", "mismatches", "gapopens", "querystart", "queryend", "substart", "subend", "evalue", "bitscore" ])

usage = """%prog [options] < external_clusterids > similarities

Output: """ + header
description = """Given a list of external clusterIDs (from cog, pfam, tigrfam...), attempts to identify proteins
in the database that are homologous (by RPSBLAST) to the profile implied by those clusters, and returns them
and their E-values."""

parser = optparse.OptionParser(usage=usage, description=description)

parser.add_option("-c", "--column", help="Column number for external cluster ID (D=1)",
                  action = "store", type="int", dest="ecc", default=1)
parser.add_option("-e", "--evalue", help="Evalue cutoff (D: 1E-5, the cutoff used to originally generate the data)", action="store",
                  type="float", dest="evalue", default=1E-5)

(options, args) = parser.parse_args()

ecc = options.ecc - 1

con = sqlite3.connect(locateDatabase())
cur = con.cursor()

q = """SELECT * FROM rpsblast_results
     INNER JOIN external_clusters ON external_clusters.cdd_id = rpsblast_results.cdd_id
     WHERE external_clusterid = ? AND evalue < ?"""

for line in fileinput.input("-"):
    spl = line.strip("\r\n").split("\t")
    externalid = spl[ecc]
    cur.execute(q, (externalid, options.evalue))
    for res in cur:
        print "\t".join([ str(s) for s in res ] )

con.close()
