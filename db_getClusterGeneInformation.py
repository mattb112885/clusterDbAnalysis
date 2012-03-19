#!/usr/bin/python

# This is a pipe command
# Pipe in a table with two columns:
# - a column with the run ID one wishes to pull out information for (default: 1st column)
# - a column with the cluster ID one wishes to pull out information for (default: second column)
#
# This command will generate information about ALL genes in the specified clusters, not
#   just the ones specified
#
# Will generate a table containing, for all genes in the specified clusters / run IDs:
# - run ID
# - cluster ID
# - gene IDs
# - organism
# - (later) if it is a pseudogene or not
# - annotation (from RAST)
# - amino acid sequence
# - nucleic acid sequence
# (note - later on I'll probably add the option to specify which of these you want, like BLAST)
#
# Optional arguments:
# -r or --rcolumn: Column number for run id (count from 1)
# -c or --ccolumn: Column number for cluster id (count from 1)

import sqlite3, fileinput, optparse

# Get input arguments                                                                                                                                                                                          
parser = optparse.OptionParser()
parser.add_option("-r", "--rcolumn", help="Column number (start from 1) for run ID", action="store", type="int", dest="runcolumn", default=1)
parser.add_option("-c", "--ccolumn", help="Column number (start from 1) for cluster ID", action="store", type="int", dest="clustercolumn", default=2)
(options, args) = parser.parse_args()
rc = options.runcolumn - 1 # Convert to Pythonic indexes                                                                                                                                                      
cc = options.clustercolumn - 1

con = sqlite3.connect("db/methanosarcina")
cur = con.cursor()

# Make a temporary table with all of the clusters that we want to actually extract
query = "CREATE TEMPORARY TABLE desiredclusters( runid VARCHAR(256), clusterid INT );"
cur.execute(query);

for line in fileinput.input("-"):
    spl = line.strip().split("\t")
    query = "INSERT INTO desiredclusters VALUES (?, ?);"
    con.execute(query, (spl[rc], spl[cc]))

# Generate information about all of those clusters (for all genes, not just the ones piped into this command)
cur.execute(""" SELECT clusters.*, processed.organism, processed.annotation, processed.aaseq, processed.nucseq
              FROM clusters
              INNER JOIN processed ON processed.geneid = clusters.geneid
              INNER JOIN desiredclusters ON desiredclusters.clusterid = clusters.clusterid
                 AND clusters.runid = desiredclusters.runid;""")

for l in cur:
    s = list(l)
    stri = "\t".join(str(t) for t in s)
    print stri

con.close()
