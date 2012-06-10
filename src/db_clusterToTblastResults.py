#!/usr/bin/python

# This is a pipe command. Pipe in a 2-column table
# containing the run ID in the first column and the
# cluster ID in the second column (DEFAULT) OR
# 
# -r [runcolumn] -c [clustercolumn]
# (starting from 1 as the first column)
# (Note - type db_clusterToTblastResults --help for help details)

import fileinput, sqlite3, optparse

# Get input arguments
usage = "%prog [options]"
description = "Given a list of run IDs and column IDs (both from stdin), returns all of the TBLASTN results between genes in the same clusters"
parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-r", "--rcolumn", help="Column number (start from 1) for run ID (default = 1)", action="store", type="int", dest="runcolumn", default=1)
parser.add_option("-c", "--ccolumn", help="Column number (start from 1) for cluster ID (default = 2)", action="store", type="int", dest="clustercolumn", default=2)

(options, args) = parser.parse_args()
rc = options.runcolumn - 1 # Convert to Pythonic indexes
cc = options.clustercolumn - 1

con = sqlite3.connect("db/methanosarcina")
cur = con.cursor()

# Make a temporary table with all of the clusters that we want to actually extract

query = "CREATE TEMPORARY TABLE allclusters( runid VARCHAR(256), clusterid INT );"
cur.execute(query);
for line in fileinput.input("-"):
    spl = line.strip('\r\n').split("\t")
    query = "INSERT INTO allclusters VALUES (?, ?);"
    con.execute(query, (spl[rc], spl[cc]))
    
# Remove duplication... especially important since we split this thing up into multiple queries now.
arglist = []
cur.execute("CREATE TEMPORARY TABLE desiredclusters AS SELECT DISTINCT * FROM allclusters;")
cur.execute("SELECT * FROM desiredclusters;")
for l in cur:
    spl = list(l)
    arglist.append( (spl[rc], spl[cc], spl[rc], spl[cc], spl[rc], spl[cc] ) )

# Get table of gene IDs
# SELECT DISTINCT because the input desiredclusters will often (usually)
# have multiple copies of the same desired cluster ID.
cur.execute("""CREATE TEMPORARY TABLE desiredgenes AS
               SELECT DISTINCT clusters.* FROM clusters
               INNER JOIN desiredclusters ON desiredclusters.runid = clusters.runid
                      AND desiredclusters.clusterid = clusters.clusterid;""")

# This performs all the necessary filters for the tblastn results
# except the requirement that both genes have to be in the same cluster / run
#
# Note - I join the desiredgenes here because that join can take a long time and I don't want to
# do it repeatedly in a loop (below).
cur.execute("""CREATE TEMPORARY TABLE joined AS
               SELECT desiredgenes.runid, desiredgenes.clusterid, 
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
               INNER JOIN desiredgenes ON desiredgenes.geneid = tblast_converted.querygene;""")

cur.execute("""CREATE INDEX joined_querygene ON joined(querygene);""")
cur.execute("""CREATE INDEX joined_targetgene ON joined(targetgene);""")
cur.execute("""CREATE INDEX joined_clusterid ON joined(clusterid);""")
cur.execute("""CREATE INDEX joined_runid ON joined(runid);""")

# Require that both the query and target gene are in the correct cluster / run
for args in arglist:
    cur.execute(""" SELECT joined.* FROM joined
                    WHERE joined.runid = ? AND joined.clusterid = ?
                          AND joined.querygene IN (SELECT geneid FROM desiredgenes WHERE desiredgenes.runid = ? AND desiredgenes.clusterid = ?)
                          AND joined.targetgene IN (SELECT geneid FROM desiredgenes WHERE desiredgenes.runid = ? AND desiredgenes.clusterid = ?)
                          ;""", args)
    for l in cur:
        s = list(l)
        stri = "\t".join(str(t) for t in s)
        print stri
