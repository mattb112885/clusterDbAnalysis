#!/usr/bin/python

# Given a list of run IDs (by default in column 1)
# generates a list of run IDs | cluster IDs | gene IDs | Annotations
# (only the first two are needed to pipe to other commands)
#
# Options: -r (--rcolumn) - column number starting from 1 for run ID (default: 1)

import sqlite3, fileinput, optparse

parser = optparse.OptionParser()
parser.add_option("-r", "--rcolumn", help="Column number (start from 1) for run ID", action="store", type="int", dest="runcolumn", default=1)
(options, args) = parser.parse_args()
rc = options.runcolumn - 1 # Convert to Pythonic indexes                                                                                                                                                      

con = sqlite3.connect("db/methanosarcina")
cur = con.cursor()

for line in fileinput.input("-"):

    spl = line.strip().split("\t")
    # Pull out organisms present in the specified clustering run.
    query = "CREATE TEMPORARY TABLE desiredorgs AS SELECT DISTINCT organism FROM distinctorgs WHERE runid = ?;"
    cur.execute(query, (spl[rc],))

    # Generates a table of all the genes for the desired organisms
    cur.execute("CREATE TEMPORARY TABLE orggenes AS SELECT * FROM processed INNER JOIN desiredorgs ON desiredorgs.organism = processed.organism;")
    
    # Generate unique combinations of run IDs, cluster IDs, and organisms for which the organism matches an organism we want
    # and the run ID matches the one we're using at the moment.
    #
    # This combined with the next command will ensure that we only get clusters that contain ALL and ONLY the organisms that we're interested in...
    cur.execute("""CREATE TEMPORARY TABLE unqclusterorgs AS
                      SELECT DISTINCT clusters.runid AS runid, clusters.clusterid AS clusterid, orggenes.organism as organism FROM clusters                                                               
                      INNER JOIN orggenes ON orggenes.geneid = clusters.geneid 
                      INNER JOIN desiredorgs on orggenes.organism = desiredorgs.organism
                      WHERE runid=?;""", (spl[rc], ))

    cur.execute("""CREATE TEMPORARY TABLE clusterids AS 
                    SELECT clusterid FROM unqclusterorgs 
                    GROUP BY clusterid HAVING count(*)=(SELECT count(desiredorgs.organism) FROM desiredorgs)
                   
                    INTERSECT
 
                    SELECT clusterid FROM clusters
                    WHERE runid=?
                    GROUP BY clusterid HAVING count(*)=(SELECT count(desiredorgs.organism) FROM desiredorgs);""", (spl[rc], ) )

    cur.execute("""SELECT clusters.runid, clusters.clusterid, clusters.geneid, processed.annotation
                      FROM clusters
                      INNER JOIN processed ON processed.geneid = clusters.geneid
                      INNER JOIN clusterids ON clusters.clusterid = clusterids.clusterid
                      WHERE clusters.runid=?;""", (spl[rc], ))

    for l in cur:
        s = list(l)
        stri = "\t".join(str(t) for t in s)
        print stri

    cur.execute("DROP TABLE clusterids;")
    cur.execute("DROP TABLE desiredorgs;")
    cur.execute("DROP TABLE orggenes;")
    cur.execute("DROP TABLE unqclusterorgs;")

con.close()
