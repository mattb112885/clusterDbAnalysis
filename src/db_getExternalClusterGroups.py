#!/usr/bin/python

import fileinput, optparse, sqlite3, sys
from locateDatabase import *

ok_databases = ["all", "cd", "cog", "pfam", "tigr", "prk", "smart"]

usage = "%prog [options] < gene_id_list > External_cluster_ids"
description = """Given a list of gene IDs, identifies the external cluster IDs
that are homologous to those gene IDs as determined by RPSBLAST. By default
returns all of the results - given a specific database to look at, will
only return results from that database."""

parser = optparse.OptionParser(usage=usage, description=description)

parser.add_option("-d", "--database", help="Return only results from the specified database. (D: return all of them). The options are: %s" %(" ".join(ok_databases)),
                  action="store", type="str", dest="database", default="all")
parser.add_option("-a", "--adddescription", help="Add external cluster description to the table (D: dont)", action="store_true", dest="adddescription", default=False)
parser.add_option("-n", "--addname", help="Add external cluster name to the table (D: dont)", action="store_true", dest="addname", default=False)
parser.add_option("-g", "--genecol", help="Column number for gene ID starting from 1 (D: 1)", action="store", type="int", dest="gc", default=1)

(options, args) = parser.parse_args()

gc = options.gc - 1

con = sqlite3.connect(locateDatabase())
cur = con.cursor()

addn = ""
addd = ""
addj = ""
addl = ""
if options.addname:
    addn = " , external_clusters.clustername"
    addj = " INNER JOIN external_clusters ON external_clusters.external_clusterid = rpsblast_results.external_clusterid "
if options.adddescription:
    addd = " , external_clusters.description"
    addj = " INNER JOIN external_clusters ON external_clusters.external_clusterid = rpsblast_results.external_clusterid "
if options.database != "all":
    lk = "\"%s%%\"" %(options.database)
    addl = " AND rpsblast_results.external_clusterid LIKE %s" %(lk)


cmd = """SELECT rpsblast_results.* %s %s FROM rpsblast_results
         %s
         WHERE rpsblast_results.querygene = ? %s """ %(addn, addd, addj, addl)

for line in fileinput.input("-"):
    spl = line.strip("\r\n").split("\t")
    geneid = spl[gc]
    cur.execute(cmd, (geneid, ) )
    for res in cur:
        print "\t".join([ str(s) for s in res ])

con.close()
