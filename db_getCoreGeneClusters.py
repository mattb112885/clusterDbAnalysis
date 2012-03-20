#!/usr/bin/python

# This is a pipe command.
#
# Pipe in a list of organisms (obtain e.g. by grepping the "organisms" file for
# the subset of organisms that you want)
#
# The function first exclude any organisms NOT in the list from the clusters
# in all of the runs. It then generates a list of core clusters according to our
# definition:
#
# Core gene cluster is defined as a cluster containing exactly one gene from each organism.
#

import sqlite3, fileinput, re

# FIXME: Option for organism column

con = sqlite3.connect("db/methanosarcina")
cur = con.cursor()

oc = 0

orgs = []
for line in fileinput.input("-"):
    spl = line.strip().split("\t")
    orgs.append(spl[oc])

cur.execute("DROP TABLE IF EXISTS myorganisms;")
cur.execute("DROP TABLE IF EXISTS clusterorgs_my;")

# Desired organism table
cur.execute("CREATE TABLE myorganisms(organism VARCHAR(256));")
for org in orgs:
    cur.execute("INSERT INTO myorganisms VALUES (?);", (org, ) )

# Filter clusterorgs by which organism we want
query = """CREATE TABLE clusterorgs_my AS
                SELECT * FROM clusterorgs
                WHERE clusterorgs.organism IN myorganisms;"""

cur.execute(query)

# Speedup
cur.execute("CREATE INDEX myorgrunidx ON clusterorgs_my(runid);")
cur.execute("CREATE INDEX myorgclusteridx ON clusterorgs_my(clusterid);")

checkPairs = []
cur.execute("SELECT DISTINCT clusterorgs_my.runid, clusterorgs_my.clusterid FROM clusterorgs_my;")
for l in cur:
    s = list(l)
    checkPairs.append( [s[0], s[1]] )

for p in checkPairs:
    cur.execute("""SELECT geneid FROM clusterorgs_my
                    WHERE clusterorgs_my.runid = ?
                     AND clusterorgs_my.clusterid = ?;""",
                (p[0], p[1]) )
    genelist = set()
    orglist =set()
    for l in cur:
        s = list(l)
        genelist.add(l[0])
        # Remove peg.#### from the names of the genes to get organism list
        orglist.add(re.sub("\.peg\.\d+", "", l[0]))

    # Definition of core gene
    if len(genelist) == len(orglist) and len(genelist) == len(orgs):
        print "%s\t%s" % (p[0], p[1])

cur.execute("DROP TABLE myorganisms;")
cur.execute("DROP TABLE clusterorgs_my;")

con.close()
