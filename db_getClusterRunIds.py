#!/usr/bin/python

# Provide a list of organisms to match [can match any portion of the organism
# so if you give it just "mazei" it will return to you a list of Methanosarcina mazei]
#
# The search organisms must match ALL of the queries (i.e. there are ANDs placed between the wildcards)
#
# Thie function prints a list of cluster IDs to stdout (which can be piped to other commands
# to generate specific data about these clusters, such as all the blast hits for these organisms
# against each other, etc.)
#
# Note - for now I throw an error if there are no cluster datas in the database containing
# exactly the organisms of interest.
#
# We could change this in the future if we only decide to use a single clustering result.

import sqlite3, sys

if len(sys.argv) < 2:
    print "ERROR: In calling db_getClusterIds.py, one must specify at least one organism pattern to match."
    raise IOError

# Change the organism names so that they are all LIKE %name%
# Remove any "%" in the organism name to prevent any issues with %% popping up...
teststr = list('%' + s + '%' for s in sys.argv[1:])

# If the user used double instead of single quotes we want to correct the escaping of the "!" character.
# Otherwise the SQL command will not work because it's trying to match backslashes!
teststr = [ s.replace("\\", "") for s in teststr ]

con = sqlite3.connect("db/methanosarcina")
cur = con.cursor()

# Generate a temporary table (desiredorgs) containing the list of
# organisms for which we want to generate a CORE
query = "CREATE TEMPORARY TABLE desiredorgs AS SELECT DISTINCT organism FROM processed WHERE "
for i in range(len(teststr)):
    if teststr[i].count("!") > 0:
        teststr[i] = teststr[i].replace("!", "")
        query = query + "processed.organism NOT LIKE ? "
    else:
        query = query + "processed.organism LIKE ? "

    if not i == len(teststr) - 1:
        query = query + "AND "
query = query + ";"
cur.execute(query, tuple(teststr))

# Identify which cluster (if any) contains exactly those organisms
# Throw an error if there is none.
#
# NOTE - if we decide to only do clustering on the entire group and not
# separately for separate clades, this query can be replaced by:
#
# SELECT DISTINCT runid FROM distinctorgs;
# (since there should be only one desired ID)
#

# Part 1 - get IDs containing all of the organisms we want (regardless of the total number)
#
# Part 2 - get IDs for clusters containing exactly n organisms
#
# Note - doing HAVING count(*) = count(desiredorgs.organism) doesn't work for some reason.
# Hence the extra weird SELECT statement...
cur.execute(""" SELECT runid FROM distinctorgs
                INNER JOIN desiredorgs ON distinctorgs.organism = desiredorgs.organism
                GROUP BY distinctorgs.runid
                HAVING count(*) = (SELECT count(desiredorgs.organism) FROM desiredorgs)

                INTERSECT

                SELECT runid FROM distinctorgs
                GROUP BY distinctorgs.runid
                HAVING count(*) = (SELECT count(desiredorgs.organism) FROM desiredorgs);""")

atleastone = False
for l in cur:
    atleastone = True
    s = list(l)
    stri = "\t".join(str(t) for t in s)
    print stri

if not atleastone:
    print "ERROR: No clusters found with specified organisms (should be one of the sets found in the 'groups' file for best results)"
    raise IOError

con.close()
