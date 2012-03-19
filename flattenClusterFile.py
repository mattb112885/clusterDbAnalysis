#!/usr/bin/python

# This is a pipe command.
# 
# Pipe in the cluster file (from MCL)
# and the flattened file comes out in stdout
#
# The "flattened file" contains only three columns:
# - A run ID (which should be unique)
# - Cluster number 
# - Gene ID
#
# This will make it easy to generate SQL queries for particular clusters
# (e.g. give me the blast results for all of the genes listed in cluster 1
# against anything else)

import fileinput, random, string

# 64 random bits should be enough to uniquely distinguish two runs...
# If we have a cluster ID we can then determine by SQL query what organisms were
#    utilized in that cluster if needed.
clusterid = ''.join(random.choice(string.ascii_uppercase + string.digits) for x in range(64))


rowIdx = 1
for line in fileinput.input("-"):
    spl = line.strip().split('\t')
    for gene in spl:
        string = str(clusterid) + "\t" + str(rowIdx) + "\t" + gene
        print string
    rowIdx = rowIdx + 1
