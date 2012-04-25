#!/usr/bin/python

# Utility script
# Pipe in a run ID(s) and specify a cluster ID (or multiple) on the command line
# 
# Outputs to stdount a list of run and cluster ID pairs.
# E.g. pipe in results of
# ./src/db_getAllClusterRuns.py | grep "max"
# (Methanosarcina_I_2.0_c_0.4_m_maxbit)
#
# and run
# ./src/clusterIdToTable "1000" "1001" "1003"
#
# results will be
# Methanosarcina_I_2.0_c_0.4_m_maxbit  1000
# Methanosarcina_I_2.0_c_0.4_m_maxbit  1001
# Methanosarcina_I_2.0_c_0.4_m_maxbit  1003

import sys
import fileinput

if len(sys.argv) < 2:
    print "ERROR: in clusterIdToTable.py - must provide at least one cluster ID"
    exit(2)

clusterIds = sys.argv[1:]

for line in fileinput.input("-"):
    for s in clusterIds:
        print "%s\t%s" %(line.strip(), s)
