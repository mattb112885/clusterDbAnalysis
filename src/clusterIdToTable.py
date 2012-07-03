#!/usr/bin/python

# Utility script
# Pipe in a run ID(s) and specify a cluster ID (or multiple) on the command line
# 
# Outputs to stdount a list of run and cluster ID pairs.
# E.g. pipe in results of
# db_getAllClusterRuns.py | grep "max"
# (Methanosarcina_I_2.0_c_0.4_m_maxbit)
#
# and run
# clusterIdToTable "1000" "1001" "1003"
#
# results will be
# Methanosarcina_I_2.0_c_0.4_m_maxbit  1000
# Methanosarcina_I_2.0_c_0.4_m_maxbit  1001
# Methanosarcina_I_2.0_c_0.4_m_maxbit  1003

import sys
import fileinput
import optparse

usage="%prog [clusterID1] [clusterID2] ... < runID > runID_clusterID_table"
description="Given a RunID from stdin and cluster IDs as arguments, makes a table duplicating the runID and adding the clusterID to each in a tab-delimited table"
parser = optparse.OptionParser(usage=usage, description=description)
(options, args) = parser.parse_args()

if len(args) == 0:
    sys.stderr.write("ERROR: in clusterIdToTable.py - must provide at least one cluster ID\n")
    exit(2)

clusterIds = args[0:]

for line in fileinput.input("-"):
    for s in clusterIds:
        print "%s\t%s" %(line.strip('\r\n'), s)
