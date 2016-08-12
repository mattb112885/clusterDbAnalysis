#!/usr/bin/env python

# Get clusters (in the same run) with the same presence/absence pattern as
# the cluster/runID pair as specified...

from __future__ import print_function
import fileinput, optparse, sys, os

usage="%prog [options] < runid_clusterid_pair > same_clusters"
description="""Generate a list of clusters with the same organism presence/absence pattern as the specified cluster/runID pair
Only looks in the same run as specified, and does not account for the number of genes present in each organism."""

parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-r", "--runcol", help="Column number for run ID starting from 1 (D=1)", action="store", type="int", dest="rc", default=1)
parser.add_option("-c", "--clustercol", help="Column number for cluster ID starting from 1 (D=2)", action="store", type="int", dest="cc", default=2)
(options, args) = parser.parse_args()

rc = options.rc - 1
cc = options.cc - 1

for line in fileinput.input("-"):
    spl = line.strip("\r\n").split("\t")
    run = spl[rc]
    cluster = spl[cc]
    # Output to stdout anyway so we don't change it here...
    # We want exactly the same pattern so we do -a -s (ALL and ONLY)
    cmd = """makeTabDelimitedRow.py "%s" "%s" | db_getOrganismsInCluster.py | db_findClustersByOrganismList.py -a -s "%s" """ %(run, cluster, run)
    os.system(cmd)
