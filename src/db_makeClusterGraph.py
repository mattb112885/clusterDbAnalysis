#!/usr/bin/env python

import fileinput
import optparse
import sqlite3
import sys

from ClusterGraph import *
from FileLocator import *

symmetric_methods, unsymmetric_methods = getValidBlastScoreMethods()
all_methods = symmetric_methods + unsymmetric_methods

usage = """ %prog -m method -u cutoff [options] < cluster_runid_pairs """
description = """From a given cluster/runID pair, calcualte and export a GML file that
can be opened in Cytoscape or similar viewers to visualize the cluster.
By default, it creates GML files with name "runid_clusterid.gml" . If you specify your own
GML file you can only provide ONE cluster/runID pair.
Currently implemented methods: %s """ %(" ".join(all_methods) )

parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-m", "--method", help="Scoring method to use (REQUIRED)", action="store", type="str", dest="method", default=None)
parser.add_option("-u", "--cutoff", help="Scoring cutoff to use (REQUIRED)", action="store", type="float", dest="cutoff", default=None)

parser.add_option("-r", "--runcol", help="Column number for Run ID starting from 1 (D:1)", action="store", type="int", dest="rc", default=1)
parser.add_option("-c", "--clusterid", help="Column number for Cluster ID starting from 1 (D: 2)", action="store", type="int", dest="cc", default=2)
parser.add_option("-n", "--blastn", help="Use BLASTN instead of BLASTP for scoring results (D: BLASTP)", action="store_true", dest="blastn", default=False)
(options, args) = parser.parse_args()

# Convert to pythonic indexes
rc = options.rc - 1
cc = options.cc - 1

if options.method is None or options.cutoff is None:
    sys.stderr.write("ERROR: both method (-m) and cutoff (-u) are required arguments\n")
    exit(2)

cluster_runs = []
for line in fileinput.input("-"):
    spl = line.strip("\r\n").split("\t")
    cluster_runs.append( ( spl[rc], spl[cc] ) )

con = sqlite3.connect(locateDatabase())
cur = con.cursor()

for cr in cluster_runs:
    runid = cr[0]
    clusterid = cr[1]
    G = getGraphForCluster( runid, clusterid, options.method, options.cutoff, cur, blastn = options.blastn )
    if options.blastn:
        blast_str = "blastn"
    else:
        blast_str = "blastp"
    exportGraphToGML(G, "%s_%s_%s_%1.2f_%s.gml" %(runid, clusterid, options.method, options.cutoff, blast_str ) )

con.close()
