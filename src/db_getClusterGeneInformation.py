#!/usr/bin/env python

# This is a pipe command
# Pipe in a table with two columns:
# - a column with the run ID one wishes to pull out information for (default: 1st column)
# - a column with the cluster ID one wishes to pull out information for (default: second column)
#
# This command will generate information about ALL genes in the specified clusters, not
#   just the ones specified
#
# Will generate a table containing, for all genes in the specified clusters / run IDs:
# - run ID
# - cluster ID
# - gene IDs
# - organism
# - (later) if it is a pseudogene or not
# - annotation (from RAST)
# - amino acid sequence
# - nucleic acid sequence
# (note - later on I'll probably add the option to specify which of these you want, like BLAST)
#
# Optional arguments:
# -r or --rcolumn: Column number for run id (count from 1)
# -c or --ccolumn: Column number for cluster id (count from 1)

import sqlite3, fileinput, optparse
from FileLocator import *
from ClusterFuncs import *

# Get input arguments                                                                  

headers = [ "organism_name", "contig_id", "start", "stop", "strand", "strandnum", "annotation", "DNA_seq", "AA_seq", "run_id", "cluster_id" ]
usage = """%prog [options] < runid_clusterid_table > cluster_gene_info

Output table: """ + " ".join(headers)

description = """Given a list of run ID / cluster ID pairs (one pair in each row of the input table), 
get a list of gene information for each gene in each of the input clusters. The results include
organism, strand, location, contig, and sequences (and the run\clusterID to which the gene belongs)""" 

parser = optparse.OptionParser(usage=usage, description=description)

parser.add_option("-r", "--rcolumn", help="Column number (start from 1) for run ID", action="store", type="int", dest="runcolumn", default=1)
parser.add_option("-c", "--ccolumn", help="Column number (start from 1) for cluster ID", action="store", type="int", dest="clustercolumn", default=2)
(options, args) = parser.parse_args()

rc = options.runcolumn - 1 # Convert to Pythonic indexes                                                                                                                                                      
cc = options.clustercolumn - 1


clusterruns = set()
for line in fileinput.input("-"):
    spl = line.strip('\r\n').split("\t")
    clusterruns.add( (spl[rc], spl[cc]) )

con = sqlite3.connect(locateDatabase())
cur = con.cursor()

for cr in clusterruns:
    geneinfo = getClusterGeneInfo(cr[0], cr[1], cur)
    for info in geneinfo:
        print("\t".join(info))

con.close()
