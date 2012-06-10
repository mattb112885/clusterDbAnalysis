#!/usr/bin/python

# Compare all of the clusters between two specified cluster runs
#
#
# Output format:
#
# Run ID 1 | Cluster ID 1 | Run ID 2 | Cluster ID 2 | Overlapping genes | Genes only in (RunID1,clusterID1) | Genes only in (RunID2,ClusterID2)
#
# Any genes that aren't in the other run at all are ignored unless -i is specified, in which case they are included as non-overlapping genes
# and clusters containing only genes not in the other run are given a file corresponding to a dummy cluster (NO_OVERLAP) in the other run ID

import optparse
import sqlite3
import sys

usage = "%prog RUNID1 RUNID2 > comparison_file"
description="Generates a file comparing each overlapping pair of clusters between the two specified run IDs. The tab-delimited file has one row for each overlapping pair of clusters and has one column with the overlapping genes and one column with the non-overlapping genes from each cluster"
parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-i", "--includeall", help="If specified, include organisms not in common between the two runs as non-overlapping genes. Default is to omit those genes", action="store_true", dest="includeall", default=False)
(options, args) = parser.parse_args()

if len(args) != 2:
    sys.stderr.write("ERROR: Must specify exactly two run IDs to compare to each other\n")
    exit(2)

##################
#### Aux functions
##################
def getClusterRunDictionaries(runid, cur):
    cur.execute("SELECT clusterid,geneid FROM clusters WHERE runid = ?;", ( runid, ) )
    #Should be one to one
    geneToCluster = {}
    #One to many
    clusterToGenes = {}
    for ln in cur:
        ls = [ str(s) for s in ln ]
        geneToCluster[ls[1]] = ls[0]
        if ls[0] in clusterToGenes:
            clusterToGenes[ls[0]].append(ls[1])
        else:
            clusterToGenes[ls[0]] = [ ls[1] ]
    return geneToCluster, clusterToGenes

# Check if there exists a gene in the specified cluster in run 1 that also exists in run2.
def overlapExists(run1_clusterid, run1_clusterToGenes, run2_geneToCluster):
    for gene in run1_clusterToGenes[run1_clusterid]:
        if gene in run2_geneToCluster:
            return True
    return False       

# Get three strings:
# - a list of genes overlapping between 
def getComparisonStrings(run1_clusterid, run2_clusterid, run1_clusterToGenes, run2_clusterToGenes, run1_geneToCluster, run2_geneToCluster, includeall):
    genes1 = set(run1_clusterToGenes[run1_clusterid])
    genes2 = set(run2_clusterToGenes[run2_clusterid])
    allgenes1 = set(run1_geneToCluster.keys())
    allgenes2 = set(run2_geneToCluster.keys())
    overlap = ";".join( list(genes1 & genes2) )
    if includeall:
        run1_only = ";".join( list(genes1 - genes2) )
        run2_only = ';'.join( list(genes2 - genes1) )
    else:
        run1_only = ";".join( list( (genes1 - genes2) & allgenes2 ) )
        run2_only = ";".join( list( (genes2 - genes1) & allgenes1 ) )
    return overlap, run1_only, run2_only


##############
# Main script
##############

con = sqlite3.connect("db/methanosarcina")
cur = con.cursor()

run1_geneToCluster, run1_clusterToGenes = getClusterRunDictionaries(args[0], cur)
run2_geneToCluster, run2_clusterToGenes = getClusterRunDictionaries(args[1], cur)

# This will be accumulated to avoid printing duplicate comparisons between clusters
donePairs = []
for cluster1 in run1_clusterToGenes:
    if not overlapExists(cluster1, run1_clusterToGenes, run2_geneToCluster):
        if options.includeall:
            overlap = ""
            run1only = ";".join(run1_clusterToGenes[cluster1])
            run2only = ""
            cluster2 = "NO_OVERLAP"
            print "%s\t%s\t%s\t%s\t%s\t%s\t%s" %(args[0], cluster1, args[1], cluster2, overlap, run1only, run2only)
            continue
        else:
            # Just skip these
            continue
    # Now we know there is at least some overlapping gene. Lets find out what that is.
    for gene in run1_clusterToGenes[cluster1]:
        # If it's not in there, we know there has to be SOMETHING overlapping so we don't worry about
        # missing any comparisons. We'll just check if we want to add them to the list of non-overlaps
        # for valid comparisons.
        if gene in run2_geneToCluster:
            cluster2 = run2_geneToCluster[gene]
            if (cluster1, cluster2) in donePairs or (cluster2, cluster1) in donePairs:
                continue
            overlap,run1only,run2only = getComparisonStrings(cluster1, cluster2, run1_clusterToGenes, run2_clusterToGenes, run1_geneToCluster, run2_geneToCluster, options.includeall)
            print "%s\t%s\t%s\t%s\t%s\t%s\t%s" %(args[0], cluster1, args[1], cluster2, overlap, run1only, run2only)
            donePairs.append( (cluster1, cluster2) )

# Do the exact same thing but the opposite direction (since we want a symmetric compare)
# This time we only need to look at things that have no overlap - everything else should've already been done above.
for cluster2 in run2_clusterToGenes:
    if not overlapExists(cluster2, run2_clusterToGenes, run1_geneToCluster):
        if options.includeall:
            overlap = ""
            run2only = ";".join(run2_clusterToGenes[cluster2])
            run1only = ""
            cluster1 = "NO_OVERLAP"
            print "%s\t%s\t%s\t%s\t%s\t%s\t%s" %(args[0], cluster1, args[1], cluster2, overlap, run1only, run2only)
            continue
        else:
            # Just skip these
            continue

con.close()
