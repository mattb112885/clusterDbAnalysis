#!/usr/bin/python

#
# Creates a cluster comparison file based on all of the cluster runs present in
# the database
#
# The comparison file looks like this:
#
# [Basis_gene] [RunID_1] [ClusterID_1] [RunID_2] [ClusterID_2] [+/-/SAME] [GeneID]
#
# Each line is a comparison of pair (RunID_1, ClusterID_1) and (RunID_2, ClusterID_2)
# that share gene [Basis_gene] in common
#
# IF two clusters are identical there will be a line printed that looks like this:
# [RunID_1] [ClusterID_1] [RunID_2] [ClusterID_2] SAME [NULL]
#
# Otherwise, for each gene that is different between the two clusters
# there will be a line like this
# [Basis_gene] [RunID_1] [ClusterID_1] [RunID_2] [ClusterID_2] [+/-] [GeneID]
# + means the gene is found in set 1 but not in set 2,
# - means the gene is found in set 2 but not in set 1.

import sqlite3

def getGenesInCluster(runid, clusterid, cursor):
    cursor.execute("SELECT * FROM clusters WHERE clusters.runid = ? AND clusters.clusterid = ?;", 
                   (runid, clusterid))
    genelist = []
    for l in cursor:
        s = list(l)
        genelist.append(s[2])
    return genelist

def getAllGenes(cursor):
    cursor.execute("SELECT geneid FROM processed;")
    genelist = []
    for l in cursor:
        s = list(l)
        genelist.append(s[0])
    return genelist

def getUniqueClustersContainingGene(geneid, cursor):
    cursor.execute("SELECT DISTINCT runid, clusterid FROM clusters WHERE clusters.geneid = ?", (geneid, ) )
    clusterList = []
    for l in cursor:
        clusterList.append(list(l))
    return clusterList

con = sqlite3.connect("db/methanosarcina")
cur = con.cursor()

genelist = getAllGenes(cur)

for gene in genelist:
    clustList = getUniqueClustersContainingGene(gene, cur)
    # Iterate over the list twice to get all pairwise comparisons between clusters.
    for i in range(0, len(clustList)):
        for j in range(i, len(clustList)):
            if clustList[i][0] == clustList[j][0] and clustList[i][1] == clustList[j][1]:
                continue
            genelist1 = getGenesInCluster(clustList[i][0], clustList[i][1], cur)
            genelist2 = getGenesInCluster(clustList[j][0], clustList[j][1], cur)
            # The meat of the code: are the two lists the same?
            isSame = True
            # "+" : Gene in list 1 but not in list 2
            for g in genelist1:
                if not g in genelist2:
                    isSame = False
                    print "%s\t%s\t%s\t%s\t%s\t%s\t%s" %( gene, clustList[i][0], clustList[i][1], clustList[j][0], clustList[j][1], "+", g )
            # "-" : Gene in list 2 but not in list 1
            for g in genelist2:
                if not g in genelist1:
                    isSame = False
                    print "%s\t%s\t%s\t%s\t%s\t%s\t%s" %( gene, clustList[i][0], clustList[i][1], clustList[j][0], clustList[j][1], "-", g )

            # SAME: The two clusters are identical
            if isSame:
                print "%s\t%s\t%s\t%s\t%s\t%s\t%s" %( gene, clustList[i][0], clustList[i][1], clustList[j][0], clustList[j][1], "SAME", "" )
