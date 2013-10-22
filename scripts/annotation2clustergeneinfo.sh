#!/bin/sh

if [ $# -ne 2 ]; then
    echo ""
    echo "USAGE: ./annotation2clustergeneinfo.sh annotation run_ID > geneinfo_file"
    echo ""
    echo "DESCRIPTION: Given a gene annotation,alias, or anything that is found in the"
    echo "annotation (NOT a gene ID as given in the database), e.g. 'MA0001',"
    echo "generates a geneinfo file for the cluster(s) containing matching genes"
    echo ""
    echo "SEE ALSO ./geneid2clustergeneinfo.sh"
    echo ""
    echo "DEPENDENCIES: The database must be loaded (with main.sh) and the specified run ID"
    echo "must be one of the possible run IDs (from db_getAllClusterRuns.py)."
    echo ""
    exit 1;
fi

GENEALIAS="$1"
RUNID="$2";

# First find genes matching the given alias...
# then get the clustering results (across all the cluster runs),
# then get a unique list of the clusters containing those genes (since the
# db_getClustersContainingGenes.py function will print one line for each gene,
# if there are multiple input genes in the same cluster the same cluster will
# appear multiple times), and finally get the cluster geneinfo.
#
# Note that this can still result in multiple clusters in the same geneinfo file,
# if multiple genes match the specified alias and they appear in different clusters
# in the specified cluster run...
OUTFILE="${GENEALIAS}_${RUNID}.geneinfo"
db_getGenesWithAnnotation.py "${GENEALIAS}" | db_getClustersContainingGenes.py | \
    grep -F "${RUNID}" | cut -f 1,2 | sort -u | db_getClusterGeneInformation.py > "${OUTFILE}"

if [ $( cat "${OUTFILE}" | cut -f 13,14 | sort -u | wc -l) -ne 1 ]; then
    echo "WARNING: The specified annotation fragment ${GENEALIAS} resulted in representatives of more than one cluster."
    echo "It may be necessary to look at the results of db_getGenesWithAnnotation and call geneid2clustergeneinfo.sh"
    echo "with the specific gene you are interested in."
else
    echo "The specified annotation gave a unique cluster."
fi