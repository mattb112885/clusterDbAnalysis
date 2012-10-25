#!/bin/sh

if [ $# -lt 2 ]; then
    echo ""
    echo "USAGE: ./geneid2clustergeneinfo.sh [gene ID] [run ID] (alternate_base) > geneinfo_file"
    echo ""
    echo "DESCRIPTION: Given a gene ID as given in the database (e.g. fig|83333.1.peg.1, the fig| is optional)",
    echo "generates a geneinfo file for the cluster(s) containing matching genes"
    echo ""
    echo "Optionally: specify an alternate base name for the file (default: the gene ID)"
    echo ""
    echo "SEE ALSO ./annotation2clustergeneinfo.sh"
    echo ""
    echo "DEPENDENCIES: The database must be loaded (with main.sh) and the specified run ID"
    echo "must be one of the possible run IDs (from db_getAllClusterRuns.py). The gene ID must"
    echo "also be the gene for one of the genes in the database. Failure results in an empty file."
    echo ""
    exit 1;
fi

GENEID="$1"
RUNID="$2";

# Unlike the annotation2clustergeneinfo it should be impossible to get more than one
# cluster for a single gene... but I'll leave the warnings in anyway in case something
# unexpected happens.
if [ $# -eq 3 ]; then
    BASENAME="$3";
else
    BASENAME="${GENEID}";
fi
OUTFILE="${BASENAME}_${RUNID}.geneinfo"
echo "${GENEID}" | db_getClustersContainingGenes.py | \
    grep -F "${RUNID}" | cut -f 1,2 | sort -u | db_getClusterGeneInformation.py > "${OUTFILE}"

if [ $( cat "${OUTFILE}" | cut -f 13,14 | sort -u | wc -l) -ne 1 ]; then
    echo "WARNING: The specified gene ID ${GENEID} was in more than one cluster. THIS SHOULD BE IMPOSSIBLE WITH MCL"
else
    echo "Results printed to ${OUTFILE}"
fi