#!/bin/sh

if [ $# -ne 4 ]; then
    echo ""
    echo "USAGE: setUpCytoscapeTable.sh [runid] [clusterid] [score_metric] [cutoff] > table"
    echo ""
    echo "DESCRIPTION: Create a BLAST score table for import into Cytoscape"
    echo "from a specified cluster run using a specified scoring metric"
    exit 1;
fi

RUNID="$1";
CLUSTERID="$2";
METRIC="$3";
CUTOFF="$4"

makeTabDelimitedRow.py "$RUNID" "$CLUSTERID" | \
    db_getGenesInClusters.py | \
    db_getBlastResultsBetweenSpecificGenes.py -g 3 | \
    makeBlastScoreTable.py -m "$METRIC" -c "$CUTOFF" -n | \
    db_addOrganismNameToTable.py -g 1 -a | \
    db_addOrganismNameToTable.py -g 2 -a

