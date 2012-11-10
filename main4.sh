#!/bin/bash

if [ $# -ne 1 ]; then
    echo " USAGE: main4.sh [NCORES]"
    echo ""
    echo " This script downloads the NCBI CDD database (if it isn't already"
    echo " on your machine) and runs RPS-BLAST to identify similarity of all of your"
    echo " proteins to existing conserved domains from CD, COG, TIGR, ..."
    echo ""
    exit 1
fi

NCORES=$1

if [ ! -f db/external_CDD ]; then
    ./computeConservedDomains.sh "CDD" db/allgenomes.faa db/external_CDD "${NCORES}";
    cat db/external_CDD | sed -r "s/^(.*?\s+.*?),/\1/g" > db/external_CDD_MOD
    mv db/external_CDD_MOD db/external_CDD
fi

sqlite3 db/DATABASE.sqlite < src/builddb_4.sql