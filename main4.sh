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

mkdir rpsblast_res 2> /dev/null
mkdir cd_db 2> /dev/null

cd cd_db
if [ ! -f cdd.tar.gz ]; then
    echo "DOWNLOADING the conserved protein domain database files..."
    wget ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/cdd.tar.gz
    # Note - this tar.gz file doesn't create another directory for all the files.
    tar xzf cdd.tar.gz
else
    echo "Conserved-protein database already present in expected location."
fi

# Gunzip (unlike tar xzf) deletes the gz file so I check for the file inside that.
if [ ! -f ../db/cddid.tbl ]; then
    echo "DOWNLOADING conserved protein database info file..."
    wget ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/cddid.tbl.gz
    gunzip cddid.tbl.gz
    mv cddid.tbl ../db/cddid.tbl
fi

cd ..

if [ ! -f db/external_CDD ]; then
#        ./computeConservedDomains.sh "CDD" db/allgenomes.faa db/external_CDD "${NCORES}";
    Rpsblast_all_vs_one.py -n "${NCORES}" -c 1E-5 "cd_db/Cdd.pn" "faa/" "rpsblast_res/"
    cat rpsblast_res/* > db/external_CDD
    cat db/external_CDD | sed -r "s/^(.*?\s+.*?),/\1/g" > db/external_CDD_MOD
    mv db/external_CDD_MOD db/external_CDD
fi

exit 0;

sqlite3 db/DATABASE.sqlite < src/builddb_4.sql