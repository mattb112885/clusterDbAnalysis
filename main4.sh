#!/bin/bash

if [ $# -ne 1 ]; then
    echo " USAGE: ./main4.sh [NCORES]"
    echo ""
    echo " This script downloads the NCBI CDD database (if it isn't already"
    echo " on your machine) and runs RPS-BLAST to identify similarity of all of your"
    echo " proteins to existing conserved domains from CD, COG, TIGR, ..."
    echo ""
    echo "This script must be run from the directory in which it is located"
    echo ""
    exit 1
fi

NCORES=$1

# Make sure we are pointing at the right version of the repo
source SourceMe.sh

mkdir rpsblast_res 2> /dev/null
mkdir cd_db 2> /dev/null

cd cd_db
if [ ! -f Cdd.pn ]; then
    echo "DOWNLOADING the conserved protein domain database files..."
    wget ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/cdd.tar.gz
    # Note - this tar.gz file doesn't create another directory for all the files.
    echo "Unzipping the conserved protein domain database..."
    tar xzf cdd.tar.gz
    rm cdd.tar.gz
else
    echo "Conserved-protein database already present in expected location."
fi

# Gunzip (unlike tar xzf) deletes the gz file so I check for the file inside that.
if [ ! -f ../db/cddid.tbl ]; then
    echo "DOWNLOADING conserved protein database info file..."
    wget ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/cddid.tbl.gz
    gunzip cddid.tbl.gz
    # I finally pinned down the issue with sqlite loading the cddid table to
    # the presence of quotation marks in the annotation.
    # SQlite interprets them as enclosing a single item to read and therefore thinks
    # that the number of columns is incorrect...
    sed s/\"//g cddid.tbl > ../db/cddid.tbl
fi

# This program has to be run from the cd_db folder - it doesn't work if you try to compile
# while you're in a different folder
#
# Note - the scale 1.0 was "strongly recommended" in one of the NCBI docs. Without it, the
# sensitivity is drastically decreased, but with it we find everything that we expect from e.g.
# searching proteins vs. the CDD website.
if [ ! -f Cdd.pn.aux ]; then
    echo "Compiling the CDD RPSBLAST database..."
    makeprofiledb -in Cdd.pn -dbtype rps -scale 1.0
fi

cd ..

Rpsblast_all_vs_one.py -n "${NCORES}" -c 1E-5 "cd_db/Cdd.pn" "faa/" "rpsblast_res/"
cat rpsblast_res/* > db/external_CDD
# Note - it's possible that this will change in the future (if NCBI decides to change the
# output file formats again). But for now all my hits look like gnl|CDD|...
#
# The ID number resulting from stripping that off SHOULD match up with the first column
# of the cddid.tbl table. So we can salvage this...
cat db/external_CDD | sed -r 's/gnl\|CDD\|([0-9]+)/\1/g' > db/external_CDD_MOD
mv db/external_CDD_MOD db/external_CDD

sqlite3 db/DATABASE.sqlite < src/internal/builddb_4.sql