#!/bin/bash

if [ $# -lt 3 ]; then
    echo "USAGE: computeConservedDomains.sh [which_family_group] [query_file] [outfile]"
    echo ""
    echo "DESCRIPTION: Automatically download, extract conserved domain profiles"
    echo "based on the NCBI data (ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/)"
    echo ""
    echo "Then compile an RPS-BLAST database with the specified group of protein families"
    echo "e.g. PFAM, COG, or CDD"
    echo ""
    echo "This script won't do anything if the tar.gz file already exists."
    echo "To update the data, just delete the cd_db directory and re-run this script."
    exit 0
fi

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

WHICHDB=$(echo "$1" | tr [:lower:] [:upper:])

if [ "${WHICHDB}" == "CDD" ]; then
    MYPN="Cdd.pn"
elif [ "${WHICHDB}" == "COG" ]; then
    MYPN="Cog.pn";
elif [ "${WHICHDB}" == "KOG" ]; then
    MYPN="Kog.pn";
elif [ "${WHICHDB}" == "PFAM" ]; then
    MYPN="Pfam.pn";
elif [ "${WHICHDB}" == "PRK" ]; then
    MYPN="Prk.pn";
elif [ "${WHICHDB}" == "SMART" ]; then
    MYPN="Smart.pn";
elif [ "${WHICHDB}" == "TIGR" ]; then
    MYPN="Tigr.pn";
else
    echo "ERROR: Unsupported database ${WHICHDB}"
    exit 1
fi

if [ ! -f "${MYPN}.phr" ]; then
    echo "Formatting the RPSBLAST database..."
    formatrpsdb -i "${MYPN}" -t "${WHICHDB}";
else
    echo "RPSBLAST database already formatted."
fi

cd ..;
echo "Runing RPSBLAST..."
rpsblast -i "$2" -d "cd_db/${MYPN}" -o "$3" -m 8;

echo "Done."