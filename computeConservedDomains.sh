#!/bin/bash

if [ $# -lt 4 ]; then
    echo "USAGE: computeConservedDomains.sh [which_family_group] [query_file] [outfile] [NUMCORES]"
    echo ""
    echo "DESCRIPTION: Automatically download, extract conserved domain profiles"
    echo "based on the NCBI data (ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/)"
    echo ""
    echo "Then compile an RPS-BLAST database with the specified group of protein families"
    echo "e.g. PFAM, COG, or CDD"
    echo ""
    echo "Finally, run RPS-BLAST on query_file and store the results in outfile"
    echo "Because RPSBLAST (unlike BLASTP) gives a LOT of hits to some clusters, we restrict it to"
    echo "only print out one of them."
    echo ""
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

# Gunzip (unlike tar xzf) deletes the gz file so I check for the file inside that.
if [ ! -f ../db/cddid.tbl ]; then
    echo "DOWNLOADING conserved protein database info file..."
    wget ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/cddid.tbl.gz
    gunzip cddid.tbl.gz
    mv cddid.tbl ../db/cddid.tbl
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
    # Note - NCBI recommends setting -S 1 although -S 100 is the default.
    # See ftp://ftp.ncbi.nih.gov/blast/documents/rpsblast.html
    # However that aws for the old version and the docs for the newer one don't say anything
    # about it...
    formatrpsdb -i "${MYPN}" -t "${WHICHDB}" -S 1;
else
    echo "RPSBLAST database already formatted."
fi

cd ..;

echo "Runing RPSBLAST..."
rpsblast -i "$2" -d "cd_db/${MYPN}" -o "$3" -m 8 -e 1E-5 -a "$4" -P 1;

echo "Done."