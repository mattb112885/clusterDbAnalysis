#!/bin/bash

if [ $# -ne 2 ]; then
    echo ""
    echo "USAGE: ./geneinfo2tree.sh geneinfo_file organism_file"
    echo ""
    echo "DESCRIPTION: Given a geneinfo file (from db_getGeneInformation.py or"
    echo "db_getClusterGeneInformation.py), uses MAFFT and FastTreeMP to compute"
    echo "a tree of all the proteins in that gene info file."
    echo "MAFFT is allowed to automatically pick the optimal algorithm for alignment"
    echo "while FastTree is used with -gamma and -wag"
    echo ""
    echo "DEPENDENCIES: MAFFT, FastTreeMP"
    echo ""
    exit 1;
fi

INFILE="$1";
ORGFILE="$2";
TMPFILE="/tmp/2802989302832.tmp";
cat "${INFILE}" | annoteSeq2Fasta.py -g 1 -a 10 -s 12 > "${TMPFILE}";
fasta2tree.sh "${TMPFILE}" "${INFILE}.nwk"
rm "${TMPFILE}";
