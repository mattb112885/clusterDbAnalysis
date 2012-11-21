#!/bin/sh

if [ $# -ne 2 ]; then
    echo ""
    echo "USAGE: fasta2tree.sh [fasta-file] [outfile]"
    echo ""
    echo "DESCRIPTION: Given an un-aligned fasta file,"
    echo "Performs the following standard steps to make"
    echo "A newick file."
    echo ""
    echo "First, transforms the fasta file so that the ID"
    echo "for the gene becomes the existing ID + all of the"
    echo "annotations, removing all special characters."
    echo ""
    echo "Then runs MAFFT --auto and gblocks lenient."
    echo "Finally the script runs the FastTree wrapper with"
    echo "100 global bootstraps."
    echo ""
    echo "IF you need more flexibility with settings,"
    echo "use this script as a guide but replace it with the"
    echo "programs you want."
    echo ""
    exit 1;
fi

INFILE="$1"
OUTFILE="$2"
# Remove special characters from header lines
# (not from non-header lines because we dont want
# to accidentally add underscores to gene sequences -
# some programs output spaces as part of the sequence
# for some reason...)
cat "${INFILE}" | sed -r '/^>/s/[^>A-Za-z0-9]/_/g' > "TEMP1"
# Run MAFFT on the INFILE
mafft --auto "TEMP1" > "TEMP2"
# Use GBlocks (lenient) to trim poor-quality parts
# of the alignment
cat "TEMP2" | Gblocks_wrapper.py -r > "TEMP3"
# Run FastTREE with 100 bootstraps (requires
# phyml's seqboot program and the conversion
# script CompareToBootstrap.pl from the FastTree
# website)
cat "TEMP3" | FastTree_wrapper.py -g -b 100 > "${OUTFILE}"

rm "TEMP1" "TEMP2" "TEMP3"