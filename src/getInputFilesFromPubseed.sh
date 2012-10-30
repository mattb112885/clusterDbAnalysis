#!/bin/bash

# Convert PUBSEED to RAST by calling appropriate server scripts, and then
# rearranging the tables to the right format.
#
# Output goes to ${ID}.txt in common format
# with the RAST files that are downloaded (and placed in the "RAW" folder)
#
# This must be run for any organism you wish to pull from pubseed for the analysis
# (as opposed to uploading from RAST)
#
# Example usage: convertToRast.sh 188937.1

if [ $# -ne 1 ]; then
    echo "Usage: convertPubseedToRast (PUBSEED genome ID)";
    echo ""
    echo "Description: Requires MYRast (svr_fasta) function"
    echo "Creates a RAW file (genomeid).rast_txt that can be used"
    echo "as input into the database."
    echo "See pubseed.theseed.org ."
    exit 0;
fi

ORGANISM="$1"

# Add nucleic acid and protein sequences from pubseed
# based on the server scripts
# First call to svr_fasta adds the nucleotide sequence for the proteins
# and second call adds the protein sequence.
URL="ftp://ftp.theseed.org/genomes/SEED/${ORGANISM}.tbl"
wget "${URL}"

if [ $? -ne 0 ]; then
    echo "ERROR: Organism table file not found on the SEED ftp server - are you sure this organism ID is correct?"
    echo "Expected location: ${URL}"
    exit 1
fi

# We need the genbank file too.
# These are VERY minimalist genbank files but they have all I need
# in them...
#
# There are bigger genbank files in the genbank.tmp directory -
# is it worth using them?
#
# Maybe if it kills the MyRAST dependency.
URL2="ftp://ftp.theseed.org/genomes/genbank/${ORGANISM}.gbk"
wget "${URL2}"

if [ $? -ne 0 ]; then
    echo "ERROR: Despite the organism table file existing for organism ${ORGANISM}, there was no genbank file."
    echo "Expected location: ${URL2}"
    echo "Unfortunately this means it is impossible to run the pipeline on this organism."
    rm "${ORGANISM}.tbl"
    exit 1
fi

# Add needed information (amino acid and nucleotide sequences) to the table file.
# This step requires the MyRAST function svr_fasta.
cat "${ORGANISM}.tbl" | svr_fasta -c 1 | svr_fasta -c 1 -protein > "${ORGANISM}.int"

# This function will convert as best it can from the tabular fromat of the pubseed
# (with sequences added) into the RAST "raw" format.
FILESIZE=$(du ${ORGANISM}.int | cut -f 1)
if [ ${FILESIZE} -eq 0 ]; then
    echo "ERROR: svr_fasta failed on organism ${ORGANISM} - no protein matches found"
    exit 1;
fi
cat "${ORGANISM}.int" | pubseed2rast.py > "${ORGANISM}.txt"

# Clean up
rm "${ORGANISM}.tbl" "${ORGANISM}.int"