#!/bin/sh

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
    echo "Usage: convertPubseedToRast [genome ID]";
    exit 0;
fi

ORGANISM="$1"

# Add nucleic acid and protein sequences from pubseed
# based on the server scripts
# First call to svr_fasta adds the nucleotide sequence for the proteins
# and second call adds the protein sequence.

wget ftp://ftp.theseed.org/genomes/SEED/${ORGANISM}.tbl

cat "${ORGANISM}.tbl" | svr_fasta -c 1 | svr_fasta -c 1 -protein > "${ORGANISM}.int"

# This function will convert as best it can from the tabular fromat of the pubseed
# (with sequences added) into the RAST "raw" format.

FILESIZE=$(du ${ORGANISM}.int | cut -f 1)
if [ ${FILESIZE} -eq 0 ]; then
    echo "Error: No proteins available for organism ${ORGANISM}"
    exit 0;
fi

cat "${ORGANISM}.int" | pubseed2rast.py > "${ORGANISM}.rast_txt"

# CLean up
rm "${ORGANISM}.tbl" "${ORGANISM}.int"