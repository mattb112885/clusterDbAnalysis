#!/bin/bash
#
# convertKbaseToRast.py
#
# Call KBase functions to obtain a "raw" file for a paritcular
# organism so that it can be imported to the database... 
#
# KBase functions can be gotten from the KBase website
# http://kbase.science.energy.gov/developer-zone/downloads/
#
# Eventually they will also be usable from CPAN.

if [ $# -ne 1 ]; then
    echo "Usage: convertKbaseToRast.py (KBase genome ID)";
    echo ""
    echo "Description: Given a KBase genome ID, converts it into"
    echo "the RAW format needed for input into the database."
    echo "It will attempt to find a SEED ID for the genome if one"
    echo "exists, and otherwise take the taxID (which must be available)"
    echo "and append .99999 ."
    echo "At the moment only SEED and MOL genomes are acceptable."
    exit 0;
fi

ORGANISM="$1"

TMPFILE="${ORGANISM}".TMP

# For now we ONLY support KBase IDs that come from either microbes online (MOL)
# or from the SEED.
#
# The reason for this is because all of the IDs in those databases are on a common
# framework (based on NCBI taxIDs), and are therefore easy to convert without worrying as
# much about namespace collisions.
#
LINE=$(echo "${ORGANISM}" | kbase_ids_to_external_ids);

SOURCE=$(echo "${LINE}" | cut -f 2);
ID=$(echo "${LINE}" | cut -f 3);
if [ "${SOURCE}" = "SEED" ]; then
    echo "Used existing SEED genome ID ${ID}";
elif [ "${SOURCE}" = "MOL:Genome" ]; then
    # The MOL ID is the tax ID by itself without a revision number.
    # I use revision number 99999 - unlikely (though not impossible) to cause conflicts with existing SEED genomes...
    echo "Original MOL genome ID ${ID}";
    ID="${ID}.99999";
    echo "Assigned MOL genome the ID ${ID}";
else
    echo "ERROR: Specified organism ${ORGANISM} either not found in the KBASE or"
    echo "was not from one of the acceptable genome sources (SEED and MOL)"
    exit 1;
fi

# Now that we've got some IDs we like, lets pull out some data.
OUTFILE="${ID}.txt"

# The final 2> /dev/null is to get rid of things with no protein sequences (which the DB
# doesn't support anyway)
echo "${ORGANISM}" | get_relationship_IsComposedOf -to id | \
    get_relationship_IsLocusFor -rel begin,len,dir -to id,function,feature_type | \
    fids_to_dna_sequences -c 6 -fasta=0 | \
    fids_to_protein_sequences -c 6 -fasta=0 2> /dev/null | kbase2rast.py "${ID}" > "${OUTFILE}"