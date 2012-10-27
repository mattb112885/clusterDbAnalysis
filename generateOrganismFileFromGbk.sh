#!/bin/bash

# The purpose of this function is to automatically generate a organism file
# from the (required) genbank file for each organism.
#
# The organism name is in a field that looks like this in both RAST and NCBI genbank files:
#
# /organism="[organism name]"
#
# And the taxon looks like this:
#
# /db_xref="Taxon:[taxid]"
#
# In RAST Genbank files the taxon ID is separated from the Taxon: field by a space, in NCBI it is not.
# The taxid should always be a number.
#
# I don't use biopython for this because RAST genbank files don't always play nicely with it.
# This is a "dumb parser".
#
# I make up an abbreviation as an incrementing number

rm organisms 2> /dev/null;

cd genbank;

ABBREV=0
for file in $(ls | grep -v "README"); do
    ORGLINE=$(cat "${file}"| grep -P -i "^\s*/organism\=\".*?\"")
    ORGANISM=$(echo "${ORGLINE}" | sed -r "s/.*?\"(.*?)\"/\1/g")
    # Do a sanity check on the genbank file name vs. its TaxID
    # The "extension" is rather arbitrary so I don't worry about that.
    # The file name must be [organismid].gbk which is [taxid].[extension].gbk
    TAXALINE=$(cat "${file}" | grep -P -i "^\s*/db_xref=\"Taxon\:\s*\d+\s*\"")
    TAXID=$(echo "${TAXALINE}" | sed -r "s/.*?\"[Tt]axon:[ ]*([0-9]+)[ ]*\"/\1/g")
    OK=$(echo "${file}" | grep -F "${TAXID}.")
    if [ $? -eq 1 ]; then
	echo "ERROR: The identified TaxID ${TAXID} from genbank file ${file} was not found in the name of that file."
	exit 0;
    fi
    # The actual organism ID will be extracted from the genbank file name since the file itself does not include the extension
    # i.e. the .1 in 83333.1
    ORGID=$(echo "${file}" | grep -o -P "\d+\.\d+")
    if [ $? -eq 1 ]; then
	echo "ERROR: The name of the genbank file must be [organismid].gbk - nothing with expected format of organism id (number.number) was found in the filename ${file}"
	exit 0;
    fi
    ABBREV=$(echo "${ABBREV} + 1" | bc)
    makeTabDelimitedRow.py "${ORGANISM}" "${ABBREV}" "${ORGID}" >> ../organisms
done

cd ..;


