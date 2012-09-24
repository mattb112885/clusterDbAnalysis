#!/bin/bash

# Given a bunch of genbank files in the genbank folder
# with names containing the organism ID in them,
# tmport the contigs into the sql database

cd genbank;
for file in *; do
    # Skip the readme file
    if [ "${file}" == "README" ]; then
	continue;
    fi
    orgid=$(echo "${file}" | grep -P -o "\d+\.\d+")
    if [ $? -ne 0 ]; then
	echo "WARNING: No organism ID could be inferred from filename ${file} so it was not converted to a contig fasta file"
	continue;
    fi
    echo ${orgid}
    genbank2nucleotides.py -f "${file}" -t -o "${orgid}" >> ../db/contigs;
done
cd ..;

sqlite3 db/DATABASE.sqlite < src/builddb_3.sql;