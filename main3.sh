#!/bin/bash

NCORES=8

# Second part of the pipeline:
#
# We have clusters, so now how do we identify the pseudogenes?
#  / curate the clusters
# The answer, tentatively, is with TBLASTN
#

# Requirements:
# Biopython (for parsing the genbank files - this can be skipped if you have
#      nucleotide fasta files already)
# BLAST+ (tblastn)

mkdir tblastn 2> /dev/null;
mkdir tblastn_cutoff 2> /dev/null;
mkdir contig_fasta 2> /dev/null;

echo "Converting genbank files to full-nucleotide fasta files...";

cd genbank;
for file in *; do
    orgid=$(echo ${file} | grep -P -o "\d+\.\d+")
    if [ $? -ne 0 ]; then
	echo "WARNING: No organism ID could be inferred from filename ${file} so it was not converted to a contig fasta file"
	continue;
    fi
    echo ${orgid}
    genbank2nucleotides.py -f ${file} -o $(echo ${file} | grep -P -o "\d+\.\d+") > ../contig_fasta/${file}.fna;
done
cd ..;

echo "Running tblastn...";
python src/tblastn_all_vs_all.py faa/ contig_fasta/ tblastn/ ${NCORES};

echo "tblastn done. processing (removing very weak hits)..."
cd tblastn;
for file in *; do
    cat ${file} | processTblastResults.py 50 > ../tblastn_cutoff/${file}.cut;
done
cd ..;

echo "Done. Concatinating tblastn results..."
cat tblastn_cutoff/* > db/tblastn_cat;

echo "Done. Adding TBLASTN information to the database... (this takes a long time)"
sqlite3 db/methanosarcina < src/builddb_3.sql