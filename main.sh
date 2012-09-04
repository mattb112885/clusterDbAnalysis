#!/bin/bash

NCORES=12;

# Note - for now the mcl parameters are located in the function
# src/specificOrganismClusterDriver.py

##########################

# Main function that generates blast results and calculates intermediate values,
# storing them in a SQLITE file
#
# This function must be run from the root install directory (so that I can find all the data files).
#
# Required input files:
# - Tab-delimited seed files (in "raw" directrory)
# - organism file (in root directory)
#
# Required to install:
#
# Python
# Ruffus (python package - for parallelizing)
# BLAST+ (including TBLASTN for pseudogene finding, and BLASTP for initial clustering...)
# SQLITE (for the database / storage / querying as needed)
#
# See README file for more details

./checkInputFormat.sh
if [ $? -ne 0 ]; then
    echo "Errors in input file format. Please fix and run again."
    exit 1
fi

# In case the user has never run this script before...
mkdir faa 2> /dev/null;
mkdir fna 2> /dev/null;
mkdir modtable 2> /dev/null;
mkdir blastres 2> /dev/null;
mkdir db 2> /dev/null;
mkdir aliases 2> /dev/null;

# Do some preprocessing required for inputs into various software tools
echo "Preprocessing raw files..."

cd raw;
for file in $(ls | grep -v "README"); do
    echo ${file}
    cat ${file} | raw2faa.py > ../faa/${file}.faa;
    cat ${file} | raw2fna.py > ../fna/${file}.fna;
    cat ${file} | raw2processed.py ../organisms > ../modtable/${file}.tsv;
done
cd ..;

# Haven't converted this one to a pipe yet...
# Run blast all vs. all organisms (takes ~ 1 hour for 18 organisms and 8 cores)
echo "Blast all vs all...";
Blast_all_v_all.py faa/ blastres/ ${NCORES};

# Concatinate results files for input into the database
# Strip titles to avoid whining about duplicate rows.
# These tables contain ALL of the organisms. The strategy is to keep everything together as much as possible
# and only to separate out organisms when needed.
#
# Also - add the gene aliases here (not above) because the annotations are ultimately built from the raw file.
cat blastres/* > db/blastres_cat;
cat modtable/* > db/mod_cat;

if [ -f aliases/aliases ]; then
    ls raw/* | grep -v "README" | xargs cat | grep -v "feature_id" | addAliasesToGeneAnnotations.py aliases/aliases > db/raw_cat;
else
    echo "WARNING: No alias file found (expected location: aliases/aliases). No aliases will be added to annotations."
    ls raw/* | grep -v "README" | xargs cat | grep -v "feature_id" > db/raw_cat;
fi

cat db/raw_cat | getNeighbors_bothStrands_rast.py  > db/neighborhoods

# Generate the first part of the SQL database
# (for use with generating the clusters for a specific group of organisms)
echo "Rebuilding database...";
rm db/DATABASE.sqlite 2> /dev/null;
sqlite3 db/DATABASE.sqlite < src/builddb_1.sql;