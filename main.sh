#!/bin/bash

if [ $# -lt 1 ]; then
    echo "usage: ./main.sh NCORES [BLASTP_cutoff] [BLASTN_cutoff]"
    echo ""
    echo "Description: runs BLAST, BLASTN all vs. all against all of the"
    echo "maintained conserved domain databases in NCBI"
    echo ""
    echo "This script MUST be run from the directory that contains it (i.e. you must call it like "
    echo "./main.sh NCORES"
    echo " or else it will not work."
    echo ""
    echo "The default E-value cutoff is 1E-5 for BLASTP and 1 for BLASTN. Use an argument for the cutoffs"
    echo "if you want to use something different. However, BE AWARE that if you change the cutoff and "
    echo "BLAST results already exist they will NOT be overwritten."
    echo ""
    exit 1
fi

./checkForDependencies.sh
if [ $? -eq 1 ]; then
    echo "ERROR: One or more required dependencies was not found";
    exit 1;
fi

NCORES=$1;

BLASTP_EVALUE="1E-5"
BLASTN_EVALUE="1"
if [ $# -ge 2 ]; then
    BLASTP_EVALUE="$2";
fi

if [ $# -ge 3 ]; then
    BLASTN_EVALUE="$3";
fi

# Make sure we are pointing at the correct folder for organisms, groups, etc
source SourceMe.sh

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

# If we have existing abbreviations we want to keep them around in the new version
# of the organisms file.
if [ -f organisms ]; then
    TM=$(date +%s)
    NEWNAME="organisms.${TM}.bk"
    mv organisms "${NEWNAME}"
    echo "WARNING: The old organism file has been backed up to ${NEWNAME}."
    echo "Automatically generating a new organisms file..."
    ./generateOrganismFileFromGbk.sh
    echo "WARNING: The new organism file has consistent abbreviations with previous one"
    echo "but any organisms not in genbank or raw folders have been removed"
#    rm "${NEWNAME}"
else
   # Automatically generate a organism file and a default groups file containing all the organisms
   # Since no previously-existing organism file is available I just assign default (unique) abbreviations
   # for each organism. 
   ./generateOrganismFileFromGbk.sh
fi

# Remove the existing "all" group from the groups file if it exists
# Any other groups are ignored in that file.
if [ -f "groups" ]; then
    cat "groups" | grep -v -P "^all\s" > "groupmod"
    mv "groupmod" "groups"
fi

./addGroupByMatch.py -n all

./checkInputFormat.sh
if [ $? -ne 0 ]; then
    echo "Errors in input file format (see printout for details). Please fix and run again."
    exit 1
fi

# In case the user has never run this script before...
mkdir faa 2> /dev/null;
mkdir fna 2> /dev/null;
mkdir modtable 2> /dev/null;
mkdir blastres 2> /dev/null;
mkdir blastn_res 2> /dev/null;
mkdir db 2> /dev/null;
mkdir db/orthofasta 2> /dev/null;
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

echo "Concatinating faa and fna files..."
echo faa/*.faa | xargs cat > db/allgenomes.faa;
echo fna/*.fna | xargs cat > db/allgenomes.fna;

# Haven't converted this one to a pipe yet...
# Run blast all vs. all organisms (takes ~ 1 hour for 18 organisms and 8 cores)
echo "Blast all vs all (BLASTP)...";
Blast_all_v_all.py -e "${BLASTP_EVALUE}" faa/ blastres/ ${NCORES};

echo "Blast all vs all (BLASTN)..."
Blast_all_v_all.py -n -e "${BLASTN_EVALUE}" fna/ blastn_res/ ${NCORES};

# Concatinate results files for input into the database
# Strip titles to avoid whining about duplicate rows.
# These tables contain ALL of the organisms. The strategy is to keep everything together as much as possible
# and only to separate out organisms when needed.
#
# Also - add the gene aliases here (not above) because the annotations are ultimately built from the raw file.
echo "Concatinating BLAST results..."
# Need to use xargs here because of too-long strings of arguments for hundreds of genomes.
echo blastres/* | xargs cat > db/blastres_cat;
echo blastn_res/* | xargs cat > db/blastnres_cat;
echo modtable/* | xargs cat > db/mod_cat;

echo "Adding aliases from the alias file to the gene annotations..."
if [ -f aliases/aliases ]; then
    ls raw/* | grep -v "README" | xargs cat | grep -v "feature_id" | addAliasesToGeneAnnotations.py aliases/aliases > db/raw_cat;
else
    echo "WARNING: No alias file found (expected location: aliases/aliases). No aliases will be added to annotations."
    ls raw/* | grep -v "README" | xargs cat | grep -v "feature_id" > db/raw_cat;
fi

echo "Computing gene neighborhoods (up to a maximum of 10)..."
cat db/raw_cat | getNeighbors_bothStrands_rast.py  > db/neighborhoods

# Generate the first part of the SQL database
# (for use with generating the clusters for a specific group of organisms)
# Note this does NOT use the FileLocator.py - which is GOOD becuase
# that means we can rebuild while keeping FileLocator.py pointing at a
# temporary copy so that people can still query it until the new copy is built.
echo "Rebuilding database...";
rm db/DATABASE.sqlite 2> /dev/null;
sqlite3 db/DATABASE.sqlite < src/internal/builddb_1.sql;