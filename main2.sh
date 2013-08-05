#!/bin/bash

# Usage: ./main2.sh [inflation] [scoremethod] [cutoff]
# Note - for now the mcl parameters are located in the function
# src/specificOrganismClusterDriver.py
if [ $# -lt 3 ]; then
    echo "Usage: ./main2.sh [inflation] [scoremethod] [cutoff]";
    echo ""
    echo "Description: This script runs MCL with the specified parameters"
    echo "This script must be run from the directory where it is located"
    echo " (e.g. as ./main2.sh )"
    exit 0;
fi

# MCL parameter -I
INFLATION=$1;
# Method parameter -method
SCOREMETHOD=$2;
# MCL parameter -rcutoff
RCUTOFF=$3;

# Make sure we are pointing at the right version of the repo
source SourceMe.sh

##########################

# This is the second main function for database building.
#
# As input arguments, provide the scoring function, cutoff, and inflation parameter
# you wish to use for clustering. These are used to make an ID that will be different 
# for any different combination of these parameters that you wish to use. 
#
# This function must be run from the root install directory!
#
# Required input files:
# - Tab-delimited seed files (in "raw" directrory)
# - organism file (in root directory)
# - groups file (in root directory)
#
# Required to install:
#
# Python
# Ruffus (python package - for parallelizing)
# BLAST+ (including TBLASTN for pseudogene finding, and BLASTP for initial clustering...)
# MCL    (clustering based on blast)
# SQLITE (for the database / storage / querying as needed)
#
# See README file for more details

# In case the user has never run this script before...
mkdir clusters 2> /dev/null;
mkdir flatclusters 2> /dev/null;

# Check the groups file consistency
db_checkGroupsFile.py;
if [ $? -ne 0 ]; then
    echo "Errors in groups file. Stop."
    exit 1
fi

# Generate clusters based on specified groups
# (if the cluster file for a particular group already exists, it is skipped over
# since the results should generally be the same.)
echo "Running MCL on organisms in the groups file...";
db_specificOrganismClusterDriver.py groups ${INFLATION} ${RCUTOFF} ${SCOREMETHOD};

# Modify the format of cluster files to be SQL-friendly.
# Each separate clustering is given a random (very highly likely unique) identifier
# which is subsequently used to identify which clustering solution we want from SQL
# given a set of organism names
cd clusters;
for file in *; do
    cat ${file} | flattenClusterFile.py -n ${file} > ../flatclusters/${file};
done
cd ..;

cat flatclusters/* > db/flat_clusters;

# Import clusters into the db
echo "Importing cluster information into database...";
sqlite3 db/DATABASE.sqlite < src/internal/builddb_2.sql;

# Once clusters are loaded, make a pre-built presence\absence table.
# This is a python script because it involves some operations that would be impossible or hard in pure SQL.
echo "Making a pre-built presence\absence table in the database..."
db_loadPresenceAbsence.py;