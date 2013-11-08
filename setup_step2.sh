#!/bin/bash

# Usage: ./setup_step2.sh inflation scoremethod cutoff
# Note - for now the mcl parameters are located in the function
# src/specificOrganismClusterDriver.py
if [ $# -lt 3 ]; then
    echo "Usage: ./setup_step2.sh inflation scoremethod cutoff";
    echo ""
    echo "Description: This script runs MCL with the specified parameters"
    echo "This script must be run from the directory where it is located"
    echo " (e.g. as ./setup_step2.sh )"
    echo ""
    echo "For an alternative option for importing clustering results see"
    echo "importAllClusters.sh and importExternalClustering.py"
    echo ""
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

# Deal with importing clusters
sh importAllClusters.sh