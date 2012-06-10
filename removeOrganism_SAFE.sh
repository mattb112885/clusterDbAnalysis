#!/bin/sh

# Need one input argument: The organism to remove
if [ $# -lt 1 ]; then
    echo "Usage: removeOrganism_SAFE.sh [organismId]";
    echo "Highly recommended to call this function before-hand to make sure all the files that will be removed are correct!"
    echo "WARNING: Removing the organism also will remove all cluster and database data (since that will all change)! Need to re-run main.sh to regenerate that data."
    exit 0;
fi

echo "This function will only list all files that will be deleted. Call the \"unsafe\" version to actually delete them"

ORGANISM=$1;

# Test if the organism exists as a whole unit
if grep -q -w ${ORGANISM} organisms; then
    echo "Organism found!";
else
    echo "Organism ${ORGANISM} Not found in the database.";
    exit 0;
fi

echo "Line to be removed from organisms file: "
grep -w ${ORGANISM} organisms;

# Try to CD into each of the data folders and remove the organism's data
# Be careful of things like xxxx.1 vs xxxx.10
if [ -d ./blastres/ ]; then
    echo "Blastres folder";
    # Beginning of filename
    cd ./blastres/;
    ls | grep -P "^${ORGANISM}\.txt";
    # Middle of filename
    ls | grep -P "_${ORGANISM}\.txt";
    cd ..
fi

if [ -d ./clusters/ ]; then
    echo "Clusters folder";
    cd ./clusters/;
    # All files would be removed
    ls;
    cd ..;
fi

if [ -d ./db/ ]; then
    echo "Database folder";
    cd ./db/;
    ls;
    cd ..;
fi

if [ -d ./faa/ ]; then
    echo "AA fasta folder";
    cd ./faa/;
    ls | grep -P "^${ORGANISM}\.txt";
    cd ..;
fi

if [ -d ./fna/ ]; then
    echo "NT fasta folder";
    cd ./fna/;
    ls | grep -P "^${ORGANISM}\.txt";
    cd ..;
fi

if [ -d ./genbank/ ]; then
    echo "Genbank folder";
    cd ./genbank;
    ls | grep -P "^${ORGANISM}\.gbk";
    cd ..;
fi

if [ -d ./genomefasta/ ]; then
    echo "genome fasta folder";
    cd ./genomefasta/;
    ls | grep -P "^${ORGANISM}\.gbk";
    cd ..;
fi

if [ -d ./modtable/ ]; then
    echo "Modified table folder";
    cd ./modtable/;
    ls | grep -P "^${ORGANISM}\.txt";
    cd ..;
fi

if [ -d ./raw/ ]; then
    echo "RAW data table folder";
    cd ./raw/;
    ls | grep -P "^${ORGANISM}\.txt";
    cd ..;
fi

if [ -d ./tblastn/ ]; then
    echo "Tblastn folder";
    cd ./tblastn/;
    # Organism at the beginning of the filename
    ls | grep -P "^${ORGANISM}\.gbk";
    # Organism at the middle of the filename
    ls | grep -P "_${ORGANISM}\.txt";
    cd ..;
fi

if [ -d ./tblasn_cutoff ]; then
    echo "Tblastn_cutoff folder";
    cd ./tblastn_cutoff;
    # Organism at the beginning of the filename
    ls | grep -P "^${ORGANISM}\.gbk";
    # Organism at the middle of the filename
    ls | grep -P "_${ORGANISM}\.txt";
    cd ..;
fi
