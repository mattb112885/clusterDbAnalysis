#!/bin/sh

######### Safety ###############################################################################
#   Comment this part out to commit to removing the files printed out by the safe version!!   #
echo "WARNING: Comment the exit 0 line below out if you actually want to use this tool."
echo "Removing an organism is irreversible!"
echo "I highly suggest running removeOrganism_SAFE first which will print out names of files to remove."
echo "Make sure those files are right before running this script!!!"
#exit 0;

###############################################################################################

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

# Remove organism from organisms file
echo "Removing line from organisms file: "
grep -v -w ${ORGANISM} organisms > organisms_mod;
mv organisms_mod organisms

# Try to CD into each of the data folders and remove the organism's data
# Be careful of things like xxxx.1 vs xxxx.10
if [ -d ./blastres/ ]; then
    echo "Blastres folder";
    # Beginning of filename
    cd ./blastres/;
    ls | grep -P "^${ORGANISM}\.txt" | xargs rm;
    # Middle of filename
    ls | grep -P "_${ORGANISM}\.txt" | xargs rm;
    cd ..
fi

if [ -d ./clusters/ ]; then
    echo "Clusters folder";
    cd ./clusters/;
    ls | xargs rm;
    cd ..;
fi

if [ -d ./db/ ]; then
    echo "Database folder";
    cd ./db/;
    ls | xargs rm;
    cd ..;
fi

if [ -d ./faa/ ]; then
    echo "AA fasta folder";
    cd ./faa/;
    ls | grep -P "^${ORGANISM}\.txt" | xargs rm;
    cd ..;
fi

if [ -d ./fna/ ]; then
    echo "NT fasta folder";
    cd ./fna/;
    ls | grep -P "^${ORGANISM}\.txt" | xargs rm;
    cd ..;
fi

if [ -d ./genbank/ ]; then
    echo "Genbank folder";
    cd ./genbank;
    ls | grep -P "^${ORGANISM}\.gbk" | xargs rm;
    cd ..;
fi

if [ -d ./genomefasta/ ]; then
    echo "genome fasta folder";
    cd ./genomefasta/;
    ls | grep -P "^${ORGANISM}\.gbk" | xargs rm;
    cd ..;
fi

if [ -d ./modtable/ ]; then
    echo "Modified table folder";
    cd ./modtable/;
    ls | grep -P "^${ORGANISM}\.txt" | xargs rm;
    cd ..;
fi

if [ -d ./raw/ ]; then
    echo "RAW data table folder";
    cd ./raw/;
    ls | grep -P "^${ORGANISM}\.txt" | xargs rm;
    cd ..;
fi

if [ -d ./tblastn/ ]; then
    echo "Tblastn folder";
    cd ./tblastn/;
    # Organism at the beginning of the filename
    ls | grep -P "^${ORGANISM}\.gbk" | xargs rm;
    # Organism at the middle of the filename
    ls | grep -P "_${ORGANISM}\.txt" | xargs rm;
    cd ..;
fi

if [ -d ./tblasn_cutoff ]; then
    echo "Tblastn_cutoff folder";
    cd ./tblastn_cutoff;
    # Organism at the beginning of the filename
    ls | grep -P "^${ORGANISM}\.gbk" | xargs rm;
    # Organism at the middle of the filename
    ls | grep -P "_${ORGANISM}\.txt" | xargs rm;
    cd ..;
fi
