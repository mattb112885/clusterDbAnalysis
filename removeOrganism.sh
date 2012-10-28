#!/bin/sh

###############################################################################################

if [ $# -lt 1 ]; then
    echo ""
    echo "Usage: removeOrganism.sh [organismId] [DELETEFLAG]";
    echo ""
    echo "DESCRIPTION: When just called with organismId (no DELETEFLAG), just prints a list of files or parts of files to be deleted"
    echo "To actually delete: run ./removeOrganism.sh [organismid] TRUE . This action is IRREVERSIBLE."
    echo "The TRUE is case sensitive and nothing else will work"
    echo ""
    echo "NOTE - for this script to work the files MUST file the standard naming conventions."
    echo "If the raw files and genbank files passed the test in checkInputFormat.sh then it should be OK."
    echo ""
    exit 0;
fi

DELETE=0
if [ $# -gt 1 ]; then
    if [ "$2" = "TRUE" ]; then
	echo "Deleting and modifying files..."
	DELETE=1
    fi
fi
if [ $DELETE -eq 0 ]; then
    echo ""
    echo "This function will only list all files that will be deleted. Pass TRUE as the DELETEFLAG to actually perform these deletions"
    echo ""
fi

ORGANISM=$1;

S=$(echo "${ORGANISM}" | grep -P -w "\d+\.\d+");
if [ $? -eq 1 ]; then
    echo "ERROR: Specified organism ID ${ORGANISM} is invalid. Stop."
    echo ""
    exit 0;
fi

# Test if the organism exists as a whole unit
if grep -q -w -F "${ORGANISM}" organisms; then
    echo "Organism found!";
    echo ""
else
    echo "Organism ${ORGANISM} Not found in the database. Stop.";
    echo ""
    exit 0;
fi

ORGNAME=$(cat organisms | grep -F -w "${ORGANISM}" | cut -f 1);

# Remove organism from organisms file
if [ $DELETE -eq 0 ]; then
    echo "Lines that would be removed from organisms file:"
    cat organisms | grep -F -w "${ORGANISM}"
    echo ""
else
    grep -F -v -w "${ORGANISM}" organisms > organisms_mod;
    mv organisms_mod organisms
fi

# Remove any groups including that organism from the groups file.
if [ $DELETE -eq 0 ]; then
    # Since the groups file now explicitly lists all the organism names we can do this.
    echo "Lines that would be removed from GROUPS file:"
    cat groups | grep -F "${ORGNAME}"
    echo ""
else
    cat groups | grep -v -F "${ORGNAME}" > groups_mod
    mv groups_mod groups
fi

# Remove aliases involving that organism from the aliases file.
if [ -f ./aliases/aliases ]; then
    if [ $DELETE -eq 0 ]; then
	echo "Lines that would be removed from the aliases file"
        # fig|organism.peg.genenum
	grep -F "|${ORGANISM}." ./aliases/aliases
	echo ""
    else
	echo "Removing lines from alias file..."
	grep -v -F "|${ORGANISM}." ./aliases/aliases > ./aliases/aliases_mod
	mv ./aliases/aliases_mod ./aliases/aliases
    fi
fi

# Remove organism-specific blast results,
# being careful not to let 83333.1 be the same as 83333.10
if [ -d ./blastres/ ]; then
    echo "blastres folder";
    cd ./blastres/;
    if [ $DELETE -eq 0 ]; then
	echo "Files that would be deleted from BLASTRES folder"
        # Beginning of filename (ORGANISM is the query)
	ls | grep -P "^${ORGANISM}\.txt";
        # Middle of filename (ORGANISM is the target)
	ls | grep -P "_${ORGANISM}\.txt";
	echo ""
    else
	ls | grep -P "^${ORGANISM}\.txt" | xargs rm;
	ls | grep -P "_${ORGANISM}\.txt" | xargs rm;
    fi
    cd ..
fi

if [ -d ./blastn_res/ ]; then
    echo "blastnres folder";
    cd ./blastn_res/;
    if [ $DELETE -eq 0 ]; then
	echo "Files that would be deleted from BLASTNRES folder"
        # Beginning of filename (ORGANISM is the query)
	ls | grep -P "^${ORGANISM}\.txt";
        # Middle of filename (ORGANISM is the target)
	ls | grep -P "_${ORGANISM}\.txt";
	echo ""
    else
	ls | grep -P "^${ORGANISM}\.txt" | xargs rm;
	ls | grep -P "_${ORGANISM}\.txt" | xargs rm;
    fi
    cd ..
fi

if [ -d ./clusters/ ]; then
    echo "Clusters folder";
    cd ./clusters/;
    if [ $DELETE -eq 0 ]; then
	echo "All clusters would be deleted (re-run main2.sh with the same parameters to get clusters without the deleted organism)"
	echo "List of files that would be deleted:"
	# Just in case we decide we want to put a README here after all...
	ls | grep -v "README"
	echo ""
    else
	ls | grep -v "README" | xargs rm;
    fi
    cd ..;
fi

if [ -d ./flatclusters/ ]; then
    echo "Flat Clusters folder";    
    cd ./flatclusters/;
    if [ $DELETE -eq 0 ]; then
	echo "All clusters would be deleted (re-run main2.sh with the same parameters to get clusters without the deleted organism)"
	echo "List of files to be deleted:"
	ls | grep -v "README"
	echo ""
    else
	ls | grep -v "README" | xargs rm;
    fi
    cd ..;
fi

if [ -d ./db/ ]; then
    echo "Database folder";
    cd ./db/;
    if [ $DELETE -eq 0 ]; then
	echo "All data in the database would be deleted (re-run main scripts to re-build without the deleted organism)"
	echo "List of files to be deleted:"
	ls | grep -v "README"
	echo ""
    else
	ls | grep -v "README" | xargs rm -r;
    fi
    cd ..;
fi

if [ -d ./faa/ ]; then
    echo "AA fasta folder";
    cd ./faa/;
    if [ $DELETE -eq 0 ]; then
	echo "List of files that would be deleted from the AA fasta folder:"
	ls | grep -P "^${ORGANISM}\.txt";
	echo ""
    else
	ls | grep -P "^${ORGANISM}\.txt" | xargs rm;
    fi
    cd ..;
fi

if [ -d ./fna/ ]; then
    echo "NT fasta folder";
    cd ./fna/;
    if [ $DELETE -eq 0 ]; then
	echo "List of files that would be deleted from the NT fasta folder:"
	ls | grep -P "^${ORGANISM}\.txt";
	echo ""
    else
	ls | grep -P "^${ORGANISM}\.txt" | xargs rm;
    fi
    cd ..;
fi

if [ -d ./genbank/ ]; then
    echo "Genbank folder";
    cd ./genbank/;
    if [ $DELETE -eq 0 ]; then
	echo "Genbank file that would be deleted:"
	ls | grep -P "^${ORGANISM}\.gbk";
	echo ""
    else
	ls | grep -P "^${ORGANISM}\.gbk" | xargs rm;
    fi
    cd ..;
fi

if [ -d ./modtable/ ]; then
    echo "Modified table folder";
    cd ./modtable/;
    if [ $DELETE -eq 0 ]; then
	echo "Modfiied table file that would be deleted:"
	ls | grep -P "^${ORGANISM}\.txt";
	echo ""
    else
	ls | grep -P "^${ORGANISM}\.txt" | xargs rm;
    fi
    cd ..;
fi

if [ -d ./raw/ ]; then
    echo "RAW data table folder";
    cd ./raw/;
    if [ $DELETE -eq 0 ]; then
	echo "RAW data table that would be deleted:"
	ls | grep -P "^${ORGANISM}\.txt";
	echo ""
    else
	ls | grep -P "^${ORGANISM}\.txt" | xargs rm;
    fi
    cd ..;
fi

