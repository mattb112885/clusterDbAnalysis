#!/bin/bash

# Check for existence and formatting of required input files.
# On success, will give no errors and retrun 0
#

STATUS=0

# Check existence of organism file. It will be automatically generated
# before calling this function but if something goes wrong this will alert us to it.
echo "Checking for existence of organisms file..."
if [ ! -f "organisms" ]; then
    echo "ERROR: organisms file not found in expected location";
    STATUS=1
fi

# Check that the organism ID is the third column
echo "Checking format of organism ID in organisms file..."
orgmatch=$(cat organisms | cut -f 3 | grep -P "^\d+\.\d+$")
if [ $? -eq 1 ]; then
    echo 'ERROR: Organism IDs not found or not in the expected format (third column of organisms file, and must have format #.# i.e. 83333.1 for E coli)';
    STATUS=1
fi

echo "Checking for uniqueness of organisms in organisms file..."
nl=$(cat organisms | wc -l)
numname=$(cat organisms | cut -f 1 | sort -u | wc -l);
numabbrev=$(cat organisms | cut -f 2 | sort -u | wc -l);
numid=$(cat organisms | cut -f 3 | sort -u | wc -l);
if [ "${nl}" != "${numname}" ]; then
    echo 'ERROR: Organism names must be unique (first column of organisms file)';
    STATUS=1;
fi
if [ "${nl}" != "${numabbrev}" ]; then
    echo 'ERROR: Organism abbreviations must be unique (second column of organisms file)';
    STATUS=1;
fi
if [ "${nl}" != "${numid}" ]; then
    echo 'ERROR: Organism IDs must be unique (third column of organisms file)';
    STATUS=1;
fi

# Check existence of groups file (for clustering)
echo "Checking existence of groups file..."
if [ ! -f "groups" ]; then
    echo 'ERROR: groups file not found in expected location';
    STATUS=1;
fi

# Check whether the organisms are found in the groups file
# (this isn't a fatal error if it fails)
echo "Checking whether organisms are found in the groups file..."
orgnames=$(cat organisms | cut -f 1);
for orgname in ${orgnames}; do
    ok=$(grep -F "${orgname}" groups)
    if [ $? -eq 1 ]; then
	echo "WARNING: Organism ${orgname} was in the organisms file but did not appear in any groups. Did you forget something?"
    fi
done

# For each orgmatch see if there is a raw file containing genes with that ID
cd raw;

echo "Checking the format of raw file names..."
for file in $(ls | grep -v "README"); do
    ok=$(echo "${file}" | grep -x -P "\d+\.\d+\.txt")
    if [ $? -eq 1 ]; then
	echo "ERROR: Raw file ${file} does not have a name in the expected format (organismid.txt)"
	STATUS=1
    fi
done

echo "Checking for existence of raw files for each organism..."
for org in ${orgmatch}; do
    # File name must exactly be [organismID].txt
    fmatch=$(ls | grep -w -F "${org}.txt");
    if [ $? -eq 1 ]; then
	echo "ERROR: No raw file match for organism ID ${org} - file name must be ${org}.txt";
	STATUS=1;
    fi
done

# Check that there arent any raw files that have no entry in the organism table
echo "Checking that all raw file organism IDs have an entry in the organisms file..."
for file in $(ls | grep -v "README"); do
    orgid=$(echo "${file}" | grep -o -P "\d+\.\d+");
    ok=$(cat ../organisms | cut -f 3 | grep -F -w "${orgid}")
    if [ $? -eq 1 ]; then
	echo "ERROR: Organism ID in raw file ${file} has no entry in the organism file"
	STATUS=1
    fi
done

echo "Checking formatting of each raw file..."
for file in $(ls | grep -v "README"); do
    # Note  - all of these check for the existence of ONE thing with the right format in each column (they don't check that ALL of the rows are the right format)
    # I dont check the following things that are still useful:
    # contig (column 1) - no specific format required
    # function (column 8) - no specific format required
    # The following columns are never used by my programs:
    # column 4 (location) - use columns 1,5,6, and 7 instead.
    # column 9 (aliases) - use the "aliases" file instead.
    # column 10 (figfam)
    # column 11 (evidence codes)
    fmatch=$(cat "${file}" | cut -f 2 | grep -o -P "^fig\|\d+\.\d+\.peg\.\d+$");
    if [ $? -eq 1 ]; then
	echo "ERROR: Gene IDs in raw file ${file} were not in expected format (fig|#.#.peg.# where the first two are the organism ID) or not in the expected place (second column in raw file)";
	STATUS=1;
    fi
    fmatch=$(cat "${file}" | cut -f 3 | grep -o -P "^peg$");
    if [ $? -eq 1 ]; then
	echo "ERROR: No objects of type peg (third column) identified in file ${file}. Only pegs (protein encoding genes) are considered in our clustering analysis!";
	STATUS=1;
    fi
    fmatch=$(cat "${file}" | cut -f 5 | grep -o -P "^\d+$");
    if [ $? -eq 1 ]; then
	echo "ERROR: Gene start location (fifth column) expected to be a number in file ${file}";
	STATUS=1;
    fi
    fmatch=$(cat "${file}" | cut -f 6 | grep -o -P "^\d+$");
    if [ $? -eq 1 ]; then
	echo "ERROR: Stop location (sixth column) expected to be a number in file ${file}";
	STATUS=1;
    fi
    fmatch=$(cat "${file}" | cut -f 7 | grep -o -P "^[+-]$");
    if [ $? -eq 1 ]; then
	echo "ERROR: Strand (seventh column) must be + or - in file ${file}";
	STATUS=1;
    fi
    # Note - NRWYMKSHBVD are ambiguous nucleotides
    # ACGT are the standard nucleotides
    # Anything that isn't one of these is an error.
    # No gaps are allowed.
    fmatch=$(cat "${file}" | cut -f 12 | grep -o -i -P "^[acgtnrwymkshbvd]+$");
    if [ $? -eq 1 ]; then
	echo "ERROR: Nucleotide sequence expected in 12th column in file ${file}";
	STATUS=1;
    fi
    # Note this wont match the header because of the "_" in aa_sequences but its a bit fragile.
    fmatch=$(cat "${file}" | cut -f 13 | grep -o -i -P "^[A-Z]+$")
    if [ $? -eq 1 ]; then
	echo "ERROR: Amino acid sequence expected in 13th column in file ${file}";
	STATUS=1;
    fi
done
cd ..;

# Check genbank files.
cd genbank;

echo "Checking the format of genbank file names..."
for file in $(ls | grep -v "README"); do
    ok=$(echo "${file}" | grep -x -P "\d+\.\d+\.gbk")
    if [ $? -eq 1 ]; then
	echo "ERROR: Genbank file ${file} does not have a name in the expected format ([organismid].gbk)"
	STATUS=1
    fi
done

echo "Checking for existence of genbank files for every organism in the organisms file..."
for org in ${orgmatch}; do
    # File name must exactly be [organismID].gbk
    fmatch=$(ls | grep -w -F "${org}.gbk");
    if [ $? -eq 1 ]; then
	echo "ERROR: No genbank file match for organism ID ${org} - file name must be ${org}.gbk and placed in the genbank folder";
	STATUS=1;
    fi
done

echo "Checking that all organism IDs in the genbank files have an entry in the organisms file..."
for file in $(ls | grep -v "README"); do
    orgid=$(echo "${file}" | grep -o -P "\d+\.\d+");
    echo "${orgid}"
    ok=$(cat ../organisms | cut -f 3 | grep -F -w "${orgid}")
    if [ $? -eq 1 ]; then
	echo "ERROR: Organism ID ${orgid} in genbank file ${file} has no entry in the organism file"
	STATUS=1
    fi
done

cd ..;

# The aliases file is optional but recommended!
if [ ! -f ./aliases/aliases ]; then
    echo "WARNING: No aliases file found - no alias subsitution will be performed for gene names"
fi

exit ${STATUS}