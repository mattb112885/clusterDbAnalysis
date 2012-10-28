#!?bin/sh

if [ $# -eq 0 ]; then
    echo ""
    echo "USAGE: addGroupByMatch NAME [match1] [match2] ..."
    echo ""
    echo "DESCRIPTION: Adds a line to the groups file with name NAME"
    echo "containing all organisms matching match1, match2, ..."
    echo ""
    echo "Matching is done against the concents of the organisms file, which must exist"
    echo "If it does not yet exist, call generateOrganismFileFromGenbank first."
    echo ""
    echo "If no matches are provided, a line is added to the groups file containing ALL of the organisms."
    echo ""
    echo "Any organism matching at least one of match1, match2, ... is added."
    echo "The name must match exactly but is not case-sensitive."
    echo ""
    echo "EAXMPLE: addGroupsByMatch e_coli \"Escherichia coli\""
    echo ""
    exit 1
fi

if [ ! -f "organisms" ]; then
    echo "ERROR: Organisms file does not exist!"
    exit 2
fi

GROUPNAME="$1";

# Check for redundancy in the group name.
if [ -f "groups" ]; then
    if cat "groups" | cut -f 1 | grep -q -F -w -i "${GROUPNAME}"; then
	echo "ERROR: Specified group name ${GROUPNAME} already exists in the groups file (note group names are not case-sensitive)."
	exit 3
    fi
fi

if [ $# -eq 1 ]; then
    echo "No match strings provided. Will make a line that matches EVERYTHING."
    MATCHLINES=$(cat "organisms" | cut -f 1)
else
    # Ignore the first argument (group name) and try to match everything else to what is in the organisms file.
    shift;
    MATCHLINES=""
    for var in "$@"; do
	MATCHLINES="$(cat "organisms" | cut -f 1 | grep -i -F "${var}");${MATCHLINES}"
    done
fi

# Remove the trailing semicolon
MATCHLINES=$(echo "${MATCHLINES}" | sed -r "s/(.*)\;$/\1/g")

if [ -z "${MATCHLINES}" ]; then
    echo "ERROR: None of the specified match strings matched organisms in the organism file."
    exit 4
fi

# Check for redundancy in organism list
if [ -f "groups" ]; then
    if cat "groups" | cut -f 2 | grep -q -F -i -w "${MATCHLINES}"; then
	echo "ERROR: The specified match arguments would produce a set of organisms already specified as a group in the groups file."
	exit 5
    fi
fi

# Make the necessary line in the groups file
GROUPLINE=$(makeTabDelimitedRow.py "${GROUPNAME}" "${MATCHLINES}")

