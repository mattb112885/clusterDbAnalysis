#!/bin/bash

# This checks for existence of bash dependencies (Python dependencies will be checked with setuptools)
# Not checking for versions of these though...
STATUS=0;

# Check for absolute requirements
command -v blastp > /dev/null 2>&1 || {
    echo "ERROR: Unable to find NCBI BLAST+ which is required";
    $STATUS=1;
}

command -v mcl > /dev/null 2>&1 || {
    echo "ERROR:Unable to find required dependency MCL (needed for clustering)";
    $STATUS=1;
}

command -v python > /dev/null 2>&1 || {
    echo "ERROR:Unable to find Python which is required.";
    $STATUS=1;
}

command -v sqlite3 > /dev/null 2>&1 || { 
    echo "ERROR:Unable to find required dependency SQLite3";
    $STATUS=1;
}

# Check for useful things
command -v FastTreeMP > /dev/null 2>&1 || {
    echo "WARNING: Unable to find FastTreeMP with the default program name - will need to specify the actual program name to"
    echo "         use the wrapper script for FastTree"
}

command -v Gblocks > /dev/null 2>&1 || {
    echo "WARNING: Unable to find Gblocks with default program name - will need to specify the actual program name to "
    echo "         use the wrapper script for GBlocks"
}

command -v 'raxmlHPC-PTHREADS' > /dev/null 2>&1 || {
    echo "WARNING: Unable to find raxml with the default program name - will need to specify the actual program name to"
    echo "         use the wrapper script for RaxML"
}

command -v mafft > /dev/null 2>&1 || {
    echo "WARNING: Unbale to find mafft - alignment wrappers depending on this will not work"
}

command -v mysql > /dev/null 2>&1 || {
    echo "WARNING: Unable to find mysql - OrthoMCL requires this for function"
}

command -v orthomclAdjustFasta > /dev/null 2>&1 || {
    echo "WARNING: Unable to find OrthoMCL - will not be able to use the orthomcl wrapper."
}

command -v pal2nal.pl > /dev/null 2&>1 || {
    echo "WARNING: Unable to find pal2nal - will be unable to use this option in the fasttree wrapper unless you install it"
}

command -v perl > /dev/null 2>&1 || {
    echo "WARNING: Unable to find Perl - some functionality (RAST integration and OrthoMCL) will not be usable"
}

command -v svr_retrieve_RAST_job > /dev/null 2>&1 || {
    echo "WARNING: Unable to find MyRAST scripts - RAST integration requires svr_retreive_RAST_job from that package"
}

exit ${STATUS};