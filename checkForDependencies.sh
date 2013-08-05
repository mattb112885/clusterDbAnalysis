#!/bin/bash

# This checks for existence of bash dependencies (Python dependencies will be checked with setuptools)
# Not checking for versions of these though...
STATUS=0;

# Check for absolute requirements
command -v blastp > /dev/null 2>&1 || {
    echo "";
    echo "ERROR: Unable to find NCBI BLAST+ which is required";
    echo "";
    $STATUS=1;
}

command -v mcl > /dev/null 2>&1 || {
    echo "";
    echo "ERROR:Unable to find required dependency MCL (needed for clustering)";
    echo "";
    $STATUS=1;
}

command -v python > /dev/null 2>&1 || {
    echo "";
    echo "ERROR:Unable to find Python which is required.";
    echo "";
    $STATUS=1;
}

command -v sqlite3 > /dev/null 2>&1 || { 
    echo "";
    echo "ERROR:Unable to find required dependency sqlite3";
    echo "";
    $STATUS=1;
}

# Check for Python packages (note this doens't do version checking)
python -c 'import Bio'
if [ $? -ne 0 ]; then
    echo ""
    echo "ERROR: Unable to find required package Biopython (Bio). This package is needed for many IO and graphics operations in ITEP"
    echo ""
    $STATUS=1;
fi

python -c 'import ete2'
if [ $? -ne 0 ]; then
    echo ""
    echo "ERROR: Unable to find required package ETE (ete2). This package is needed for any tree manipulation in ITEP"
    echo ""
    $STATUS=1;
fi

python -c 'import ruffus'
if [ $? -ne 0 ]; then
    echo ""
    echo "ERROR: Unable to find required package ruffus (ruffus). This package is needed for parallelizing BLAST and RPSBlast calls"
    echo ""
    $STATUS=1;
fi

python -c 'import scipy'
if [ $? -ne 0 ]; then
    echo ""
    echo "ERROR: Unable to find required package scipy. This package is needed for various analyses in ITEP."
    echo ""
    $STATUS=1;
fi

# Check for useful things
command -v CompareToBootstrap.pl > /dev/null 2>&1 || {
    echo "";
    echo "WARNING: Unable to find CompareToBootstrap.pl - needed to do gloabl boostraps with FastTree";
    echo "If you want to use it, download it from http://www.microbesonline.org/fasttree/treecmp.html";
    echo "";
}

command -v FastTreeMP > /dev/null 2>&1 || {
    echo "";
    echo "WARNING: Unable to find FastTreeMP with the default program name - will need to specify the actual program name to";
    echo "use the wrapper script for FastTree";
    echo "";
}

command -v Gblocks > /dev/null 2>&1 || {
    echo "";
    echo "WARNING: Unable to find Gblocks with default program name - will need to specify the actual program name to ";
    echo "use the wrapper script for GBlocks";
    echo "If you don't have it, it can be found at http://molevol.cmima.csic.es/castresana/Gblocks.html";
    echo "";
}

command -v mafft > /dev/null 2>&1 || {
    echo "";
    echo "WARNING: Unbale to find mafft - alignment wrappers depending on this will not work";
    echo "";
}

command -v mysql > /dev/null 2>&1 || {
    echo "";
    echo "WARNING: Unable to find mysql - OrthoMCL requires this for function";
    echo "";
}

command -v orthomclAdjustFasta > /dev/null 2>&1 || {
    echo "";
    echo "WARNING: Unable to find OrthoMCL - will not be able to use the orthomcl wrapper.";
    echo "";
}

command -v pal2nal.pl > /dev/null 2&>1 || {
    echo "";
    echo "WARNING: Unable to find pal2nal - will be unable to use this option in the alignment wrappers unless you install it";
    echo "It can be found at http://www.bork.embl.de/pal2nal/";
    echo "";
}

command -v perl > /dev/null 2>&1 || {
    echo "";
    echo "WARNING: Unable to find Perl - some functionality (RAST integration and OrthoMCL) will not be usable";
    echo "";
}

command -v 'raxmlHPC-PTHREADS' > /dev/null 2>&1 || {
    echo "";
    echo "WARNING: Unable to find raxml with the default program name - will need to specify the actual program name to";
    echo "use the wrapper script for RaxML";
    echo "";
}

command -v svr_retrieve_RAST_job > /dev/null 2>&1 || {
    echo "";
    echo "WARNING: Unable to find MyRAST scripts - RAST integration requires svr_retreive_RAST_job from that package";
    echo "If you don't have MyRast you can download it at:";
    echo "http://blog.theseed.org/servers/installation/distribution-of-the-seed-server-packages.html";
    echo "";
}

python -c 'import networkx'
if [ $? -ne 0 ]; then
    echo ""
    echo "WARNING: Unable to find the Python package networkx. You will need to install this to generate GML files for visualization"
    echo "of clusters in other programs like Cytoscape."
    echo ""
fi

python -c 'import reportlab'
if [ $? -ne 0 ]; then
    echo ""
    echo "WARNING: Unable to find the Python package reportlab. You will need this package to use some Biopython graphics capabilities."
    echo ""
fi

python -c 'import matplotlib'
if [ $? -ne 0 ]; then
    echo ""
    echo "WARNING: Unable to find the Python package matplotlib. This package is needed for some of ITEP's plotting functionality."
    echo ""
fi

exit ${STATUS};