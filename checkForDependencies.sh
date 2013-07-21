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
    echo "ERROR:Unable to find required dependency SQLite3";
    echo "";
    $STATUS=1;
}

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

exit ${STATUS};