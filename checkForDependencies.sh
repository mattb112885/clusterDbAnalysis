#!/bin/bash

# This checks for existence of bash dependencies (Python dependencies will be checked with setuptools)
# Not checking for versions of these though...
STATUS=0;

echo ""
echo "looking for existence of dependencies on your machine..."
echo ""
echo ""
echo "******************************************************************"
echo "* IMPORTANT NOTE - If you are using Ubuntu many packages can     *" 
echo "* also be found via aptitude (Ubuntu's package manager) and it   *"
echo "* is often easier to install them that way if they are missing   *"
echo "******************************************************************"
echo ""
echo ""

# Check for absolute requirements
command -v blastp > /dev/null 2>&1 || {
    echo "";
    echo "ERROR: Unable to find NCBI BLAST+ which is required";
    echo "It can be found at ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/"
    echo "Note that the \"blast\" in aptitude is NOT BLAST+ - you must have e.g. the blastp executable"
    echo "";
    $STATUS=1;
}

command -v mcl > /dev/null 2>&1 || {
    echo "";
    echo "ERROR:Unable to find required dependency MCL (needed for clustering)";
    echo "It can be found at http://micans.org/mcl/"
    echo "";
    $STATUS=1;
}

command -v python > /dev/null 2>&1 || {
    echo "";
    echo "ERROR:Unable to find Python which is required.";
    echo "Please install python 2.6 or 2.7"
    echo "";
    $STATUS=1;
}

command -v sqlite3 > /dev/null 2>&1 || { 
    echo "";
    echo "ERROR:Unable to find required dependency sqlite3";
    echo "We recommend getting the latest version from https://sqlite.org/ (some older versions have bugs that break ITEP)"
    echo "";
    $STATUS=1;
}

# Check for Python packages (note this doens't do version checking)
python -c 'import Bio'
if [ $? -ne 0 ]; then
    echo ""
    echo "ERROR: Unable to find required package Biopython (Bio). This package is needed for many IO and graphics operations in ITEP"
    echo "It can be downloaded from http://biopython.org or installed (via setuptools) with sudo easy_install -f http://biopython.org/DIST/ biopython"
    echo ""
    $STATUS=1;
fi

python -c 'import ete2'
if [ $? -ne 0 ]; then
    echo ""
    echo "ERROR: Unable to find required package ETE (ete2). This package is needed for any tree manipulation in ITEP"
    echo "It can be downloaded at http://ete.cgenomics.org/ or installed (via setuptools) with sudo easy_install -U ete2"
    echo ""
    $STATUS=1;
fi

python -c 'import ruffus'
if [ $? -ne 0 ]; then
    echo ""
    echo "ERROR: Unable to find required package ruffus (ruffus). This package is needed for parallelizing BLAST and RPSBlast calls"
    echo "It can be downlaoded at http://www.ruffus.org.uk/ or installed (via setuptools) with sudo easy_install -U ruffus"
    echo ""
    $STATUS=1;
fi

python -c 'import scipy'
if [ $? -ne 0 ]; then
    echo ""
    echo "ERROR: Unable to find required package scipy. This package is needed for various analyses in ITEP."
    echo "It can be downloaded at http://www.scipy.org/ (follow directions on their website to install)."
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
    echo "It can be downloaded from http://mafft.cbrc.jp/alignment/software/"
    echo "";
}

command -v mysql > /dev/null 2>&1 || {
    echo "";
    echo "WARNING: Unable to find mysql - OrthoMCL requires this for function";
    echo "It can be downloaded from https://dev.mysql.com/downloads/mysql/ "
    echo "";
}

command -v orthomclAdjustFasta > /dev/null 2>&1 || {
    echo "";
    echo "WARNING: Unable to find OrthoMCL - will not be able to use the orthomcl wrapper.";
    echo "It can be downloaded from http://orthomcl.org/common/downloads/ (download 2.0 or greater)"
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
    echo "If you dont have RaxML it can be downloaded at http://www.exelixis-lab.org/ or checked out of github"
    echo "using git clone git@github.com:stamatak/standard-RAxML (requires a github account and SSH key)."
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
    echo "It can be found at http://networkx.github.io/ or installed (using pip) via pip install networkx"
    echo ""
fi

python -c 'import reportlab'
if [ $? -ne 0 ]; then
    echo ""
    echo "WARNING: Unable to find the Python package reportlab. You will need this package to use some Biopython graphics capabilities."
    echo "It can be found at http://www.reportlab.com/ (note you might also have to install fonts - see http://www.reportlab.com/software/installation/ for details)"
    echo ""
fi

python -c 'import matplotlib'
if [ $? -ne 0 ]; then
    echo ""
    echo "WARNING: Unable to find the Python package matplotlib. This package is needed for some of ITEP's plotting functionality."
    echo "It can be found at http://matplotlib.org/"
    echo ""
fi

exit ${STATUS};