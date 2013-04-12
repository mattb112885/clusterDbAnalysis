#!/bin/bash

# Source this file to set up your environment.
# To make this happens every time you boot up,
# type 
#
# $ pwd
#
# and then add the following line to your .bashrc file (~/.bashrc):
#
# source [results_from_pwd]/SourceMe.sh

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

echo ""
echo "***********************************"
echo "SourceMe.sh notes from ITEP:"
echo ""
echo "NOTE - SourceMe.sh is intended to be sourced, not run."
echo "If you have run this file as a command (e.g. ./SourceMe.sh)"
echo "you should run source this file instead (we recommend adding"
echo "the source command to your .bashrc file)"
echo ""
echo "If you have sourced SourceMe.sh, you can safely ignore"
echo "this message."
echo ""
echo "The current ITEP database is located at: ${DIR}"
echo "***********************************"
echo ""

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export PYTHONPATH=$PYTHONPATH:$DIR"/lib"
export PATH=$PATH:$DIR"/src"
export PATH=$PATH:$DIR"/src/internal"
export PATH=$PATH:$DIR"/src/utilities"
export PATH=$PATH:$DIR"/scripts"
