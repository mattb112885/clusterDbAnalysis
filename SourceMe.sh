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

# Do not echo from things sourced in .bashrc.
# Doing so breaks sftp, scp and who knows what else.
#echo ""
#echo "***********************************"
#echo "SourceMe.sh notes from ITEP:"
#echo ""
#echo "NOTE - SourceMe.sh is intended to be sourced, not run."
#echo "If you have run this file as a command (e.g. ./SourceMe.sh)"
#echo "you should run source this file instead (we recommend adding"
#echo "the source command to your .bashrc file)"
#echo ""
#echo "If you have sourced SourceMe.sh, you can safely ignore"
#echo "this message."
#echo ""
#echo "The current ITEP database is located at: ${DIR}"
#echo "***********************************"
#echo ""

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export PYTHONPATH="$DIR/lib":"$PYTHONPATH"
export PATH="$DIR/src":"$PATH"
export PATH="$DIR/src/internal":"$PATH"
export PATH="$DIR/src/utilities":"$PATH"
export PATH="$DIR/scripts":"$PATH"
