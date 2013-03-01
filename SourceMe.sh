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
export PYTHONPATH=$PYTHONPATH:$DIR"/lib"
export PATH=$PATH:$DIR"/src"
export PATH=$PATH:$DIR"/scripts"