#!/bin/sh

# This is a script to run phylip
# seqboot non-interactively

if [ $# -ne 2 ]; then
    echo "USAGE: phylipSeqbootScript.sh [phylip_file] [out_file]";
    exit 0;
fi

infile=$1;
outfile=$2

# SEQBOOT asks for a filename first
str1=${infile}
# Are the settings acceptable?
str2="Y"
# Random number seed
str3="12345"

echo ${str1}'\n'${str2}'\n'${str3}'\n' | phylip seqboot

# Move the "outfile" to the desired location
mv outfile ${outfile}