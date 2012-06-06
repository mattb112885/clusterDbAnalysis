#!/bin/sh

# This is a script to run phylip
# seqboot non-interactively

if [ $# -ne 3 ]; then
    echo "USAGE: phylipSeqbootScript.sh [phylip_file] [out_file] [N_replicates]";
    exit 0;
fi

infile=$1;
outfile=$2
reps=$3

# SEQBOOT asks for a filename first (this script assumes you did NOT name the file "infile" like it looks for by default)
str1=${infile}
# Are the settings acceptable? - No, we want to specify the number of replicates.
str2="R"
# This is the number we want
str3=${reps}
# The settings are now acceptable
str4="Y"
# Random number seed
str5="12345"

# Redirect to /dev/null because this will be called as part of a pipe command.
echo ${str1}'\n'${str2}'\n'${str3}'\n'${str4}'\n'${str5} | phylip seqboot > /dev/null

# Move the "outfile" to the desired location
mv outfile ${outfile}