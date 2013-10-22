#!/usr/bin/env python

# Intended to be called from the root directory.
#
# Read the file group file (specified in stdin) 
# and for each group, set up and call a MCL command
# clustering these groups.
#
# First calls the python command to connect to the SQL database and collect
# all of the blast results for the organisms of interest.
# Then, pipe those results to MCL via the mcxdeblast command
#
# Put the results in a file with the names of all the parts of the group
# separated by "_"
#
# Note - this function itself doesn't specifically call the database but depends on two functions that do.

import os, sys
import optparse
from FileLocator import *

usage = "%prog [groupfile] [Inflation] [cutoff] [scoremethod]"
description = """This file is intended to be run as part of Main1.sh. 
Run MCL clustering on the organisms specified in the file 'groupfile' 
with specified parameters (Automatically dumps result file into clusters/ folder)

Alternatively, if scoremethod=orthomcl, runs OrthoMCL with specified inflation and 
(log E-value) cutoff"""

parser = optparse.OptionParser(usage=usage, description=description)
(options, args) = parser.parse_args()

if not len(args) == 4:
    sys.stderr.write("Input arguments missing - try -h for usage details!\n")
    exit(2)

groupfile=args[0]
inflation=args[1]
cutoff=args[2]
method=args[3]

for line in open(groupfile, "r"):
    # Skip comment lines
    if line.startswith("#"):
        continue

    spl = line.strip('\r\n').split('\t')

    # In case the desired organisms have spaces in them (e.g. Methanosarcina acetivorans), 
    # we need quotes in order to have a functional command.
    safepart = ['\'' + k + '\'' for k in spl[1].split(";")]
    safewhole = " ".join(s for s in safepart)

    # Generate a filename - no quotes needed or wanted here but we do need to replace spaces
    # The file name is derived based on the group's name (first column of the groups file)
    foutname = "clusters/%s_I_%s_c_%s_m_%s" %(spl[0], inflation, cutoff, method)
    
    # Make sure that file doesn't exist already. If it does, skip over it... mcl takes a long time.
    # and orthoMCL takes longer than a long time.
    try:
        fid = open(foutname, "r")
        fid.close()
    except IOError:
        # Generate and run an MCL command
        # FIXME - or orthoMCL!!!
        cmd = ("db_getBlastResultsBetweenSpecificOrganisms.py -s " + safewhole + 
               " | makeBlastScoreTable.py -m " + str(method) + " -c " + str(cutoff) +
               " | mcl - --abc -I " + str(inflation) + " -o " + foutname + ";")
        sys.stderr.write("%s\n" %(cmd) )
        os.system(cmd)
