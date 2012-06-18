#!/usr/bin/python

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

usage = "%prog [groupfile] [Inflation] [cutoff] [scoremethod]"
description = "This file is intended to be run as part of Main1.sh. Run MCL clustering on the organisms specified in the file 'groupfile' with specified parameters (Automatically dumps result file into clusters/ folder)"
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
    # In case the desired organisms have spaces in them (e.g. Methanosarcina acetivorans), we need quotes in order to have a functional command.
    # Use single quotes becuase double quotes cause the ! character to be escaped and mess up the SQL command...
    safepart = ['\'' + k + '\'' for k in spl]
    safewhole = " ".join(s for s in safepart)

    # Generate a filename - no quotes needed or wanted here but we do need to replace spaces
    orgstr = "_".join(s.replace(" ", "_").replace("!","NOT_") for s in spl)
    foutname = "clusters/%s_I_%s_c_%s_m_%s" %(orgstr, inflation, cutoff, method)

    # Check for filename length. If the file name is too long (255 character limit in modern windows, OSX, and unix), cut it off and warn the user
                     
    if len(foutname) > 255:
        sys.stderr.write("WARNING: resulting file name %s is too long and will be truncated to 255 characters. This could cause naming conflicts between run IDs...\n" %(foutname))
        foutname = foutname[0:255]

    # Make sure that file doesn't exist already. If it does, skip over it... mcl takes a long time.
    try:
        fid = open(foutname, "r")
        fid.close()
    except IOError:
        # Generate and run an MCL command
        cmd = ("python src/db_getBlastResultsBetweenSpecificOrganisms.py " + safewhole + 
               " | python src/db_makeBlastScoreTable.py -m " + str(method) + " -c " + str(cutoff) +
               " | mcl - --abc -I " + str(inflation) + " -o " + foutname + ";")
        sys.stderr.write("%s\n" %(cmd) )
        os.system(cmd)
