#!/usr/bin/env python

import optparse
import os
import random
from Bio import SeqIO
from GenbankHandler import *

usage = "%prog genbank_file raw_file output_file"
description = """Takes an existing Genbank file and a tab-delimited raw file (from the raw/ 
folder) and adds to it a db_xref for 'ITEP' that contains ITEP IDs"""

parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-t", "--truncate", help="""Replace contig IDs with IDs less than 16 characters so that biopython can write the genbank file.
The original ID will be saved as a db_xref.
(D: Dont - but the function will fail if any contig IDs are larger than 16 characters""", 
                  action="store_true", dest="truncateContigIds", default=False)
parser.add_option("-b", "--tblfile", help="Specify this option if instead of a RAW file you are using the 'tbl' file from RAST tarball ($basedir/Features/peg/tbl)",
                  action="store_true", dest="tblfile", default=False)
(options, args) = parser.parse_args()

gbk_obj = SeqIO.parse(open(args[0], "r"), "genbank")

tbl_obj = []
ignoreseq = False

if options.tblfile:
    # The 'tbl' file is in the following format:
    # RASTID      Contig_start_stop
    # We need to turn this into raw format to as good an approximation as possible
    # given our limited knowledge...
    import re
    regex = re.compile(r"^(.*)_(\d+)_(\d+)$")
    rastmatch = re.compile(r"fig\|\d+\.\d+\.peg\.\d+")
    tbl_obj = []
    for line in open(args[1], "r"):
        spl = line.strip("\r\n").split("\t")
        rast_id = spl[0]
        matchobj = rastmatch.match(rast_id)
        if matchobj is None:
            raise IOError("Incorrectly formatted input file - Rast IDs should be in first column (fig\|\d+\.\d+\.peg\.\d+).")

        matchobj = regex.search(spl[1])
        if matchobj is None:
            raise IOError("Incorrectly formatted input file - second column should be in format CONTIG_STARTLOC_STOPLOC")

        contig = matchobj.group(1)
        start = matchobj.group(2)
        stop = matchobj.group(3)
        tbl_obj.append( [ contig, rast_id, "peg", "", start, stop, "", "", "", "", "", "", "" ] )
    # There is no sequence info in this file so we don't use it
    ignoreseq = True
else:
    tbl_obj = [ line.strip("\r\n").split("\t") for line in open(args[1], "r") ]
    ignoreseq = False

modified_gbk, newToOriginalName = addItepGeneIdsToGenbank(gbk_obj, tbl_obj, ignoreseq=ignoreseq)

# Make a temporary file containing replaced IDs so biopython will let us write the file.
fname = str(random.randint(0,1E10))
tmp_fid = open(fname, "w")
SeqIO.write(modified_gbk, tmp_fid, "genbank")
tmp_fid.close()

# Hacky way to get the original IDs back
tmp_fid = open(fname, "r")
output_fid = open(args[2], "w")
replaceTemporaryIdsWithOriginalIds(tmp_fid, newToOriginalName, output_fid)
output_fid.close()
tmp_fid.close()
os.remove(fname)
