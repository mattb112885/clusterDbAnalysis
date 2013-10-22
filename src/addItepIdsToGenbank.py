#!/usr/bin/env python

import optparse
import os
import random
from Bio import SeqIO
from GenbankHandler import *

usage = "addItepIdsToGenbank [genbank_file] [raw_file] [output_file]"
description = """Takes an existing Genbank file and a tab-delimited raw file (from the raw/ 
folder) and adds to it a db_xref for 'ITEP' that contains ITEP IDs"""

parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-t", "--truncate", help="""Replace contig IDs with IDs less than 16 characters so that biopython can write the genbank file.
The original ID will be saved as a db_xref.
(D: Dont - but the function will fail if any contig IDs are larger than 16 characters""", 
                  action="store_true", dest="truncateContigIds", default=False)
(options, args) = parser.parse_args()

gbk_obj = SeqIO.parse(open(args[0], "r"), "genbank")

tbl_obj = [ line.strip("\r\n").split("\t") for line in open(args[1], "r") ]
modified_gbk, newToOriginalName = addItepGeneIdsToGenbank(gbk_obj, tbl_obj)

# Make a temporary file containing replaced IDs
fname = str(random.randint(0,1E10))
tmp_fid = open(fname, "w")
SeqIO.write(modified_gbk, tmp_fid, "genbank")
tmp_fid.close()

tmp_fid = open(fname, "r")
output_fid = open(args[2], "w")
replaceTemporaryIdsWithOriginalIds(tmp_fid, newToOriginalName, output_fid)
output_fid.close()
tmp_fid.close()
os.remove(fname)
