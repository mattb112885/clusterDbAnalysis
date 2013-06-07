#!/usr/bin/python

import optparse
from Bio import SeqIO
from GenbankHandler import *

usage = "addItepIdsToGenbank [genbank_file] [raw_file] [output_file]"
description = """Takes an existing Genbank file and a tab-delimited raw file (from the raw/ 
folder) and adds to it a db_xref for 'ITEP' that contains ITEP IDs"""

parser = optparse.OptionParser(usage=usage, description=description)
(options, args) = parser.parse_args()

gbk_obj = SeqIO.parse(open(args[0], "r"), "genbank")
tbl_obj = [ line.strip("\r\n").split("\t") for line in open(args[1], "r") ]
modified_gbk = addItepGeneIdsToGenbank(gbk_obj, tbl_obj)

SeqIO.write(modified_gbk, args[2], "genbank")
