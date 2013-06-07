#!/usr/bin/python

import optparse
from Bio import SeqIO
from GenbankHandler import *

usage = "addItepIdsToGenbank [genbank_file] [raw_file] [output_file]"
description = """Takes an existing Genbank file and a tab-delimited raw file (from the raw/ 
folder) and adds to it a db_xref for 'ITEP' that contains ITEP IDs"""

parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-t", "--truncate", help="Truncate contig IDs to 16 characters so that biopython can write the genbank file (could lead to non-uniqueness of contig IDs so be careful) (D: Dont - but the function will fail if any contig IDs are larger than 16 characters", action="store_true", dest="truncateContigIds", default=False)
(options, args) = parser.parse_args()

gbk_obj = SeqIO.parse(open(args[0], "r"), "genbank")
tbl_obj = [ line.strip("\r\n").split("\t") for line in open(args[1], "r") ]
modified_gbk = addItepGeneIdsToGenbank(gbk_obj, tbl_obj, truncateContigIds=options.truncateContigIds)

SeqIO.write(modified_gbk, args[2], "genbank")
