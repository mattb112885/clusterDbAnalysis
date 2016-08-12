#!/usr/bin/env python
#
# Make a alignment for a particular cluster and add it to the database
#
# I do all of these steps at once because we have all the inputs required here and otherwise
# the user has to manually re-add all of these things to get everything
# into the db correctly.

from __future__ import print_function
import fileinput, optparse, os, random, sqlite3, sys
from Bio import AlignIO
from FileLocator import *

#### MAFFT methods
# Linsi: Highest accuracy and slowest (recommended for most clusters since they're small)
# Einsi
# Ginsi
# Default: Lowest accuracy and fastest [mafft in > out]
#### CLUSTALW methods
# Default: Just cluster them with default parameters (not as accurate as mafft)
valid_methods = ['mafft_linsi', 'mafft_einsi', 'mafft_ginsi', 'mafft_default', 
	'clustalw_default']

usage="%prog -m Method (-n|-t|-g|-p) [options] < Cluster_RunID > Final_alignment"
description="""Generates a new alignment from a piped-in pair of run\cluster ID """
parser = optparse.OptionParser(usage=usage, description=description)
### Input options
parser.add_option("-r", "--runcol", help="Column number for run ID starting from 1 (D=1)", action="store", type="int", dest="rc", default=1)
parser.add_option("-c", "--clustercol", help="Column number for cluster ID starting from 1 (D=2)", action="store", type="int", dest="cc", default=2)
### Alignment options
parser.add_option("-m", "--alnmethod", help="""Method used to generate alignment (no default)
Valid methods: %s""" %(",".join(valid_methods)), action="store", type="str", dest="method", default=None)
### Trimming options
parser.add_option("-n", "--notrim", help="Do not trim alignment (no default)", action="store_true", dest="notrim", default=False)
parser.add_option("-t", "--trimtomedian", help="Use svr script to trim to median ends (no default)", action="store_true", dest="median", default=False)
parser.add_option("-g", "--gblocks_stringent", help="Use gblocks stringent criteria to trim (no default)", action="store_true", dest="gblocks_stringent", default=False)
parser.add_option("-p", "--gblocks_permissive", help="Use gblocks permissive criteria to trim (no default)", action="store_true", dest="gblocks_permissive", default=False)
### Output options
parser.add_option("-k", "--noclean", help="Do not clean up temporary files (default: Delete them)", action="store_true", dest="noclean", default=False)
parser.add_option("-a", "--addtodb", help="Add new alignment to the database (default: Do not add it)", action="store_true", dest="usedb", default=False)
(options, args) = parser.parse_args()
rc = options.rc - 1
cc = options.cc - 1

if not options.method in valid_methods:
	sys.stderr.write("ERROR: Provided method was not in the available valid methods. Valid methods include:\n")
	sys.stderr.write("\t".join(valid_methods) + "\n")
	exit(2)

# We can only have one trimming method specified
trim_methods = (options.notrim, options.median, options.gblocks_stringent, options.gblocks_permissive,)
sm = sum(trim_methods)
if not sm == 1:
	sys.stderr.write("ERROR: Exactly one of -n, -t, -g, -p (trimming options) must be specified\n")
	exit(2)

#############################
##### Make initial alignment
#############################

# Make a FASTA file containing all the protiens in the specified cluster / run IDs
rann = random.randint(0, 2**30)
fname = "%d.tmp" %(rann)
fid = open(fname, "w")
n = 0
for line in fileinput.input("-"):
	n += 1
	fid.write(line)
fid.close()

if n > 1:
	sys.stderr.write("ERROR: Must provide EXACTLY ONE cluster\run id pair as input\n")
	os.system("rm %s" %(fname))
	exit(2)

fasta = "%d.fasta" %(rann)
cmd = "cat %s | db_getClusterGeneInformation.py -r %d -c %d | annoteSeq2Fasta.py > %s" %(fname, rc+1, cc+1, fasta)

sys.stderr.write("%s\n" %(cmd) )
os.system(cmd)

# Run alignments with specified method
faa = "%d.faa" %(rann)
alncmd = ""
if options.method == 'mafft_linsi':
	alncmd = "mafft --maxiterate 1000 --localpair %s > %s 2> /dev/null" %(fasta, faa)
if options.method == 'mafft_einsi':
	alncmd = "mafft --maxiterate 1000 --genafpair %s > %s 2> /dev/null" %(fasta, faa)
if options.method == 'mafft_ginsi':
	alncmd = "mafft --maxiterate 1000 --globalpair %s > %s 2> /dev/null" %(fasta, faa)
if options.method == 'mafft_default':
	alncmd = "mafft %s > %s 2> /dev/null" %(fasta, faa)
if options.method == 'clustalw_default':
	alncmd = "clustalw -type=PROTEIN -infile=%s -outfile=%s -output=FASTA -align > /dev/null" %(fasta, faa)

os.system("%s" %(alncmd) )

############################
# TRIM ALIGNMENT as specified..
############################

# FIXME - for consistency, I should just change this into the same method as the alignment method
# (i.e. ask for a string and yell if the string is invalid)
trimcmd = ""
trimaln = "%d.faa.trim" %(rann)
trimmethod = "None"
if options.notrim:
	# Just copy it over and don't trim
	trimcmd = "cp %s %s" %(faa, trimaln)
if options.median:
	trimcmd = "cat %s | svr_trim_ali -m > %s " %(faa, trimaln)
	trimmethod = "median"
if options.gblocks_stringent:
	trimcmd = "cat %s | Gblocks_wrapper.py -s > %s " %(faa, trimaln)
	trimmethod = "gblocks_stringent"
if options.gblocks_permissive:
	trimcmd = "cat %s | Gblocks_wrapper.py -r > %s " %(faa, trimaln)
	trimmethod = "gblocks_permissive"

sys.stderr.write("Trimming with command:\n%s\n" %(trimcmd))
os.system(trimcmd)

# Since this is intended to be a pipe command lets get the final alignment to the screen
os.system("cat %s" %(trimaln) )

# Now... lets insert this thing into our lovely database.
if options.usedb:
	con = sqlite3.connect(locateDatabase())
	cur = con.cursor()
	# Note - the ID is autoincremented so we don't need to insert into that column.
	cmd1 = "INSERT INTO alignments VALUES (?,?,?);"
	cmd2 = "INSERT INTO alignedseqs VALUES (?,?,?);"
	cmd3 = "INSERT INTO clusterlists VALUES (?,?,?);"
	cur.execute(cmd1, (rann, options.method, trimmethod) )
	# Read alignment and put into the database.
	recs = AlignIO.read(open(trimaln, "r"), "fasta")
	for rec in recs:
		cur.execute(cmd2, (rann,str(rec.id), str(rec.seq)))
	# Read our temp file back up to get the list of cluster/runIDs
       	for line in open(fname, "r"):
		spl = line.strip("\r\n").split("\t")
		cur.execute(cmd3, (rann, spl[rc], spl[cc]) )
	# Save changes
	con.commit()
	con.close()

# Clean up temporary files
if not options.noclean:
	os.system("rm %d.*" %(rann) )
