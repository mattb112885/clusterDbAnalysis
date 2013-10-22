#!/usr/bin/env python
#
# Use MAFFT to make a quick and dirty alignment between a (small) set of genes of interest.

import fileinput, optparse, os, random, sqlite3, sys
from FileLocator import *

usage="%prog [options] < gene_list > alignment"
description="""Get a quick and dirty (whole-gene) alignment between 
a set of genes of interest piped from stdin. Runs the alignment using
mafft with the --auto and --reorder flags.
Does not do any trimming. To get a nucleotide alignment you can run MAFFT with
the nucleotides directly (-n) or align the proteins and then uses the
PAL2NAL package available at https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1538804/ (-l or --pal2nal)
Aligning the proteins first should give more accuracy and allow calculations of DN\DS and other
such parameters. """

parser = optparse.OptionParser(usage=usage, description=description)
parser.add_option("-g", "--genecol", help="Column number for gene starting from 1 (D=1)", action="store", type="int", dest="gc", default=1)
parser.add_option("-p", "--phylip", help="Output alignment as PHYLIP (D: FASTA)", action="store_true", dest="phylip", default=False)
parser.add_option("-k", "--keep", help="Keep temporary files (D: Delete them)", action="store_true", dest="keep", default=False)
parser.add_option("-n", "--nucleotide", help="Make a nucleotide rather than an amino acid alignment. Runs MAFFT directly on nucleotide sequenced. (D: Amino acids)", 
                  action="store_true", dest="nucleotide", default=False)
parser.add_option("-l", "--pal2nal", help="Make a nucleotide alignment by first making an AA alignment then reverse-translating using PAL2NAL (D: Amino acids)",
                  action = "store_true", dest="pal2nal", default=False)
(options, args) = parser.parse_args()

gc = options.gc - 1

con = sqlite3.connect(locateDatabase())
cur = con.cursor()

id2aaseq = {}
id2nucseq = {}
for line in fileinput.input("-"):
    spl = line.strip("\r\n").split("\t")
    gene = spl[gc]
    cur.execute("SELECT aaseq,nucseq FROM processed WHERE processed.geneid = ?", (gene, ) )
    for res in cur:
        id2aaseq[gene] = str(res[0])
        id2nucseq[gene] = str(res[1])

con.close()

# Make a temporary fasta file for both nucleotides and amino acids
rn = random.randint(0, 2**30)
aa_fname = "%d.faa" %(rn)
aa_fid = open(aa_fname, "w")
nuc_fname = "%d.fna" %(rn)
nuc_fid = open(nuc_fname, "w")
for myid in id2nucseq:
    nuc_fid.write(">%s\n%s\n" %(myid, id2nucseq[myid]))
    aa_fid.write(">%s\n%s\n" %(myid, id2aaseq[myid]))
aa_fid.close()
nuc_fid.close()

# Run MAFFT
if options.nucleotide:
    cmd = "mafft --auto %s > %d.alignment 2> /dev/null" %(nuc_fname, rn)
else:
    cmd = "mafft --auto %s > %d.alignment 2> /dev/null" %(aa_fname, rn)

os.system(cmd)

# If we want to run PAL2NAL.pl we do that now.
if options.pal2nal:
    cmd = "pal2nal.pl -nostderr -output fasta -codontable 11 %d.alignment %s > %d.tmp" %(rn, nuc_fname, rn)
    os.system(cmd)
    cmd = "mv %d.tmp %d.alignment" %(rn, rn)
    os.system(cmd)

# If we want a PHYLIP file I'll give you one but it will have gene IDs replaced...
if options.phylip:
    sys.stderr.write("WARNING: Phylip file requested - gene IDs replaced with 10-character IDs\n")
    os.system("cat %d.alignment | fastaToPhylip.py -c %d.conv" %(rn, rn))
else:
    os.system("cat %d.alignment" %(rn))

# Remove temporary files (including the alignment file)
if not options.keep:
    os.system("rm %d.*" %(rn) )
