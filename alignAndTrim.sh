#!/bin/bash

# Standardized function calls for aligning and trimming alignments
# Call: ./alignAndTrim.sh [FASTA file]

alignprogram="mafft"

# Run alignment
mafft --auto $1 > /tmp/$1.aln

# Trim alignment
# This uses the "relaxed selection" from the gblocks paper
# which they said works better for shorter alignments.
#
# b1 = "Minimum number of sequences for a conserved position"
# b2 = "Minimum number of sequences for a flank position"
# b3 = "Maximum number of contiguous nonconserved positions"
# b4 = "Minimum length of a block"
# b5 = "Allowed gap positions" (n, h, or a - none, half, or any number)
#
# -p=y: Yes I want you to actually do blocking instead of giving me a menu
# Note that it automatically outputs the same format that it reads (in this case, FASTA)
# the FASTA files it outputs contain extra spaces, which isn't a problem for FASTTREE
# but could be for some other parsers...
Gblocks /tmp/$1.aln -b1=9 -b2=9 -b3=10 -b4=5 -b5=h -p=y

# In case I need it, this is the same thing but with the "stringent" criteria.
#Gblocks /tmp/$1.aln -b1=9 -b2=13 -b3=8 -b4=10 -b5=n -p=y