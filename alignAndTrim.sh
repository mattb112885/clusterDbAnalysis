#!/bin/bash

# Standardized function calls for aligning and trimming alignments

if [ $# -lt 2 ]; then
    echo "Usage: ./src/alignAndTrim.sh [Input file] [Output file]"
    echo "Both input file and output file should be FASTA files"
    echo "Any folder location should be OK for either of them..."
    exit 0;
fi

alignprogram="mafft"

# Get just the file name out of the path
base=$(basename $1)

# Run alignment (automatically choose an appropriate accuracy level)
mafft --auto $1 > /tmp/${base}.aln

# Wrapper for GBlocks that doesnt' whine about long names (which is really just
# long annotations...)
#
# For now the parameters are buried in here - some time I would like to make them
# more explicit.
python src/Gblocks_wrapper.py /tmp/${base}.aln $2