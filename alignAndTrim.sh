#!/bin/bash

# Standardized function calls for aligning and trimming alignments
# Accepts input from stdin and outputs to stdout.

# Run alignment (automatically choose an appropriate accuracy level)
# For now the parameters for GBlocks are buried in the wrapper but
# I should bring them out sometime!
mafft --auto -  | \
python src/Gblocks_wrapper.py