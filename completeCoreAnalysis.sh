#!/bin/sh

# Script to completely do a core genome analysis based on a specific run id

if [ $# -lt 2 ]; then
    echo "USAGE: completeCoreAnalysis.sh [runid] [rootorg] [NTHREADS (optional)]";
    echo "rootorg should have spaces replaced by underscores"
    exit 0;
fi

runid=$1;
root=$2;

if [ $# -gt 2 ]; then
    NTHREADS="$3";
else
    NTHREADS="1";
fi

# Make directories if needed

mkdir "core_fasta_${runid}" 2> /dev/null
mkdir "core_genelabels_${runid}" 2> /dev/null
mkdir "core_aln_${runid}" 2> /dev/null
mkdir "core_aln_trim_${runid}" 2> /dev/null
mkdir "core_subsets_${runid}" 2> /dev/null

# get geninfo
if [ ! -f "core_geneinfo_${runid}" ]; then
    echo "${runid}" | db_getOrganismsInClusterRun.py | db_findClustersByOrganismList.py -a -u "${runid}" | db_getClusterGeneInformation.py > "core_geneinfo_${runid}";
fi

# Get fasta files
cat "core_geneinfo_${runid}" | getClusterFastas.py "core_fasta_${runid}";

# Align
cd "core_fasta_${runid}";
for file in *; do
    if [ ! -f "../core_aln_${runid}/${file}.aln" ]; then
	mafft --thread ${NTHREADS} --auto "${file}" > "../core_aln_${runid}/${file}.aln";
    fi
done

# Trim
cd "../core_aln_${runid}/";
for file in *; do
    if [ ! -f "../core_aln_trim_${runid}/${file}.trim" ]; then
	cat "${file}" | Gblocks_wrapper.py -r > "../core_aln_trim_${runid}/${file}.trim"
    fi
done
cd ..;

# Cat
if [ ! -f "core_cat_aln_${runid}" ]; then
    catAlignments.py -k "${runid}" "core_aln_trim_${runid}" > "core_cat_aln_${runid}";
fi

# Tree
if [ ! -f "core_cat_aln_fasttree_${runid}" ]; then
    cat "core_cat_aln_${runid}" | FastTreeMP -gamma -wag > "core_cat_aln_fasttree_${runid}";
fi

# Reroot tree
if [ ! -f "core_cat_fasttree_rerooted_${runid}" ]; then
    db_displayTree.py -n -b "core_cat_fasttree_rerooted_${runid}" -r "${root}";
fi

# Get subsets
#makeLeafListFromEachNode.py -f "core_subsets_${runid}" "core_cat_fasttree_rerooted_${runid}"

# TODO:
# Get core gene lists for each of the subsets
# ALL, ALL/ONLY, and NONE

# TODO: 
# Get gene info