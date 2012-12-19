#!/bin/sh

# Script to completely do a core genome analysis based on a specific run id

if [ $# -lt 2 ]; then
    echo ""
    echo "USAGE: completeCoreAnalysis.sh [runid] [rootorg] [NTHREADS (optional)]";
    echo "rootorg should have spaces replaced by underscores"
    echo ""
    echo "DESCRIPTION: Performs common tasks needed to compute a core tree,"
    echo "saving all intermediate files."
    echo ""
    echo "Completed tasks:"
    echo "1. Computation of core gene clusters. Core gene clusters for purposes of"
    echo "   this script are those with exactly one representative in each organism"
    echo "   in the cluster run."
    echo "2. Generation of (protein) FASTA files for each of the core gene clusters."
    echo "3. Alignment with MAFFT (--auto)"
    echo "4. Trimming with GBlocks (with the relaxed settings)"
    echo "5. Concatination of the alignments"
    echo "6. Treeing with FASTTREE (with local bootstraps as in the default)"
    echo ""
    echo "This script can be run from anywhere and will create a lot of files so I suggest"
    echo "running it from its own designated folder."
    echo ""
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
    cat "core_cat_aln_${runid}" | FastTree_wrapper.py -g -b 100 > "core_cat_aln_fasttree_${runid}";
fi

# Reroot tree
if [ ! -f "core_cat_fasttree_rerooted_${runid}".nwk ]; then
    db_displayTree.py -n -b "core_cat_fasttree_rerooted_${runid}" -r "${root}" "core_cat_aln_fasttree_${runid}";
fi

# Get subsets
#makeLeafListFromEachNode.py -f "core_subsets_${runid}" "core_cat_fasttree_rerooted_${runid}"

# TODO:
# Get core gene lists for each of the subsets
# ALL, ALL/ONLY, and NONE

# TODO: 
# Get gene info