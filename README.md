Welcome to the landing page for ITEP, the Integrated Toolkit for the Exploration of microbial Pan-genomes! ITEP is a toolkit for comparative genomics of multiple related
microbial genomes. The capabilities of the tool include:

* de novo protein family prediction by clustering, 
* ortholog detection, 
* analysis of functional domains,
* identification of core and variable genes and gene regions, 
* analysis of gene gain and loss patterns according to the protein family predictions,
* sequence alignments, alignment curation and tree building, 
* curation of annotations and gene calls,
* integration of cross-genome analysis and metabolic networks for study of metabolic network evolution

The tool has been designed to be modular, and as such additional capabilities will be added as needed.

Tutorials for use of this toolkit are available on the [wiki page](https://github.com/mattb112885/clusterDbAnalysis/wiki/) for this project. You can also download them locally 
using:

$ git submodule init
$ git submodule update

(If you have already initialized before you only need to update).

Help text for all Python scripts is available using the "-h" command. In addition we have dumped all of the help text to /doc/help_texts to make for easy searching.

## Setting up paths

Make sure you source the SourceMe.sh file for the repository before trying to run any ITEP commands. This will add appropriate paths. 
To switch between different ITEP databases on the same machine, just source the SourceMe.sh file in the root repository for the database you want to use.

## Database setup scripts

These scripts should be run in the following order to set up the database, once all of the input data is placed correctly (see the 
[wiki](https://github.com/mattb112885/clusterDbAnalysis/wiki) for complete directions). All of these scripts must be run from the
root directory of the repository (as e.g. ./setup_step1.sh).

* setup_step1.sh

    Backs up existing organisms file, generates a new one from the genbank files, and reconciles abbreviations. 

    Checks formatting of input files and throws an error if something is inconsistent.

    Adds organisms to the raw files in $ROOT/raw/ using the data in the "organisms" file.

    (Optionally) adds aliases in $ROOT/aliases/aliases file to the annotations if that file exists.
    
    Set up fasta files and BLAST databases from the raw files in $ROOT/raw (NOT from the genbank files).

    Runs BLASTP and BLASTN all vs. all 

    Computes gene neighborhoods

    Dumps the results into a sqlite database $ROOT/db/Database.SQLite

* setup_step2.sh

    Runs clustering with the specified parameters (cutoff, inflation value and method\metric for clustering) for EVERY group of organisms in the "groups" file
    
    Processes the clustering results into a flat table format for input into the database.

    Pre-computes a presence-absence table for every gene and every organism for each cluster run.

    Imports the results into the sqlite database

    NOTE - If you wish, rather than running MCL on blast results, you can also have the option of importing your own clustering results that are generated using any
    other method you want. The only requirements are that you format the input files correctly and that the gene IDs are ITEP IDs (matching fig\|\d+\.\d+\.peg\.\d+). See the wiki
    and help text for importExternalClustering.py for details.

* setup_step3.sh 

    Parses Genbank files in the genbank/ folder to get whole-genome nucleotide sequences.

    Adds the sequences to the database (WARNING: Do not do this for human-sized genomes!)

* setup_step4.sh

    Downloads a copy of the NCBI CDD if it doesn't already exist

    Runs RPSBLAST for each of your genomes against the CDD

    Formats output for consistency

    Imports results into the database

* setup_step5.sh (EXPERIMENTAL)

    Import user-specified gene calls and cluster information into the database. At the moment this is only used to amend the existing presence-absence tables but other uses are likely forthcoming.

## Other files included in this folder

* addGroupByMatch.py 

    Given a key or a set of keys, searches through the organisms file for all of the organism names matching at least one of the provided keys, 
    and adds all the matching organisms to a new group in the groups file.

    Example: ./addGroupByMatch.py -o organisms -g groups -n Escherichia_and_shigella "Escherichia" "Shigella"
 
    would make a new group called "Escherichia_and_shigella" containing all organisms matching either "Escherichia" or
    "Shigella" in part of their names.

* checkInputFormat.sh 

    Check existence, consistency, and formatting of input files

* dumpDocumentation.sh 

    Automatically generate documentation files for scripts in src/ and in scripts/ and save it to doc/

* generateOrganismFileFromGbk.sh 

    This function is meant for internal use. It will delete the existing organism file and replace 
    with one automatically generated from the genbank files (it looks for the field /organism="[organismname]"
    and pulls out the organism name from that, and gets the organism ID from the file name)

* importAllClusters.sh

    This function performs tasks needed to import all of the clusters found in clusters/ and flatclusters/ into the database
    and generate a presence-absence table.

* removeOrganism.sh 

    Remove all traces of the specified organism, including BLAST results, raw data and genbank file, all clusters 
    containining the specified gene, and all aliases associated with it. Deletes the database file as well.

    By default, this function lists what files and\or lines from data files would be deleted but does NOT delete anything.
    If you specify "TRUE" as an argument, it actually performs the deletion.

## Folders included with this installation

* This folder ($ROOT): Contains scripts needed to set up and maintain the database
* $ROOT/aliases/: Folder in which alias tables should be placed if available
* $ROOT/doc/: Software documentation (including installation directions)
* $ROOT/genbank/: Location for all GENBANK files (see /doc/INSTALL for details)
* $ROOT/raw/: Location for all RAW files (see /doc/INSTALL for details)
* $ROOT/scripts/: Contains "dead-end" \ convenient wrapper scripts to do common analysis tasks.
* $ROOT/src/: Contains modules that can be executed after loading up the data using setup scripts.

## Folders created by the setup scripts

* $ROOT/blastn_res/: BLASTN all vs. all results
* $ROOT/blastres/: BLASTP all vs. all results
* $ROOT/clusters/: Cluster files outputted by MCL
* $ROOT/db/: Contains the SQLite database and tables used as sources for that database
* $ROOT/faa/: Protein FASTA files (generated from raw files for each organism) and their compiled BLAST databases.
* $ROOT/flatclusters: "Flattened" cluster files - runID and clusterID assigned to each cluster and put into a database-friendly format.
* $ROOT/fna/: Nucleotide (gene) FASTA files (generated from raw files for each organism) and their compiled BLAST databases.
* $ROOT/modtable/: Raw tables reformatted for input into database (intended for internal use)
* $ROOT/rpsblast_res/: RPS-BLAST results.