#!/usr/bin/python

'''

This script contains UI functions helpful for building
(relatively quick and dirty) phylogenetic trees for organisms
and analyzing them.


Usage:
ui = PhyloUI()
ui.start()

'''

import easygui
import operator
import os
import sqlite3
import sys

from GuiBase import *
from ClusterFuncs import *
from CoreGeneFunctions import *

class PhyloUI(GuiBase):
    def __init__(self, cur):
        self.sqlite_cursor = cur
        return
    def start(self):
        tasks = [ 'Build an organism tree',
                  'Analyze an organism tree' ]
        msg = "Please choose a task" 
        task = easygui.choicebox(msg, 'Select a task', tasks)
        
        if task == 'Build an organism tree':
            tb = TreeBuilder(self.sqlite_cursor)
            tb.buildOrganismTree()
            pass
        elif task == 'Analyze an organism tree':
            pass

class TreeBuilder(PhyloUI):
    # Init: TreeBuilder(sqlite3_cursor)
    def __init__(self, cur):
        PhyloUI.__init__(self, cur)
    def _get_runid(self):
        ''' This is a slightly different function from the one in singleGeneAnalysis since we aren't limiting it to runs containing a specific gene
        Maybe at some point I'll merge them. '''
        sys.stderr.write("WARNING: _get_runid is not yet implemented. I hardcoded something.\n")
        return "orthomcl_Methanosarcina_typestrains_I_1.5000_c_-5_50"
    def buildOrganismTree(self):
        '''
        Build an organism tree for a specific cluster run.
        Optionally pick a subset of (core) clusters to built it based on.
        '''

        # Set up results locations
        results_location = self._get_directory()
        if results_location is None:
            raise GuiError("You must specify a location to save results files!")
        if not os.path.exists(results_location):
            os.makedirs(results_location)
        fasta_dir = os.path.join(results_location, "fasta_files")
        if not os.path.exists(fasta_dir):
            os.makedirs(fasta_dir)
        alignment_dir = os.path.join(results_location, "alignments")
        if not os.path.exists(alignment_dir):
            os.makedirs(alignment_dir)
        trimmed_dir = os.path.join(results_location, "trimmed_alignments")
        if not os.path.exists(trimmed_dir):
            os.makedirs(trimmed_dir)

        # Now we need a list of core genes.
        runid = self._get_runid()
        sys.stderr.write("Getting a list of core genes in the specified cluster run. Please be patient...\n")
        orglist = getOrganismsInClusterRun(runid, self.sqlite_cursor)
        # ALL and UNIQUE - the genes must have representatives in every organism in the cluster run and have exactly one copy per organism.
        cluster_runs = findGenesByOrganismList(orglist, runid, cl = None, all_org = True, uniq_org = True, pct_cutoff = None, outgroup = None)

        # The user can choose to only make a tree from a subset of the clusters if he or she wants.
        # Subset is the default because it is way faster to compute.
        msg = "Do you want to make a tree from ALL core genes or only a SUBSET of them? (Note that using ALL takes a lot of time)"
        choices = ("Subset", "All")
        choice = easygui.buttonbox(msg=msg, title='', choices=choices)
        if choice is None:
            raise UserCancelError("User cancelled the operation")

        if choice == "Subset":
            # Make a subset
            #     Build a list of the clusters with annotation
            #     User selects which one(s) he or she wants
            ann2id = {}
            sys.stderr.write("Populating list of clusters from which to choose...\n")
            for cr in cluster_runs:
                runid = cr[0]
                clusterid = cr[1]
                ann = findRepresentativeAnnotation(runid, clusterid, self.sqlite_cursor)
                ann += "(ID %s)" %(clusterid)
                ann2id[ann] = (runid, clusterid)
            msg = """Please select one or more of the following clusters to include (click to select and click OK when finished)"""
            choices = ann2id.keys()
            chosen = easygui.multchoicebox(msg=msg, choices=choices)
            final_cluster_set = []
            for choice in chosen:
                final_cluster_set.append(ann2id[choice])
            pass
        elif choice == "All":
            # No: Use all of them.
            final_cluster_set = cluster_runs

        # Make fasta files
        sys.stderr.write("Creating FASTA files for your chosen clusters...\n")
        for cr in final_cluster_set:
            geneinfo = getClusterGeneInfo(cr[0], cr[1], self.sqlite_cursor)
            filename = os.path.join(fasta_dir, cr[0] + "_" + cr[1] + ".faa")
            fid = open(filename, "w")
            for info in geneinfo:
                fid.write(">%s %s\n%s\n" %(info[0], info[9], info[11]))
            fid.close()
        # Align the fasta files.
        sys.stderr.write("Creating protein alignments using MAFFT...\n")
        fasta_files = [ f for f in os.listdir(fasta_dir) if os.path.isfile(os.path.join(fasta_dir,f)) ]
        for f in fasta_files:
            fasta_path = os.path.join(fasta_dir, f)
            # Same name as the input file but different location.
            alignment_path = os.path.join(alignment_dir, f)
            cmd = "mafft --auto \"%s\" > \"%s\"" %(fasta_path, alignment_path)
            print cmd
            os.system(cmd)

        # Trim the alignments.
        alignment_files = [ f for f in os.listdir(alignment_dir) if os.path.isfile(os.path.join(alignment_dir,f)) ]

        # For now I use conservative trimming.
        sys.stderr.write("Trimming alignments (using Gblocks)...\n")
        for f in alignment_files:
            alignment_path = os.path.join(alignment_dir, f)
            trimmed_path = os.path.join(trimmed_dir, f)
            cmd = "cat \"%s\" | Gblocks_wrapper.py -r > \"%s\"" %(alignment_path, trimmed_path)
            print cmd
            os.system(cmd)

        # Concatenate the trimmed alignments.
        sys.stderr.write("Concatenating the trimmed alignments...\n")
        cat_aln_path = os.path.join(results_location, "Concatenated_alignment.faa")
        cmd = "catAlignments.py \"%s\" > \"%s\"" %(alignment_dir, cat_aln_path)
        print cmd
        os.system(cmd)

        # Build a tree with organism IDs.
        sys.stderr.write("Building a tree from the concatenated alignment...\n")
        tree_path = os.path.join(results_location, "Concatenated_protein_tree.nwk")
        cmd = "cat \"%s\" | FastTreeMP -wag -gamma > \"%s\"" %(cat_aln_path, tree_path)

        # Done.
        sys.stderr.write("Done. Use the UI or replaceOrgsWithAbbrev.py to put organism names on the tree leaves.\n")

        return tree_path

if __name__ == "__main__":
    con = sqlite3.connect(locateDatabase())
    cur = con.cursor()
    ui = PhyloUI(cur)
    ui.start()
    con.close()
