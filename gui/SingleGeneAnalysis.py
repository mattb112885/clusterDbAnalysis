#!/usr/bin/env python

import easygui

import operator
import os
import sqlite3
import sys
import tempfile

from ClusterFuncs import *
from FileLocator import *
from BioPythonGraphics import *

# The user probably doesn't want to see another box if they cancelled it themselves.
class UserCancelError(Exception):
    pass

# This is the base class for GUI errors displaying error messages.
# You can inherit from this class if you want to give a more descriptive name to your errors.
class GuiError(Exception):
    def __init__(self, errormsg):
        msg = "The program encountered the following error:\n\n%s\n\nPress OK to terminate the program.\n" %(errormsg)
        easygui.msgbox(msg=msg)

class NoGeneError(GuiError):
    pass

# codebox makes the font monospaced so we can actually prettyprint (use codebox for tables and alignments).
# But textbox is useful when we want word wrapping.

class ITEPGui:
    # Utilities
    def _geneInfoHeader(self):
        header = [ 'ITEP_geneID', 'Organism', 'Organism ID', '(Ignore)', 'ITEP Contig ID', 'Start', 'Stop', 'Strand', 'Strandsign', 'Function', 'NT sequence', 'AA sequence' ]
        return header
    def _blastHeader(self):
        header = [ 'Query gene', 'Target gene', 'Percent identity', 'HSP length', 'Percent mismatch', 
                   'Gap opens', 'Query start', 'Query end', 'Target start', 'Target end', 'E-value',
                   'Bit score', 'Query self-bit score', 'Target self-bit score' ]
        return header
    def _createTemporaryFile(self, delete=False):
        f = tempfile.NamedTemporaryFile(delete=delete)
        fname = f.name
        return (f, fname)
    def _getClusterId(self):
        # Get the cluster in which the chosen gene is found in the chosen cluster run.
        # Put into its own function because it's so ugly.
        return self.accumulated_data['run_to_cluster'][self.accumulated_data['runid']]
    def _get_run_id(self):
        msg = ''' 
Please choose one of the following methods for defining gene families.

OrthoMCL runs are useful for identifying orthologs (genes likely to share a function)

maxbit runs are useful for identifying broader
gene families. c_xxx in the following list means xxx was used as a cutoff. Higher 
cutoffs mean more stringent similarity to define a family of related genes.

Note that only the groups of organisms that contain your gene are listed here.
'''
        valid_choices = self.accumulated_data['run_to_cluster'].keys()

        if len(valid_choices) == 0:
            easygui.msgbox('The chosen gene is not found in any clustering results!')
            return True

        runid = easygui.choicebox(msg, 'Select a cluster run', valid_choices)

        # Canceling from here - just go back to the other menu
        if runid is None:
            return runid

        self.accumulated_data['runid'] = runid
        return runid
    def _print_readable_table(self, rows, header=True, separator = '|'):
        '''Print a readable table from an array of arrays. '''
        finaltext = ''
        # What is the maximum length of each column?
        numcols = len(rows[0])
        lens = []
        for i in range(numcols):
            col = map(operator.itemgetter(i), rows)
            maxlen = max(map(len, col))
            lens.append(maxlen)
        # Print the header row if there is one.
        if header:
            headers = rows[0]
            header_elements = []
            sep_elements = []
            for i in range(numcols):
                diff = lens[i] - len(headers[i])
                header_elements.append( ' ' + headers[i] + ' '*diff + ' ')
                sep_elements.append( '-'*( lens[i] + 2 ) )
            finaltext += separator.join(header_elements) + "\n" + separator.join(sep_elements) + "\n"
            rows = rows[1:]
        # Print the rest of the rows.
        formats = []
        for line in rows:
            row_elements = []
            for i in range(numcols):
                diff = lens[i] - len(line[i])
                row_elements.append( ' ' + line[i] + ' '*diff + ' ')
            finaltext += separator.join(row_elements) + '\n'
        return finaltext
    # Analyses on a single gene.
    def _get_nucleotide_fasta(self):
        geneinfo = self.accumulated_data['geneinfo']
        text = '>%s %s\n%s\n' %(geneinfo[0], geneinfo[9], geneinfo[10])
        easygui.textbox(text=text)
        return True
    def _get_amino_acid_fasta(self):
        geneinfo = self.accumulated_data['geneinfo']
        text = '>%s %s\n%s\n' %(geneinfo[0], geneinfo[9], geneinfo[11])
        easygui.textbox(text=text)
        return True
    def _get_gene_neighborhood(self):
        self._get_run_id()
        diagram = makeSingleGeneNeighborhoodDiagram(self.accumulated_data['ITEP_id'], self.accumulated_data['runid'], self.sqlite_cursor)
        os.system("display %s" %(diagram))
        return True
    # Analysis Related to getting related genes
    def _get_cluster_blast(self):
        clusterid = self._getClusterId()
        genelist = getGenesInCluster(self.accumulated_data['runid'], clusterid, self.sqlite_cursor)
        blast = getBlastResultsBetweenSpecificGenes(genelist, self.sqlite_cursor)
        blast.insert(0, self._blastHeader())
        text = self._print_readable_table(blast, header=True)
        easygui.codebox(text=text)
        return True
    def _get_cluster_geneinfo(self):
        clusterid = self._getClusterId()
        genelist = getGenesInCluster(self.accumulated_data['runid'], clusterid, self.sqlite_cursor)
        geneinfo = getGeneInfo(genelist, self.sqlite_cursor)
        geneinfo.insert(0, self._geneInfoHeader())
        text = self._print_readable_table(geneinfo, header=True)
        easygui.codebox(text=text)
        return True
    def _get_cluster_fasta(self, amino=True):
        r2c = self.accumulated_data['run_to_cluster']
        clusterid = self._getClusterId()
        genelist = getGenesInCluster(self.accumulated_data['runid'], clusterid, self.sqlite_cursor)
        geneinfo = getGeneInfo(genelist, self.sqlite_cursor)
        if amino:
            idx = 11
        else:
            idx = 10
        text = ''
        for gi in geneinfo:
            text += '>%s %s\n%s\n'%(gi[0], gi[9], gi[idx])
        easygui.textbox(text=text)
        return True
    def _get_presence_absence_table(self):
        (pa_file, pa_fname) = self._createTemporaryFile()
        cluster = self._getClusterId()
        cmd = 'db_getPresenceAbsenceTable.py -r %s -c %s > %s 2> /dev/null' %(self.accumulated_data['runid'], cluster, pa_fname)
        print cmd
        os.system(cmd)
        pa_table = [ line.strip('\r\n').split('\t') for line in pa_file ] 
        text = self._print_readable_table(pa_table)
        easygui.codebox(text=text)
        return True
    def _make_crude_alignment(self):
        (aln_file, aln_fname) = self._createTemporaryFile()
        cluster = self._getClusterId()
        cmd = 'makeTabDelimitedRow.py %s %s | db_makeClusterAlignment.py -m mafft_linsi -n | Gblocks_wrapper.py | db_replaceGeneNameWithAnnotation.py -a -o > %s 2> /dev/null' \
            %(self.accumulated_data['runid'], cluster, aln_fname)
        print cmd
        os.system(cmd)
        text = ''.join( [ line for line in aln_file ] )
        easygui.codebox(text=text)
        return True
    def _make_crude_tree(self):
        # Create tree and make human-readable.
        (nwk_file, nwk_fname) = self._createTemporaryFile()
        cluster = self._getClusterId()
        cmd = 'makeTabDelimitedRow.py %s %s | db_makeClusterAlignment.py -m mafft_linsi -n | Gblocks_wrapper.py | FastTreeMP -wag -gamma | db_replaceGeneNameWithAnnotation.py -a -o > %s 2> /dev/null' \
            %(self.accumulated_data['runid'], cluster, nwk_fname)
        print cmd
        os.system(cmd)
        text = ''.join( [ line for line in nwk_file ] )
        easygui.textbox(text=text)
        return True
    def _display_crude_neighborhood_tree(self):
        # Create tree
        (nwk_file, nwk_fname) = self._createTemporaryFile()
        cluster = self._getClusterId()
        cmd = 'makeTabDelimitedRow.py %s %s | db_makeClusterAlignment.py -m mafft_linsi -n | Gblocks_wrapper.py | FastTreeMP -wag -gamma > %s 2> /dev/null' \
            %(self.accumulated_data['runid'], cluster, nwk_fname)
        print cmd
        os.system(cmd)
        # View tree with neighborhoods
        second_cmd = 'db_makeNeighborhoodTree.py -p %s -r %s -d' %(nwk_fname, self.accumulated_data['runid'])
        print second_cmd
        os.system(second_cmd)
        return True
    # Loop - study genes related to the starting gene.
    def _get_related_genes(self):
        self._get_run_id()
        ok = True
        while ok:
            ok = self._handle_cluster_run_options()
        return True
    def _handle_cluster_run_options(self):
        valid_choices = [ 'Make Amino acid FASTA file', 
                          'Make nucleotide FASTA file', 
                          'Make a crude AA alignment', 
                          'Make a crude Newick tree from AA alignment',
                          'Display a crude tree with neighborhoods attached',
                          'Get a presence and absence table',
                          'Get information on related genes',
                          'Get blast support for a protein family']
        option = easygui.choicebox("What do you want to do with it?", "Choose an analysis", valid_choices)        
        if option is None:
            return False
        if option == 'Make Amino acid FASTA file':
            self._get_cluster_fasta(amino=True)
        elif option == 'Make nucleotide FASTA file':
            self._get_cluster_fasta(amino=False)
        elif option == 'Make a crude AA alignment':
            self._make_crude_alignment()
        elif option == 'Make a crude Newick tree from AA alignment':
            self._make_crude_tree()
        elif option == 'Get a presence and absence table':
            self._get_presence_absence_table()
        elif option == 'Display a crude tree with neighborhoods attached':
            self._display_crude_neighborhood_tree()
        elif option == 'Get information on related genes':
            self._get_cluster_geneinfo()
        elif option == 'Get blast support for a protein family':
            self._get_cluster_blast()
        return True
    # Setup
    def _setUpClusterInfo(self):
        clusterrun_list = getClustersContainingGenes( [ self.accumulated_data['ITEP_id'] ], self.sqlite_cursor)
        run_to_cluster = {}
        for cr in clusterrun_list:
            run_to_cluster[cr[0]] = cr[1]
        self.accumulated_data['run_to_cluster'] = run_to_cluster
    def _setUpGeneInfo(self, alias):
        # Try ITEP ID first
        geneinfo = getGeneInfo( [ alias ], self.sqlite_cursor)
        if len(geneinfo) == 0:
            alias_file = locateAliasesFile()
            alias2gene = {}
            for line in open(locateAliasesFile()):
                spl = line.strip("\r\n").split("\t")
                alias2gene[spl[1]] = spl[0]
            if alias not in alias2gene:
                raise NoGeneError("Sorry, we could not find gene ID %s in the database or in our aliases file. It might not be in this database.\n" %(alias))
            itep_id = alias2gene[alias]
            geneinfo = getGeneInfo( [ itep_id ], self.sqlite_cursor)
        else:
            # ITEP ID was provided
            itep_id = alias

        geneinfo = geneinfo[0]
        self.accumulated_data['alias'] = alias
        self.accumulated_data['ITEP_id'] = itep_id
        self.accumulated_data['geneinfo'] = geneinfo        
        return True
    def __init__(self, cur):
        self.valid_choices = [ 'Nucleotide FASTA', 'Amino acid FASTA', 'Gene neighborhood', 'Related genes in other organisms']
        self.sqlite_cursor = cur
        self.accumulated_data = {}
        return
    # Interface
    def getGeneId(self):
        gene_alias = easygui.enterbox("Please enter the locus tag or ITEP ID of the gene you wish to study.")
        if gene_alias is None:
            raise UserCancelError('User cancelled the operation.')
        self._setUpGeneInfo(gene_alias)
        self._setUpClusterInfo()
        return gene_alias
    def askForChoice(self):
        # Display some information about the gene.
        geneinfo = self.accumulated_data['geneinfo']
        alias = self.accumulated_data['alias']
        msg = '''
You selected %s. Here is some basic information about this gene.

ITEP gene ID: %s
Organism: %s
Organism ID: %s
Contig ID: %s
Start location: %s
Stop location: %s
Strand: %s
Annotated Function: %s

What do you want to know about this gene?
''' %(alias, geneinfo[0], geneinfo[1], geneinfo[2], geneinfo[4], geneinfo[5], geneinfo[6], geneinfo[7], geneinfo[9])
    
        choice = easygui.choicebox(msg, 'Select an analysis.', gui.valid_choices)
    
        if choice is None:
            raise UserCancelError('User clicked CANCEL. No action taken.')
        return choice
    def runChosenAnalysis(self, choice):
        if choice == 'Nucleotide FASTA':
            self._get_nucleotide_fasta()
        elif choice == 'Amino acid FASTA':
            self._get_amino_acid_fasta()
        elif choice == 'Gene neighborhood':
            self._get_gene_neighborhood()
        elif choice == 'Related genes in other organisms':
            self._get_related_genes()

        return True


if __name__ == "__main__":
    print "WARNING! This is highly experimental and will probably break in strange and wonderful ways."

    # Initialization
    con = sqlite3.connect(locateDatabase())
    cur = con.cursor()
    gui = ITEPGui(cur)

    # Lets get a focus gene to study.
    gui.getGeneId()

    # What do you want to do with it?
    while 1:
        choice = gui.askForChoice()
        gui.runChosenAnalysis(choice)

    con.close()
