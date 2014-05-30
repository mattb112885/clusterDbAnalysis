#!/usr/bin/env python

import easygui

import operator
import os
import shutil
import sqlite3
import subprocess
import sys
import tempfile

from ClusterFuncs import *
from ClusterGraph import *
from FileLocator import *
from BioPythonGraphics import *
from GuiBase import *
from sanitizeString import *

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

# NOTE: codebox makes the font monospaced so we can actually prettyprint (use codebox for tables and alignments).
# textbox is useful when we want word wrapping.
class ITEPGui(GuiBase):
    # Utilities
    def _geneInfoHeader(self):
        header = [ 'ITEP_geneID', 'Organism', 'Organism ID', '(Ignore)', 'ITEP Contig ID', 'Start', 'Stop', 'Strand', 'Strandsign', 'Function', 'NT sequence', 'AA sequence' ]
        return header
    def _blastHeader(self):
        header = [ 'Query gene', 'Target gene', 'Percent identity', 'HSP length', 'Percent mismatch', 
                   'Gap opens', 'Query start', 'Query end', 'Target start', 'Target end', 'E-value',
                   'Bit score', 'Query self-bit score', 'Target self-bit score' ]
        return header
    def _tblastnHeader(self):
        header = ["Query gene", "Length of query gene", "Target contig", "Target organism", 
                  "Start of target sequence", "End of target sequence", "HSP length", 
                  "Percent of query overlapping", "E-value", "Bit score", "Frame of hit",
                  "StrandedString", "Called gene in hit region", "Annotation of called gene", "Length of called gene",
                  "Percent of called gene overlapping with HSP", "tBLASTn hit ID" ]      
        return header
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
    # Analyses on a single gene.
    def _get_nucleotide_fasta(self):
        geneinfo = self.accumulated_data['geneinfo']
        text = '>%s %s\n%s\n' %(geneinfo[0], geneinfo[9], geneinfo[10])
        easygui.textbox(text=text)
        output_file = self._save_file_dialogs(extension="fasta")
        if output_file is not None:
            self._save_text(text, output_file)
            self._success_dialog(output_file)
        return True
    def _get_amino_acid_fasta(self):
        geneinfo = self.accumulated_data['geneinfo']
        text = '>%s %s\n%s\n' %(geneinfo[0], geneinfo[9], geneinfo[11])
        easygui.textbox(text=text)
        output_file = self._save_file_dialogs(extension="fasta")
        if output_file is not None:
            self._save_text(text, output_file)
            self._success_dialog(output_file)
        return True
    def _get_gene_neighborhood(self):
        self._get_run_id()
        diagram = makeSingleGeneNeighborhoodDiagram(self.accumulated_data['ITEP_id'], self.accumulated_data['runid'], self.sqlite_cursor, labeltype = 'aliases')
        os.system("display %s" %(diagram))
        output_file = self._save_file_dialogs(extension="png")
        if output_file is not None:
            shutil.copyfile(diagram, output_file)
            self._success_dialog(output_file)
        return True
    def _get_similar_genes(self, blastn=False):
        blastres = getBlastResultsContainingGenes( [ self.accumulated_data['ITEP_id'] ], self.sqlite_cursor, blastn=blastn )
        blastres.insert(0, self._blastHeader())
        text = self._print_readable_table(blastres, header=True)
        easygui.codebox(text=text)
        output_file = self._save_file_dialogs(extension="txt")
        if output_file is not None:
            self._save_text(text, output_file)
            self._success_dialog(output_file)
        return True
    def _run_tblastn(self):
        self._get_run_id()
        # Organism file has to be closed to be used by the tblastn.
        # So we have to clean it up manually. The results file though we can
        # have cleaned up automatically.
        (orgf, orgfname) = self._createTemporaryFile(delete=False)
        (resf, resfname) = self._createTemporaryFile(delete=False)
        orglist = getOrganismsInClusterRun(self.accumulated_data['runid'], self.sqlite_cursor)
        orgids = []
        for org in orglist:
            orgid = organismNameToId(org, self.sqlite_cursor)
            orgf.write(orgid + "\n")
        orgf.close()
        sys.stderr.write("Now running the tblastn wrapper (this could take some time...)\n")
#        f1 = subprocess.Popen([ "echo", self.accumulated_data['ITEP_id'] ], stdout=subprocess.PIPE)
#        f2 = subprocess.Popen(["db_TBlastN_wrapper.py", "-f", orgfname, "-r", "1" ], stdin=f1.stdout, stdout=resf)
#        f1.stdout.close()
#        f2.communicate()
        cmd = "echo '%s' | db_TBlastN_wrapper.py -f %s -r 1 > %s" %(self.accumulated_data['ITEP_id'], orgfname, resfname)
        print cmd
        os.system(cmd)
        tblastn_results = [ line.strip("\r\n").split("\t") for line in resf ]
        tblastn_results.insert(0, self._tblastnHeader())
        text = self._print_readable_table(tblastn_results, header=True)
        easygui.codebox(text=text)
        output_file = self._save_file_dialogs(extension="txt")
        if output_file is not None:
            self._save_text(text, output_file)
            self._success_dialog(output_file)      
        # Clean up temp files
        os.remove(orgfname)
        os.remove(resfname)
        return True
    def _get_conserved_domains(self):
        tmpdir = tempfile.mkdtemp()       
        f1 = subprocess.Popen(["echo", self.accumulated_data['ITEP_id'] ], stdout=subprocess.PIPE )
        # Limiting hits to 8 is probably more reasonable than applying a generic Evalue cutoff.
        f2 = subprocess.Popen(["db_displayExternalClusterHits.py", "-o", tmpdir, "--maxhits", "15", "--showevalue" ], stdin=f1.stdout)
        f1.stdout.close()  # Allow echo to receive a SIGPIPE if the second command exits before the first.
        f2.communicate()
        outfiles = [ f for f in os.listdir(tmpdir) if os.path.isfile(os.path.join(tmpdir,f)) ]
        if len(outfiles) == 0:
            easygui.msgbox(msg = "Sorry, the specified gene had no external cluster hits")
        else:
            for outfile in outfiles:
                f3 = subprocess.Popen(["display", os.path.join(tmpdir, outfile)])
                f3.communicate()

        output_file = self._save_file_dialogs(extension="svg")
        if output_file is not None:
            f4 = subprocess.Popen(["cp", os.path.join(tmpdir, outfile), output_file])
            f4.communicate()
        shutil.rmtree(tmpdir)
        return True
    # Analysis Related to getting related genes
    def _get_cluster_blast(self):
        clusterid = self._getClusterId()
        genelist = getGenesInCluster(self.accumulated_data['runid'], clusterid, self.sqlite_cursor)
        blast = getBlastResultsBetweenSpecificGenes(genelist, self.sqlite_cursor)
        blast.insert(0, self._blastHeader())
        text = self._print_readable_table(blast, header=True)
        easygui.codebox(text=text)
        output_file = self._save_file_dialogs(extension="txt")
        if output_file is not None:
            self._save_text(text, output_file)
            self._success_dialog(output_file)
        return True
    def _get_cluster_geneinfo(self):
        clusterid = self._getClusterId()
        genelist = getGenesInCluster(self.accumulated_data['runid'], clusterid, self.sqlite_cursor)
        geneinfo = getGeneInfo(genelist, self.sqlite_cursor)
        geneinfo.insert(0, self._geneInfoHeader())
        text = self._print_readable_table(geneinfo, header=True)
        easygui.codebox(text=text)
        output_file = self._save_file_dialogs(extension="txt")
        if output_file is not None:
            self._save_text(text, output_file)
            self._success_dialog(output_file)
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
        output_file = self._save_file_dialogs(extension="fasta")
        if output_file is not None:
            self._save_text(text, output_file)
            self._success_dialog(output_file)
        return True
    def _get_presence_absence_table(self):
        (pa_file, pa_fname) = self._createTemporaryFile(delete=True)
        cluster = self._getClusterId()
        cmd = 'db_getPresenceAbsenceTable.py -r %s -c %s > %s 2> /dev/null' %(self.accumulated_data['runid'], cluster, pa_fname)
        print cmd
        os.system(cmd)
        pa_table = [ line.strip('\r\n').split('\t') for line in pa_file ] 
        text = self._print_readable_table(pa_table)
        easygui.codebox(text=text)
        output_file = self._save_file_dialogs(extension="txt")
        if output_file is not None:
            self._save_text(text, output_file)
            self._success_dialog(output_file)
        return True
    def _make_crude_alignment(self):
        (aln_file, aln_fname) = self._createTemporaryFile(delete=True)
        cluster = self._getClusterId()
        cmd = 'makeTabDelimitedRow.py %s %s | db_makeClusterAlignment.py -m mafft_linsi -n | Gblocks_wrapper.py | db_replaceGeneNameWithAnnotation.py -a -o > %s 2> /dev/null' \
            %(self.accumulated_data['runid'], cluster, aln_fname)
        print cmd
        os.system(cmd)
        text = ''.join( [ line for line in aln_file ] )
        easygui.codebox(text=text)
        output_file = self._save_file_dialogs(extension="fasta")
        if output_file is not None:
            self._save_text(text, output_file)
            self._success_dialog(output_file)
        return True
    def _make_crude_tree(self, replacename = True):
        # Create tree and make human-readable.
        (nwk_file, nwk_fname) = self._createTemporaryFile(delete=True)
        cluster = self._getClusterId()
        if replacename:
            cmd = 'makeTabDelimitedRow.py %s %s | db_makeClusterAlignment.py -m mafft_linsi -n | Gblocks_wrapper.py | FastTreeMP -wag -gamma | db_replaceGeneNameWithAnnotation.py -a -o > %s 2> /dev/null' \
                %(self.accumulated_data['runid'], cluster, nwk_fname)
        else:
            cmd = 'makeTabDelimitedRow.py %s %s | db_makeClusterAlignment.py -m mafft_linsi -n | Gblocks_wrapper.py | FastTreeMP -wag -gamma > %s 2> /dev/null' \
                %(self.accumulated_data['runid'], cluster, nwk_fname)           
        print cmd
        os.system(cmd)
        text = ''.join( [ line for line in nwk_file ] )
        easygui.textbox(text=text)
        output_file = self._save_file_dialogs(extension="nwk")
        if output_file is not None:
            self._save_text(text, output_file)
            self._success_dialog(output_file)
        return True
    def _display_crude_neighborhood_tree(self):
        # Unlike other commands we need to know if we are saving the results BEFORE we run it.
        output_file = self._save_file_dialogs(extension="png")
        if output_file is not None:
            output_file = output_file[0:len(output_file)-4]


        # Create tree
        (nwk_file, nwk_fname) = self._createTemporaryFile(delete=True)
        cluster = self._getClusterId()
        cmd = 'makeTabDelimitedRow.py %s %s | db_makeClusterAlignment.py -m mafft_linsi -n | Gblocks_wrapper.py | FastTreeMP -wag -gamma > %s 2> /dev/null' \
            %(self.accumulated_data['runid'], cluster, nwk_fname)
        print cmd
        os.system(cmd)

        # View tree with neighborhoods
        second_cmd = 'db_makeNeighborhoodTree.py -p %s -r %s -d -l' %(nwk_fname, self.accumulated_data['runid'])

        if output_file is not None:
            second_cmd += " -o %s --png" %(output_file)

        print second_cmd
        os.system(second_cmd)
        return True
    def _make_cluster_gml_file(self):
        output_file = self._save_file_dialogs(extension="gml")
        if output_file is None:
            return True
        clusterid = self._getClusterId()
        genelist = getGenesInCluster(self.accumulated_data['runid'], clusterid, self.sqlite_cursor)
        blastres = getBlastResultsBetweenSpecificGenes(genelist, self.sqlite_cursor)
        graph = makeNetworkObjectFromBlastResults(blastres, "maxbit", 0.1, self.sqlite_cursor )
        exportGraphToGML(graph, output_file)
        easygui.msgbox(msg="GML file saved to %s. Import into Cytoscape or similar programs to view. A VizMapper file is available at lib/ITEP_vizmapper.props." %(output_file))
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
                          'Make a crude Newick tree with ITEP IDs',
                          'Display a crude tree with neighborhoods attached',
                          'Get a presence and absence table',
                          'Get information on related genes',
                          'Make a GML file to import into Cytoscape',
                          'Get blast support for a protein family']
        option = easygui.choicebox("What do you want to do with the related genes?", "Choose an analysis", valid_choices)        
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
        elif option == 'Make a crude Newick tree with ITEP IDs':
            self._make_crude_tree(replacename=False)
        elif option == 'Get a presence and absence table':
            self._get_presence_absence_table()
        elif option == 'Display a crude tree with neighborhoods attached':
            self._display_crude_neighborhood_tree()
        elif option == 'Get information on related genes':
            self._get_cluster_geneinfo()
        elif option == 'Get blast support for a protein family':
            self._get_cluster_blast()
        elif option == 'Make a GML file to import into Cytoscape':
            self._make_cluster_gml_file()
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
        # Support either sanitized or unsanitized versions.
        itep_id = unsanitizeGeneId(alias)
        geneinfo = getGeneInfo( [ itep_id ], self.sqlite_cursor)
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

        geneinfo = geneinfo[0]
        self.accumulated_data['alias'] = alias
        self.accumulated_data['ITEP_id'] = itep_id
        self.accumulated_data['geneinfo'] = geneinfo        
        return True
    def __init__(self, cur):
        self.valid_choices = [ 'Nucleotide FASTA', 
                               'Amino acid FASTA', 
                               'Gene neighborhood',
                               'Get similar genes by BLASTP',
                               'Get similar genes by BLASTN',
                               'Run tBLASTn against a group of organisms',
                               'Related genes in other organisms',
                               'Show conserved domain hits',
                               ]
        self.sqlite_cursor = cur
        self.accumulated_data = {}
        return
    # Interface
    def getGeneId(self):
        gene_alias = easygui.enterbox("Please enter the locus tag or ITEP ID (sanitized or not) of the gene you wish to study.")
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
        elif choice == 'Run tBLASTn against a group of organisms':
            self._run_tblastn()
        elif choice == 'Get similar genes by BLASTP':
            self._get_similar_genes(blastn=False)
        elif choice == 'Get similar genes by BLASTN':
            self._get_similar_genes(blastn=True)
        elif choice == 'Show conserved domain hits':
            self._get_conserved_domains()

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
