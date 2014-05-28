#!/usr/bin/python

''' This file contains some functions to make it easier to automate fetching genomes from Entrez '''

import easygui
import optparse
import sys
from Bio import Entrez

from GuiBase import *

class NcbiError(GuiError):
    pass

class OrganismImporter(GuiBase):
    ''' 
    You must specify an email so NCBI can nag you if you use too much of their resources.

    This class contains functions for users to identify Genbank files to download.
    The class functions download the selected files and places them in the specified location.
    '''
    def __init__(self, email):
        Entrez.email = email
        return
    def _checkForEntrezError(self, record):
        ''' Raise specific errors if Entrez raised them '''
        if "ErrorList" in record:
            errorstr = ""
            for key in record["ErrorList"]:
                if len(record["ErrorList"][key]) > 0:
                    errorstr += "Error %s: %s\n" %(key, record["ErrorList"][key])
            raise NcbiError("Error querying database: \n\n %s \n " %(errorstr))
        return
    def searchForOrganism(self, orgname):
        ''' Look for organisms matching the specified organism string in the Genome database. '''
        RETMAX=100
        handle = Entrez.esearch(db="genome", term="%s[organism]" %(orgname), retmax=RETMAX)
        record = Entrez.read(handle)
        self._checkForEntrezError(record)
        print record

        idList = record["IdList"]
        # Note - we could run into trouble here if idList is too long.
        # should break it up.
        handle = Entrez.efetch(db="genome", id=idList, rettype = "docsum", retmax=RETMAX)
        records = Entrez.read(handle)
        results_tuples = []
        for record in records:
            self._checkForEntrezError(record)
            genome_id = record["Id"]
            assembly_id = record["AssemblyID"]
            organism_name = record["Organism_Name"]
            results_tuples.append( (genome_id, assembly_id, organism_name) )
        return results_tuples
    def selectOrganism(self, results_tuples):
        ''' Ask user for which of the found organisms he or she wants. '''
        choices = []
        for tup in results_tuples:
            stri = "%s || %s || %s" %(tup[2], tup[0], tup[1])
            choices.append(stri)
        chosen_list = easygui.multchoicebox(msg="Select one or more organisms.\n Format is organism name || Genome ID || Assembly ID. \n", 
                                             title = "Select organisms. ", choices=choices)
        if chosen_list is None or len(chosen_list) == 0:
            raise UserCancelError("User cancelled the operation or did not select any organisms!")

        genome_ids = []
        for choice in chosen_list:
            spl = choice.split("||")
            # We only need the assembly IDs
            genome_id = spl[1].replace(" ", "")
            genome_ids.append(genome_id)
        return genome_ids
    def downloadGenbank(self, genome_ids):
        ''' Download Genbank files for a genome with specified Genome ID '''
        # Ask where are we going to save the files?
        # Get nucleotide IDs (elink to nuccore)
        # Efetch them [Note - will need to do batching here]. Get docsum first.
        # Check the ID ("caption" field) and match to the appropriate database.
        # Make sure "ReplacedBy" is empty.
        # Efetch again but get genbank (gb) this time.
        # Save the files with name = a sanitized version of "Title"

        
if __name__ == "__main__":

    usage = "%prog organism_name email"
    description = """A UI for searching for and downloading Genbank files for a single organism or group of organisms
 from the NCBI Genome database. Note that the search is done by species, and that the Genbank files for individual strains  
will need to be separated and concatenated after this script is done running."""

    parser = optparse.OptionParser(usage=usage, description=description)
    (options, args) = parser.parse_args()

    # FIXME: These should be inputted into the UI.
    if len(args) < 2:
        sys.stderr.write("ERROR: Email address and organism name are required arguments.\n")
        exit(2)

    orgstring = args[0]
    email = args[1]
    importer = OrganismImporter(email)
    found_organisms = importer.searchForOrganism(orgstring)
    selected_organisms = importer.selectOrganism(found_organisms)
