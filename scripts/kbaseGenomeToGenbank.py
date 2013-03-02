#!/usr/bin/python

import json
import sys

# CDMI.py is from the KBase - we need it to get the Taxon ID
# Download it at http://kbase.science.energy.gov/developer-zone/downloads/
from CDMI import CDMI_EntityAPI
URL="https://www.kbase.us/services/cdmi_api/"

# This is from the lib/ directory of this repo
from sanitizeString import *

# These are from Biopython
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def getFieldFromRelationship(seedRelationship, fieldName, objtype):
    '''
    INPUTS:
    seedRelationship: The result of one of the get_relationship_xxx functions
    fieldName: The field you want to extract from the object.
    objtype: "TO", "REL", or "FROM"

    OUTPUTS:
    A list (in the same order as the list from the get_relationship function)
    of the values with the specified field name.

    The get_relationship_xxx functions return lists of lists.
    The first list is a list of all the links
    The second list has three dictionaries in it: the TO dictionary, the REL dictionary and the FROM dictionary
    describing properties on either end of the link and of the link itself...

    If you want to  maintain duplicicate relationships (many-to-one, one-to-many, many-to-many), this function should be called at
    least twice (once on each end of the relationship, or once on an end and once in the middle)..
    '''
    if seedRelationship is None:
        sys.stderr.write("INTERNAL ERROR: The provided relationship was None - usually this means you were searching for something that doesn't exist in the database.\n")
        raise ValueError
    objidx = None
    if objtype.lower() == "from":
        objidx = 0
    elif objtype.lower() == "rel":
        objidx = 1
    elif objtype.lower() == "to":
        objidx = 2
    else:
        sys.stderr.write("INTERNAL ERROR: In getFieldFromRelationship - objtype must be TO, REL, or FROM\n")
        raise ValueError
    if not isinstance(seedRelationship, list):
        sys.stderr.write("INTERNAL ERROR: getFieldFromRelationship expects a list - perhaps you meant to call getFieldFromEntity?\n")
        raise ValueError
    # Unravel
    f = []
    for entry in seedRelationship:
        # TO CHECK: Is it possible to have one of the links lead to nothing?
        # Check field name validity - if it links to something there has to be the data request there
        # or else something is wrong.
        if fieldName not in entry[objidx]:
            sys.stderr.write("INTERNAL ERROR: Field name %s not found in provided relationship\n" %(fieldName))
            raise ValueError
        f.append(entry[objidx][fieldName])
    return f

def kbaseGenomeToGenbank(genome_object):
    '''Convert a KBase genome object into a Genbank file incorporating as much info as we can
    as found in the NCBI genbank files.

    Note - the genome object (not to be confused with a ModelSEED "annotation" object) has both annotations / translations
    AND the DNA sequence. It's obtained by calling annotate_genome on an object that only has the DNA sequence.

    Hopefully they won't change this otherwise I'll have to do more cross-referencing and ask for two files. Sigh...'''

    organism_name = genome_object["scientific_name"]
    organism_domain = genome_object["domain"]
    organism_id = genome_object["id"]
    organism_genetic_code = genome_object["genetic_code"]

    # Get the TaxID
    # I  have no idea why the CDMI  isn't working...
    # This is SO FRUSTRATING because the command-line scripts work just fine.
    #cdmi_entity = CDMI_EntityAPI(URL)
    #reldict = cdmi_entity.get_relationship_IsInTaxa(organism_id, [], [], ["id"])
    #taxidlist = getFieldFromRelationship(reldict, "id", "to")
    #taxid = taxidlist[0]

    annotations = { 'source': organism_name, 'organism': organism_name }

    contig_to_sequence = {}
    contig_to_feature_data = {}
    for contig in genome_object["contigs"]:
        contig_to_sequence[contig["id"]] = contig["dna"]
        qualifiers = {}
        qualifiers["organism"] = organism_name
        qualifiers["mol_type"] = "Genomic DNA"
        #qualifiers["db_xref"] = "taxon:%s" %(taxid)
        # We also need to create a "source" feature that will hold the TaxID and some other stuff.
        # We need the TaxID to built our genome IDs in at least a reasonably-consistent manner.
        feature = SeqFeature(FeatureLocation(0, len(contig["dna"])), strand=1, type="source", qualifiers=qualifiers)
        contig_to_feature_data[contig["id"]] = [ feature ]

    # The contig references are inside the feature definitions in the Genome object file, but
    # in a genbank file the features in a contig must all be separated.
    # Therefore I have to keep track of them in one step and then create the SeqRecord objects
    # in a separate step.

    for feature in genome_object["features"]:
        # FIXME - What do I do with things that have more than one location?
        assert(len(feature["location"]) == 1)

        # First lets Deal with start and stop locations...
        # I verified against Pubseed that these semantics and calcualtions are correct, at least
        # for the proteins I checked that are the same between pubseed and KBase...
        loc = feature["location"][0]
        contig = loc[0]
        start = int(loc[1])
        strandstr = loc[2]
        if strandstr == "-":
            strand = -1
        else:
            strand = 1
        featurelen = loc[3]

        if strand == -1:
            stop = start - featurelen + 1
        else:
            stop = start + featurelen - 1
        # Now I need to convert these into Python slicing indexes...because that is what FeatureLocation wants.
        # This includes making the start always less than stop and offsetting the stop by 1 because slide [a,b] only goes up to position b-1
        seqstart = min(start, stop) - 1
        seqstop = max(start, stop)

        feature_id = feature["id"]
        feature_type = feature["type"]

        qualifiers = {}
        # Unfortunately there are features including proteins in the genome objects that have no function (not even "hypothetical protein")
        # Thankfully this isn't a required field in the Genbank file
        if "function" in feature:
            qualifiers["product"] = feature["function"]
        if feature_type == "CDS" or feature_type == "peg":
            qualifiers["translation"] = feature["protein_translation"]
            qualifiers["protein_id"] = feature_id

        # "RNA" is not an official type in a GENBANK file.
        # We attempt to figure out based on the annotation whether it is a tRNA, rRNA, or other (misc_RNA) RNA.
        # These are the offiial RNA types (aside from mRNA but those don't have special fields in the Genome object)
        if feature_type == "rna":
            rRNA_finders = [ "rRNA", "ribosomal", "5S", "16S", "23S", "5.8S", "28S", "18S" ]
            tRNA_finders = [ "tRNA", "transfer" ]
            for finder in rRNA_finders:
                if finder in feature["function"]:
                    feature_type = "rRNA"
            for finder in tRNA_finders:
                if finder in feature["function"]:
                    feature_type = "tRNA"
            if feature_type == "rna":
                feature_type = "misc_RNA"

        # I checked that the above formulas give the correct positions in the genbank file (or at least, the same as the PubSEED genabnk files).
        # FIXME: We still need to check that the program that loads it in and turns them back into start and stop in the
        # raw files gets it right.
        feature = SeqFeature(FeatureLocation(seqstart, seqstop), strand=strand, type=feature_type, id=feature_id, qualifiers=qualifiers)

        if contig in contig_to_feature_data:
            contig_to_feature_data[contig].append(feature)
        else:
            contig_to_feature_data[contig] = [ feature ]

    records = []
    for contig in contig_to_feature_data:
        seq = Seq(contig_to_sequence[contig], IUPAC.ambiguous_dna)
        record = SeqRecord(seq, id=sanitizeString(contig, False), description = "%s contig %s" %(organism_name, contig), name=contig, features=contig_to_feature_data[contig], annotations=annotations)
        records.append(record)
    SeqIO.write(records, sys.stdout, "genbank")

    return

if __name__ == '__main__':

    import json
    import optparse

    usage = "%prog [Genome_object] > Concatinated_genbank_file"
    description = '''Converts a KBase genome JSON object into a genbank file (each contig gets its own file and then they are concatinated)'''
    parser = optparse.OptionParser(usage=usage, description=description)
    (options, args) = parser.parse_args()
    if len(args) < 1:
        sys.stderr.write("ERROR: The JSON genome object is a required argument...\n")
        exit(2)

    genome_object = json.load(open(args[0], "r"))
    kbaseGenomeToGenbank(genome_object)
