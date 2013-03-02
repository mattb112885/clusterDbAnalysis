#!/usr/bin/python

import json
import sys

from sanitizeString import *

from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

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

    annotations = { 'source': organism_name, 'organism': organism_name }

    contig_to_sequence = {}
    for contig in genome_object["contigs"]:
        contig_to_sequence[contig["id"]] = contig["dna"]


    # The contig references are inside the feature definitions in the Genome object file, but
    # in a genbank file the features in a contig must all be separated.
    # Therefore I have to keep track of them in one step and then create the SeqRecord objects
    # in a separate step.

    contig_to_feature_data = {}
    for feature in genome_object["features"]:
        # FIXME - What do I do with things that have more than one location?
        assert(len(feature["location"]) == 1)
        loc = feature["location"][0]
        contig = loc[0]

        # Deal with start and stop locations...
        start = int(loc[1])
        strandstr = loc[2]
        if strandstr == "-":
            strand = -1
        else:
            strand = 1
        featurelen = loc[3]

        # I verified against Pubseed that these semantics and calcualtions are correct, at least
        # for the proteins I checked that are the same between pubseed and KBase...
        if strand == -1:
            stop = start - featurelen + 1
        else:
            stop = start + featurelen - 1            

        # Now I need to convert these into Python slicing indexes...
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

        # Always start < stop in a FeatureLocation...
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

'''    strand = "+"
    neargeneid = "X"
    start = 1
    stop = 10
    seq = Seq("MAAGVKQLADDRTLLMAGVSHDLRTPLTRIRLATEMMSEQDGYLAESINKDIEE", IUPAC.protein)
    qualifiers = {}
    qualifiers["translation"] = "MAAGVKQLADDRTLLMAGVSHDLRTPLTRIRLATEMMSEQDGYLAESINKDIEE"
    qualifiers["product"] = "Phosphofructokinase"
    features = [ SeqFeature(FeatureLocation(start, stop), type="CDS", id="NP_149993.1", qualifiers=qualifiers) ]
    # The SeqRecord holds the CONTIG. We make a separate record for each contig and then spit them all out as a big concatinated
    # genbank file which is one genbank for each contig.
    #
    # Note we have to sanitize the contig names in the ID or else it will think the "." indicates a version number and it will fail.
    annotations = { 'source': 'Escherichia coli K12', 'organism': 'Escherichia coli K12' }
    record = SeqRecord(seq, id="kb_g_0_c_1", description="Escherichia coli K12 contig kb_g_0_c_1", name="kb|g.0.c.1", features=features, annotations=annotations)
    record.source = "Escherichia coli K12"
    records = [record]
    SeqIO.write(record, sys.stdout, "genbank")'''

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
