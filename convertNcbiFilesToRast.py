#!/usr/bin/python

#####
#
# NITS / Things I need to fix. PLEASE READ THIS CAREFULLY BEFORE RUNNING.
#
# You MUST concatinate all plasmids, chromosomes, contigs, etc. before running this
# or the gene IDs we make in the expected format will not be incremented correctly.
#
#####
#
# THIS FUNCTION WILL NOT WORK WITH EUKARYOTES since it does not yet deal with spliced genes properly. 
#     I would need to find all of the CDS for each labeled gene and stitch them together.
#     FOr example, I see things like this even for the archaea:
#
# NC_003552.1     RefSeq  gene    974836  974945  .       +       .       ID=gene847;Name=MAt4696;Dbxref=GeneID:3362139;gbkey=Gene;locus_tag=MAt4696
# NC_003552.1     RefSeq  tRNA    974836  974873  .       +       .       ID=rna12;Dbxref=GeneID:3362139;gbkey=tRNA;product=tRNA-Met
# NC_003552.1     RefSeq  tRNA    974910  974945  .       +       .       ID=rna12;Dbxref=GeneID:3362139;gbkey=tRNA;product=tRNA-Met
#     [They even have the same ID!] 
# 
# The assoiated location string in the .frn file is:
# :974836-974873,974910-974945
#
# Note that the pubseed doesn't have tab-delimited format available for eukaryotes either, since it isn't very amenable for storing splice sites.
#
#####
#
# THIS FUNCTION DOES NOT YET DEAL WITH TRNA / RRNA. IT RETURNS PROTEINS ONLY
#
# These are in the frn file (not the ffn file), and the location formatting is different.
# Unlike the proteins, there is no "c" for - strand genes, and the DNA sequence
# in the RNA FASTA file is not complemented
# For example, here's a methionine tRNA gene located on the "-" strand:
#
# NC_003552.1     RefSeq  gene    67779   67852   .       -       .       ID=gene55;Name=MAt4684;Dbxref=GeneID:3362135;gbkey=Gene;locus_tag=MAt4684
# NC_003552.1     RefSeq  tRNA    67779   67852   .       -       .       ID=rna0;Dbxref=GeneID:3362135;gbkey=tRNA;product=tRNA-Lys
# NC_003552.1     RefSeq  exon    67779   67852   .       -       .       ID=id1;Parent=rna0;Dbxref=GeneID:3362135;gbkey=tRNA;product=tRNA-Lys
#
# The corresponding entry in the frn file is:
#
# ref|NC_003552|:67779-67852|Lys tRNA| [locus_tag=MAt4684]
#
######
#
# Convert NCBI genome files into
# "RAST" format for use with the rest of my pipeline (or to make your life easier)
#
# The three files you need are:
#
# 1) the "gff" file, which contains this:
#
# contig_id  source   type   start   stop   strand   offset  feature_id_list
#
# 2) The faa file, which contains one of the IDs in the feature_id_list in the header
# 
# 3) The ffn file, which contains the start and stop locations for coding sequences and the resulting sequence.
#
# If there are more than one chromosome, there will be multiple gff, fna, and  faa files. This function
# should be run SEPARATELY for each set of those (so that the start/stop locations of the fna file
# match up correctly), and then concatinated together.
#
# RAST format (for reference) looks like this:
#
# contig_id       feature_id      type    location        start   stop    strand  function        aliases figfam  evidence_codes  nucleotide_sequence     aa_sequence
#
# Many of these will be left blank.
#
#
# Also... the ffn file contains only coding sequences for proteins, it does not contain
# RNAs. Those are found in the frn file and are coded differently.

from Bio import SeqIO
import sys
import re

if not len(sys.argv) == 4:
    print "Usage: convertGenbankToRast.py [gff file] [faa file] [ffn file] > [RAST_tab_delmited_file]"
    exit(2)

aaHeaderToSeq = {}
faa_recs = SeqIO.parse(open(sys.argv[2], "r"), "fasta")
for record in faa_recs:
    aaHeaderToSeq[record.description] = str(record.seq)


nucHeaderToSeq = {}
ffn_recs = SeqIO.parse(open(sys.argv[3], "r"), "fasta")
for record in ffn_recs:
    nucHeaderToSeq[record.description] = str(record.seq)

gffIdFinder = re.compile("protein_id=([^;]*)")
rnaFunctionFinder = re.compile("product=([^;]*)")
taxIdFinder = re.compile("Dbxref=taxon:(\d*)")

# I haven't dealt with RNAs in the pubseed converter either, and
# it is not at all clear to me how the frn file (which contains the 
# nucleotide sequences) works.
# 
# Until I understand that, this function will deal with proteins only.
acceptableTypes = [ "CDS", "region" ]

feature_counter = 1

# We will parse this out of the "region"-tagged line, which I hope is the first one...
# If there are more than one region line they should at LEAST all link to the same
# organism, and in my experience it's always the first line.
org_id = "" 

for line in open(sys.argv[1], "r"):
    # Discard comment lines
    if line.startswith("#"):
        continue

    spl = line.strip('\r\n').split("\t")

    # We want the actual coding regions only.
    if not spl[2] in acceptableTypes:
        continue

    # Only the region line has a link to the TaxID.
    # We need to pull the TaxID out because it gives us a stable id consistent with the PubSEED
    # After pulling it out once we don't need to do it again.
    if spl[2] == "region" and org_id == "":
        org_id = "%s.%d" %(taxIdFinder.search(spl[8]).group(1), 1)
        continue
    elif spl[2] == "region":
        continue

    # Save data we actually need from the gff file
    # and spaces for other stuff we need
    contig = spl[0]
    feature_type = ""
    if spl[2] == "CDS":
        feature_type = "peg"
    location = "%s_%s_%s" %(spl[0], spl[3], spl[4])

    # RAST uses a different convention for "-" strand genes start and stop than
    # NCBI.
    #
    # NCBI forces start to be < stop for all genes including + and - strand genes.
    # RAST makes start < stop for + strand genes and stop < start [backwards] for - strand genes.
    #
    # Lets deal with that here.
    strand = spl[6]
    if strand == "+":
        start = spl[3]
        stop  = spl[4]
    elif strand == "-":
        start = spl[4]
        stop = spl[3]
    else:
        sys.stderr.write("ERROR: Unrecognized strand type %s in the following line: \n %s \n " %(strand, line.strip('\r\n')))
        exit(2)

    # This is just an insanity check. If offset is not 0 I don't know what to do and won't know unless I see an example of it.
    # For now I just look for that case and throw an error...
    _offset = spl[7]
    if not _offset == "0":
        sys.stderr.write("WARNING: Offset is not equal to 0 in the following line: \n %s \n " %(line.strip('\r\n')) )

    idstring = spl[8]

    # Get the ID we will look for within each of the protein AA file keys.
    # Then look for the AA sequence with that ID and save it - will be using it soon!!!
    feature_id = gffIdFinder.search(idstring)
    protseq = ""
    function = ""
    if feature_id == None:
        # for RNAs we need still need an ID number - it'll just be "rnaxx" but whatever. They're arbitrary anyway...
        # But we don't bother looking for an AA sequence since it's an RNA...
        function = rnaFunctionFinder.search(idstring).group(1)
        if org_id == "":
            sys.stderr.write("ERROR: I could not find an organism ID before trying to write protein IDs!\n")
            exit(2)
        feature_id = "fig|%s.rna.%d" %(org_id, feature_counter)
    else:
        # What we actually want to match is the ID, not
        # the "protein_id=" or the ";"
        feature_id = feature_id.group(1)
        # Now lets find which protein has that sequence
        # This is a terribly inefficient method but we'll see how slow it is
        # I imagine it isn't that bad given there's onyl ~4000 proteins...
        done = False
        for header in aaHeaderToSeq:
            if feature_id in header:
                # Check for duplicates
                if done:
                    sys.stderr.write("ERROR: Multiple sequences found for a single protein ID: %s\n" %(protId) )
                    exit(2)
                protseq = aaHeaderToSeq[header]
                function = header
                done = True
        if not done:
            sys.stderr.write("WARNING: No sequence found in the fasta file for protein ID %s ... \n" %(protId) )
            continue
        # Since many scripts I wrote expect a very specific format for the protein IDs
        # I go ahead and assign one.
        if org_id == "":
            sys.stderr.write("ERROR: I could not find an organism ID before attempting to write protein IDs!\n")
            exit(2)
        feature_id = "fig|%s.peg.%d" %(org_id, feature_counter)

    # Nucleotide sequences must be matched up by location on the chromosome.
    # This information is saved in a particular format in the ffn file...

    # The NCBI format is
    # :[start]-[end]
    # for + strand genes, and
    # :c[end]-[start] (as they define it) or :c[start]-[end] (as RAST defines it, and as defined above)
    # for - strand genes.
    optional = ""
    if strand == "-":
        optional = "c"
    stri = ":%s%s-%s" %(optional, start, stop)

    # Search for that location...
    nucseq = ""
    done = False
    for header in nucHeaderToSeq:
        if stri in header:
            if done:
                sys.stderr.write("ERROR: Multiple sequences found matching a single location string... %s\n" %(stri))
            nucseq = nucHeaderToSeq[header]
            done = True
    if not done:
        # Note - this will happen if we try to load up RNAs. This is because the RNAs are in a separate file and the
        # format for encoding their location is not the same.
        sys.stderr.write("WARNING: No sequence found in the fasta file for location: %s ... \n" %(stri) )
        continue

    # Now that we have success, lets increase our feature counter.
    feature_counter += 1

    # Generate a new line to print out
    ln = "\t".join( [ contig, feature_id, feature_type, location, start, stop, strand, function, "", "", "", nucseq, protseq])
    print ln

sys.stderr.write("IMPORTANT: Add the following organism ID to the organisms file: %s \n" %(org_id))

