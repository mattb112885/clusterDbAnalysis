#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Libraries needed
import fileinput, optparse
import sys, csv, getpass, socket

from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein

'''
Standards set by Matt, 10/30/2012
0: All organism IDs must match the regex \d+\.\d+
1: All raw files must be in the raw/ directory and named [organismid].txt
2: All genbank files must be in the genbank/ directory and named [organismid].gbk
3: All genbank files must have a corresponding raw file.
4: All raw files must have a corresponding genbank file.
5: If there are multiple genbank files (for multiple contigs) they should be concatinated and placed in the genbank folder as a single file with the name [organismid].gbk
6: The raw files must have a specific format (this is to be dealt with automatically if using standardized input functions):
    - Column order as specified in the help files.
    - Gene IDs are in the format fig\|\d+\.\d+\.peg\.\d+ where the first \d+\.\d+ is the organism ID of the corresponding genes
    - Strand must be + or -
    - "Start" is the nucleotide position of the first transcribed base (i.e. for - strand genes, the start position will be bigger than the stop).
'''



#functions
def lookupStrainID(accession):
    #search NCBI for all IDs (the unique numerical key, not just accession)
    genbank = Entrez.read(Entrez.esearch(db="nucleotide", term=accession))
    IDlist = genbank['IdList']
    strainIDs=[]
    #there can me multiple records returned
    for thisID in IDlist:
        summaries = Entrez.read(Entrez.esummary(db="nucleotide", id=thisID))
        for s in summaries:
            strainIDs.append(s['TaxId'])
    strainIDs = list(set(strainIDs))
    assert len(strainIDs)==1, "Note: more than one ID returned from query of NCBI % s" % accession
    return strainIDs[0]


def info_from_genbank(gb_seqrec):
    info = {}
    info["id"]= gb_seqrec.id
    info["gi"]= gb_seqrec.annotations['gi']
    if gb_seqrec.name:
        info["gb_name"] = gb_seqrec.name
    if gb_seqrec.description:
        info["gb_description"] = gb_seqrec.description
    info["taxon"] = lookupStrainID(gb_seqrec.id)
    ##this is another way to get the taxon information, but not as reliable
    #for dbref in gb_seqrec.dbxrefs:
    #    field, value = dbref.split(' ')
    #    if field =="taxon":
    #        info["taxon2"] = value
    info["number_of_features"] = len(gb_seqrec.features)
    numcds = len([f for f in gb_seqrec.features if (f.type =='CDS')])
    info["number_of_cds"] = numcds
    if gb_seqrec.annotations:
        info.update([("gb_annotation: "+k,v) for k, v in gb_seqrec.annotations.items()])
    return info

def info_from_feature(feature):
    info = {}
    info["aa_sequence"] = feature.qualifiers['translation'][0]
    info["aliases"] = feature.qualifiers['protein_id'][0]
    info["function"] = feature.qualifiers['product'][0]
    #need to change strand encoding so that
    info["figfam"] = ""
    info["evidence_codes"] = ""
    if feature.type =='CDS':
        info["type"] = 'peg'
    #add this to preserve other types
    #info["type"] = feature.type
    #Must add one to biopython's 0 indexed to get ohe original genbank one indexed counting
    info["start"] = int(feature.location.start) + 1
    info["stop"] = int(feature.location.end) + 1
    if feature.strand == +1:
        info["strand"] = str("+")
    if feature.strand == -1:
        info["strand"] = str("-")
        #invert the numbers
        info["start"], info["stop"] = info["stop"], info["start"]
    return info

def info_from_record(record):
    info = {}
    info["nucleotide_sequence"] = record.seq.tostring()
    #info["Nucleotide ID"] = record.id
    #info["Nucleotide Description"] = record.description
    #if record.dbxrefs:
    #    info["Database cross-references"] = ";".join(record.dbxrefs)
    return info

def genbank_exstract(name, writefastas=True):
    #set up filenames
    gbin_filename = name + '.gbk'
    nfasta_filename = name + '_genes.fna'
    pfasta_filename = name + '_genes.faa'
    seqfasta_filename = name + '.faa'
    #get data
    gb_seqrec_multi = SeqIO.parse(open(gbin_filename,"r"), "genbank")
    #loop, as there may be multple genbanks together (although this is non-conical)
    #lists to store exstracted seqs
    nfastas = []
    pfastas = []
    genes = []
    for gb_seqrec in gb_seqrec_multi:
        orginfo = info_from_genbank(gb_seqrec)
        for feature in gb_seqrec.features:
            if feature.type =="CDS":
                #check there is only one translation and get info
                assert len(feature.qualifiers['translation'])==1
                geneinfo = {}
                #get aa info
                geneinfo.update(info_from_feature(feature))
                #get na info
                record = feature.extract(gb_seqrec)
                geneinfo.update(info_from_record(record))
                #build output with custom fields (and add them to info list)
                xrefdict = dict([xref.split(':') for xref in feature.qualifiers['db_xref']])
                geneid = "fig|" + str(orginfo["taxon"]) + ".88888" + ".peg." + xrefdict['GI'] #the 1 is arbitrary
                geneinfo["feature_id"] = geneid
                geneinfo["location"] = orginfo["gi"]
                geneinfo["contig_id"] = orginfo["id"]
                geneinfo["source_description"] = orginfo["gb_description"]
                genename = geneinfo["aliases"]
                genedesc = geneinfo["function"] + " " + orginfo["gb_description"]
                geneinfo["gene_description"] = genedesc
                #Save in list to writeout
                precord = SeqRecord(Seq(geneinfo["aa_sequence"], generic_protein),
                                    name=genename,
                                    id = geneid,
                                    description = genedesc
                                    )
                pfastas.append(precord)
                nrecord = SeqRecord(record.seq,
                                    name=genename,
                                    id = geneid,
                                    description = genedesc
                                    )
                nfastas.append(nrecord)
                genes.append(geneinfo)
    if writefastas:
        #writeout  full sequence, and lists
        SeqIO.write(gb_seqrec, open(seqfasta_filename,"w"), "fasta")
        SeqIO.write(nfastas, open(nfasta_filename,"w"), "fasta")
        SeqIO.write(pfastas, open(pfasta_filename,"w"), "fasta")
    return orginfo, genes

def fasta_to_fastas(name, table, i):
    #TODO: This function is for inputing organisms that have ONLY FASTA FILES
    #it is based on an earlier version of the code for Virus files
    #IT NEEDS SERIOUS WORK

    #set up filenames
    seqfastain_filename = name + '_orig.faa'
    seqfasta_filename = name + '.faa'
    nfastain_filename = name + '_genes_orig.fna'
    nfasta_filename = name + '_genes.fna'
    pfasta_filename = name + '_genes.faa'

    #get data and writeout with standard headers
    seq = SeqIO.read(open(seqfastain_filename,"r"), "fasta")
    orginfo = {}
    orginfo["gb_description"] = seq.id
    orginfo["id"] = table['accession'][i]
    seq.id = orginfo["id"]
    seq.name = name
    SeqIO.write(seq, open(seqfasta_filename,"w"), "fasta")
    #lists to store exstracted seqs
    nfastas = []
    pfastas = []
    genes = []
    seqrec = SeqIO.parse(open(nfastain_filename,"r"), "fasta")
    for j, feature in enumerate(seqrec.features):
        geneinfo = {}
        #get aa info
        geneinfo["aa_sequence"] = feature.seq.translate(11).tostring()
        geneinfo["aliases"] = "XX%06i" % j
        geneinfo["function"] = feature.id
        #other ways to do it
        #geneinfo["start"] = int(feature.location.start)
        #geneinfo["stop"] = int(feature.location.end)
        #geneinfo["strand"] = feature.strand
        #get na info
        geneinfo["nucleotide_sequence"] = feature.seq.tostring()
        #build output with custom fields (and add them to info list)
        geneinfo["location"] = orginfo["id"]
        geneinfo["source_description"] = orginfo["gb_description"]
        geneid = "fig|" + orginfo["id"] + ".peg." + geneinfo["aliases"]
        geneinfo["id"] = geneid
        genename = geneinfo["aliases"]
        genedesc = geneinfo["function"] + " from " + orginfo["gb_description"]
        geneinfo["gene_description"] = genedesc
        #Save in list to writeout
        precord = SeqRecord(Seq(geneinfo["aa_sequence"], generic_protein),
                            name=genename,
                            id = geneid,
                            description = genedesc
                            )
        pfastas.append(precord)
        nrecord = SeqRecord(feature.seq,
                            name=genename,
                            id = geneid,
                            description = genedesc
                            )
        nfastas.append(nrecord)
        genes.append(geneinfo)
    #writeout sequences and lists
    SeqIO.write(nfastas, open(nfasta_filename,"w"), "fasta")
    SeqIO.write(pfastas, open(pfasta_filename,"w"), "fasta")
    return orginfo, genes

def findtypes(info):
    return dict([(k, type(v)) for k, v in info.items()])

def fasta2Table(name, table, i):
    #TODO: This function is for inputing organisms that have ONLY FASTA FILES
    #it is based on an earlier version of the code for Virus files
    #IT NEEDS SERIOUS WORK
    orginfo, genes = fasta_to_fastas(name, table, i)
    info = dict(zip(table.dtype.names, table[i]))
    orginfo.update(info)
    orginfos.append(orginfo)
    geneinfos.append(genes)

#Entrez requires an email, if not set namually, this guesses one
email = None
guessemail = getpass.getuser() + '@' + socket.getfqdn()
if email == None:
    email = guessemail
Entrez.email = email

#if run from command line
if __name__ == '__main__':
    #main script function
    usage="%prog < taxonid_list 1> organism_information 2> messages_&_errors"
    description='''
    Given a set of taxonids (from stdin) it will make a table for all genebank files for that taxonid, as well as output a table of all the organism information'''
    parser = optparse.OptionParser(usage=usage, description=description)
    (options, args) = parser.parse_args()

    inputlist = [taxonID.strip() for taxonID in list(fileinput.input("-"))]
    #set up objects to hold information
    orginfos=[]
    info = {}.fromkeys(inputlist)
    #collect information for files
    for i, name in enumerate(inputlist):
        sys.stderr.write("Processing %s" % name)
        geneinfos = []
        orginfo, genes = genbank_exstract(str(name), writefastas=False)
        orginfo.update(info)
        orginfos.append(orginfo)
        geneinfos.append(genes)
        #write out all genes, with collumns in correct order.  the comments are where they originate from
        names = ["contig_id",          #from gi of gb file
                 "feature_id",         #taxonID from NCBI lookup of acention, geneID from the db_xref GI in feature
                 "type",               #CDS = peg
                 "location",           #ignored by iTEP, source assention here
                 "start",              #+1 index and rev if on neg strand
                 "stop",               #+1 index and rev if on neg strand
                 "strand",             #converted to + or -
                 "function",           #from feature product
                 "aliases",            #ignored by iTEP, from gene asention
                 "figfam",             #ignored by iTEP, empty for all
                 "evidence_codes",     #ignored by iTEP, empty for all
                 "nucleotide_sequence",#as recorded in feature
                 "aa_sequence"]        #as recorded in feature, not translated manually
        geneout_filename = str(name)+".1.txt"
        geneout_file = open(geneout_filename, 'w')
        sys.stderr.write(",saved as %s\n" % geneout_filename)
        geneout = csv.DictWriter(geneout_file, fieldnames = names, delimiter="\t")
        #geneout.writeheader()
        for orggenes in geneinfos:
            for gene in orggenes:
                geneout.writerow(dict([(n, gene[n]) for n in names]))
        geneout_file.close()

    #write out all organism, with the original data collumns followed by collumns in alphabetical order
    orgout_file = sys.stdout
    names = orginfos[0].keys()
    names.sort()
    orgout = csv.DictWriter(orgout_file, fieldnames = names)
    orgout.writeheader()
    for org in orginfos:
        orgout.writerow(org)
