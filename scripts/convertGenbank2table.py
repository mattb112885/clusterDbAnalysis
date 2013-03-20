#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Libraries needed
import fileinput, optparse
import os, sys, csv, getpass, socket, shutil, re

from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein

from FileLocator import *

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
    # Not all Genbank files (i.e. those from RAST fall into this category) actually
    # assign the contig ID in the ID field but they call it a "name" instead...
    # Biopython automatically gives those an ID of "unknown".
    if gb_seqrec.id == "unknown":
        info["id"] = gb_seqrec.name
    else:
        info["id"] = gb_seqrec.id
    
    info["number_of_features"] = len(gb_seqrec.features)
    numcds = len([f for f in gb_seqrec.features if (f.type =='CDS')])
    info["number_of_cds"] = numcds

    if gb_seqrec.description:
        info["gb_description"] = gb_seqrec.description
    else:
        info["gb_description"] = "unknown"

    # The "try" part will only work for Entrez-derived GBK files but we want to be more
    # generic. If we can't find it from tehre we try to get it from the dbxref and if it fails there too
    # there's nothing we can do.
    try:
        info["gi"]= gb_seqrec.annotations['gi']
        if gb_seqrec.name:
            info["gb_name"] = gb_seqrec.name
        info["taxon"] = lookupStrainID(gb_seqrec.id)
    except:
        for record in gb_seqrec.features:
            if record.type == "source":
                for xref in record.qualifiers["db_xref"]:
                    if "taxon" in xref:
                        taxid = re.search("(\d+)", xref).group(1)
                        info["taxon"] = taxid
                        break
            if "taxon" in info:
                break

    if "taxon" not in info:
        sys.stderr.write("ERROR: Could not identify taxonID from querying Entrez or from the Genbank file itself...\n")
        exit(2)

    if gb_seqrec.annotations:
        info.update([("gb_annotation: "+k,v) for k, v in gb_seqrec.annotations.items()])

    return info

def info_from_feature(feature):
    info = {}
    # An extra space in the amino acid sequences exists which 
    # throws off the validator and may cause downstream problems.
    info["aa_sequence"] = feature.qualifiers['translation'][0].strip()
    if "protein_id" in info:
        info["aliases"] = feature.qualifiers['protein_id'][0]
    else:
        info["aliases"] = ""
    if "product" in feature.qualifiers:
        info["function"] = feature.qualifiers['product'][0]
    else:
        info["function"] = "NO_FUNCTION_ASSIGNED"
    info["figfam"] = ""
    info["evidence_codes"] = ""
    if feature.type =='CDS':
        info["type"] = 'peg'
    #add this to preserve other types
    #info["type"] = feature.type
    #Must add one to biopython's 0 indexed to get the original genbank one indexed counting
    #However the stop is slice notation so it is already 1 higher than it "should" be.
    # Thus we only add 1 to the start.
    # I verified that this gives the correct answer relative to what is in the PubSEED
    # for some + and - strand genes in E. coli.
    info["start"] = int(feature.location.start) + 1
    info["stop"] = int(feature.location.end)
    if feature.strand == +1:
        info["strand"] = str("+")
    if feature.strand == -1:
        info["strand"] = str("-")
        # invert the start and stop locations...
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

def genbank_extract(ptr, version_number):
    '''Extract data from a genbank file...returns some organism-specific data,
    a dictionary from gene ID to a list of aliases to that gene ID (including the
    original IDs from the genbank file), and the data needed to build the raw data
    table required by ITEP.'''
    #get data
    gb_seqrec_multi = SeqIO.parse(ptr, "genbank")
    #loop, as there may be multple genbanks together (although this is non-canonical)
    #lists to store extracted seqs
    genes = []
    geneidToAlias = {}
    counter = 0
    for gb_seqrec in gb_seqrec_multi:
        orginfo = info_from_genbank(gb_seqrec)
        for feature in gb_seqrec.features:
            if feature.type =="CDS":
                #check there is only one translation and get info
                #TODO - Should we attempt to translate manually if no translations are preent in the gbk file?
                if 'translation' not in feature.qualifiers:
                    sys.stderr.write("WARNING: CDS found with no translation\n")
                    sys.stderr.write("Qualifiers:\n")
                    for key in feature.qualifiers:
                        sys.stderr.write("%s\t%s\n" %(key, feature.qualifiers[key]))
                    continue
                assert len(feature.qualifiers['translation'])==1
                counter += 1

                geneinfo = {}
                #get aa info
                geneinfo.update(info_from_feature(feature))
                #get na info
                record = feature.extract(gb_seqrec)
                geneinfo.update(info_from_record(record))
                #build output with custom fields (and add them to info list)
                geneid = "fig|" + str(orginfo["taxon"]) + "." + str(version_number) + ".peg." + str(counter)
                geneinfo["feature_id"] = geneid
                geneinfo["contig_id"] = orginfo["id"]
                geneinfo["source_description"] = orginfo["gb_description"]
                if "protein_id" in feature.qualifiers:
                    geneinfo["location"] = feature.qualifiers["protein_id"][0]
                else:
                    geneinfo["location"] = ""
                genename = geneinfo["aliases"]
                genedesc = geneinfo["function"] + " " + orginfo["gb_description"]
                geneinfo["gene_description"] = genedesc               

                # Add locus tag and existing feature_id to list of aliases
                aliases = []
                if "protein_id" in feature.qualifiers:
                    aliases.append(feature.qualifiers["protein_id"][0])
                if "locus_tag" in feature.qualifiers:
                    aliases.append(feature.qualifiers["locus_tag"][0])
                if "gene" in feature.qualifiers:
                    aliases.append(feature.qualifiers["gene"][0])
                geneidToAlias[geneid] = aliases
                genes.append(geneinfo)

    return orginfo, genes, geneidToAlias

def findtypes(info):
    return dict([(k, type(v)) for k, v in info.items()])

#Entrez requires an email, if not set namually, this guesses one
email = None
guessemail = getpass.getuser() + '@' + socket.getfqdn()
if email == None:
    email = guessemail
Entrez.email = email

if __name__ == '__main__':
    usage="%prog [options] -g genbank_file"
    description='''    The purpose of this script is to take in a Genbank file (.gbk) 
    with multiple contig genbank files concatinated and 
    automatically organize three pieces of information required for input into ITEP:

    1: A genbank file with the coorrect name reflecting the taxID and a version number
    2: A tab-delimited file with the required format (see documentaion in the header of this file)
    3: An augmented aliases file containing locus tag and geneId information from the genbank file if available.

    You should NOT put your original genbank file in the $ROOT/genbank/ folder - instead put it somewhere else and call this
    script. This script will automatically reformat it and put the reformatted file into the /genbank/ folder.

    The genbank file is placed in $ROOT/genbank/, the tab-delimited file in $ROOT/raw/ and the aliases
    file is appended to any existing aliase file in the location $ROOT/aliases/aliases

    '''
    parser = optparse.OptionParser(usage=usage, description=description)
    parser.add_option("-g", "--genbank_file", help="Input genbank file concatinated across contigs (REQUIRED)", action="store", type="str", dest="genbank_file", default=None)
    parser.add_option("-o", "--org_file", help="(OPTIONAL) a file to which to dump organism data.", action="store", type="str", dest="org_file", default=None)
    parser.add_option("-r", "--replace", 
                      help="If specified, replaces old data related to the derived organism ID with the new data (D: Throws an error if an organism already is present with the derived organism ID)",
                      action="store_true", dest="replace", default=False)
    parser.add_option("-v", "--version_number", help="The second number in the \d+\.\d+ format of the organism ID - use this to distinguish between multiple genbank files with the same taxID (D:88888)",
                      action="store", type="int", dest="version_number", default=88888)
    (options, args) = parser.parse_args()
    
    if options.genbank_file is None:
        sys.stderr.write("ERROR: Genbank_file (-g) is a required argument\n")
        exit(2)

    orginfo, genes, aliases = genbank_extract(options.genbank_file, options.version_number)

    rootdir = locateRootDirectory()
    organism_id = str(orginfo["taxon"]) + "." + str(options.version_number)
    geneout_filename = os.path.join(rootdir, "raw", "%s.txt" %(organism_id))
    genbank_filename = os.path.join(rootdir, "genbank", "%s.gbk" %(organism_id))
    alias_filename = os.path.join(rootdir, "aliases", "aliases")

    # Try to prevent conflics between multiple organisms with the same taxID
    if os.path.exists(geneout_filename):
        if options.replace:
            # Note - to make this really robusst I should probably add a timestamp or something...
            bkgeneout_filename = os.path.join(rootdir, "%s.txt.bk" %(organism_id))
            sys.stderr.write("""WARNING: Backing up original gene output file %s to location %s in case something went wrong\n""" %(geneout_filename, bkgeneout_filename))
            shutil.copyfile(geneout_filename, bkgeneout_filename)
        else:
            sys.stderr.write("""ERROR: Gene output file %s already exists!
This could indicate a conflict in taxIDs between multiple organisms. Use -r to override this error and 
replace the existing file with a new one and remove the existing file, or use a different version number (-v)
if the genomes are really different.\n""" %(geneout_filename))
            exit(2)

    if os.path.exists(genbank_filename) and not options.replace:
        if options.replace:
            # Note - to make this really robusst I should probably add a timestamp or something...
            bkgenbank_filename = os.path.join(rootdir, "%s.txt.bk" %(organism_id))
            sys.stderr.write("""WARNING: Backing up original gene output file %s to location %s in case something went wrong\n""" %(genbank_filename, bkgenbank_filename))
            shutil.copyfile(geneout_filename, bkgenbank_filename)
        else:
            sys.stderr.write("""ERROR: Genbank output file %s already exists!
This could indicate a conflict in taxIDs between multiple organisms. Use -r to override this error and 
replace the existing file with a new one and remove the existing file, or use a different version number (-v
if the genomes are really different.\n""" %(genbank_filename))
            exit(2)

    if os.path.exists(alias_filename):
        # Lets check to see if we can find our organism in the existing alias file
        bk_file = os.path.join(rootdir, "aliases", "aliases.bk")
        shutil.copyfile(alias_filename, bk_file)
        alias_ptr = open(alias_filename, "w")
        # The entire organism ID must match and we must avoid subsets.
        tosearch = "|" + organism_id + "."
        for line in open(bk_file, "r"):
            if tosearch in line:
                if not options.replace:
                    sys.stderr.write("""ERROR: Entries with gene ID %s already exist in the aliases file.
This could indicate a conflict in taxIDs between multiple organisms. Use -r to override this error and
replace the existing entries in the aliases file with new entries or use a different version number (v) if the
genomes are really different.""")
                    alias_ptr.close()
                    os.rename(bk_file, alias_filename)
                    exit(2)
                else:
                    continue
            else:
                alias_ptr.write(line)
        alias_ptr.close()

    # Lets re-write the genbank file first.
    shutil.copyfile(options.genbank_file, genbank_filename)

    # Now we generate a tab-delimited file with the following fields:
    names = ["contig_id",          #from gi of gb file
    "feature_id",         #taxonID from NCBI lookup of accession, geneID from the db_xref GI in feature
    "type",               #CDS = peg
    "location",           #ignored by iTEP, source accession here
    "start",              #+1 index and rev if on neg strand
    "stop",               #+1 index and rev if on neg strand
    "strand",             #converted to + or -
    "function",           #from feature product
    "aliases",            #ignored by iTEP, from gene accession
    "figfam",             #ignored by iTEP, empty for all
    "evidence_codes",     #ignored by iTEP, empty for all
    "nucleotide_sequence",#as recorded in feature
    "aa_sequence"]        #as recorded in feature, not translated manually

    geneout_file = open(geneout_filename, 'w')
    geneout = csv.DictWriter(geneout_file, fieldnames = names, delimiter="\t")
        #geneout.writeheader()

    for gene in genes:
        geneout.writerow(dict([(n, gene[n]) for n in names]))
    geneout_file.close()
    sys.stderr.write("Text file saved as %s\n" % geneout_filename)

    # IMPORTANT: Append, don't use "w" here (we don't want to blow away all the different organism entries...)
    alias_file = open(alias_filename, "a")
    for geneid in aliases:
        for alias in aliases[geneid]:
            alias_file.write("%s\t%s\n" %(geneid, alias))
    alias_file.close()

    #Write out the organism data
    if options.org_file is not None:
        names = orginfo.keys()
        names.sort()
        orgout_file = open(options.orgout_file, "w")
        orgout = csv.DictWriter(orgout_file, fieldnames = names)
        orgout.writeheader()
        for org in orginfos:
            orgout.writerow(org)
        orgout_file.close()
