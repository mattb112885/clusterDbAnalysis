#!/usr/bin/env python

'''
This file contains functions intended to deal with modifications to GENBANK files
for increased compatibility between ITEP and other tools.

'''

from __future__ import print_function
import sys

def addItepGeneIdsToGenbank(multi_gbk_object, tbl, ignoreseq=False):
    '''
    Given a RAW table (list of tuples) and a genbank object,
    attempt to match up every object in the RAW table with every
    gene in the table and add the ITEP IDs as a db_xrer (called ITEP)
    to the genbank file.
    '''
    # Index search on the table
    contigidx = 0
    ididx = 1
    startidx = 4
    stopidx = 5
    seqidx = 11
    indexed_array = {}
    for ln in tbl:
        if ln[startidx] == 'start':
            continue

        # Do both of these because of possible incompatibilities in how the orderinf of
        # start and stop works...
        contig = ln[contigidx]
        start = int(ln[startidx])
        stop = int(ln[stopidx])
        seq = ln[seqidx]
        myid = ln[ididx]
        indexed_array[(contig, start, stop)] = (seq, myid)
        indexed_array[(contig, stop, start)] = (seq, myid)

    # First we need this to not be a generator so we can actually modify the darn thing
    multi_gbk_object = list(multi_gbk_object)

    # I match up by the following:
    # Start location
    # Stop location
    # DNA sequence (due to possible differences in translation)
    # If these all match with a given element of the table then its good.
    newToOriginalName = {}
    for ii in range(len(multi_gbk_object)):
        # This needs to match what we do to make the raw file table in the first place.
        if multi_gbk_object[ii].id == "unknown":
            original_name = multi_gbk_object[ii].name
        else:
            original_name = multi_gbk_object[ii].id
        # We make a temporary contig ID that is less than the 16 character limit from Biopython.
        # This will be OK as long as the number of contigs is less than 10 million or so.
        new_name = "REP_CTG%08d" %(ii)
        multi_gbk_object[ii].name = new_name
        newToOriginalName[new_name] = original_name

        # Add ITEP IDs to the genbank files.
        for jj in range(len(multi_gbk_object[ii].features)):
            if multi_gbk_object[ii].features[jj].type == "source":
                if "db_xref" in multi_gbk_object[ii].features[jj].qualifiers:
                    multi_gbk_object[ii].features[jj].qualifiers["db_xref"].append("originalContig:%s" %(original_name))
                else:
                    multi_gbk_object[ii].features[jj].qualifiers["db_xref"] = [ "originalContig:%s" %(original_name) ]
            # We don't want to modify things that aren't coding sequences...
            if multi_gbk_object[ii].features[jj].type != "CDS":
                continue
            record = multi_gbk_object[ii].features[jj].extract(multi_gbk_object[ii])
            seq = str(record.seq)
            # Get locations to match up (note they are splice indexes and start from 0 hence the +1 to start and nothing added to end)
            location = multi_gbk_object[ii].features[jj].location
            featurestart = int(location.start) + 1
            featureend = int(location.end)
            querytup = (original_name, featurestart, featureend)
            if querytup in indexed_array and ( ignoreseq or seq.lower() == indexed_array[querytup][0].lower() ):
                if "db_xref" in multi_gbk_object[ii].features[jj].qualifiers:
                    alreadyHasItep = False
                    itep_id = "ITEP:%s" %(indexed_array[querytup][1])
                    for ref in multi_gbk_object[ii].features[jj].qualifiers["db_xref"]:
                        if itep_id in ref:
                            alreadyHasItep = True
                            break
                    if not alreadyHasItep:
                        multi_gbk_object[ii].features[jj].qualifiers["db_xref"].append("ITEP:%s" %(indexed_array[querytup][1]))
                else:
                    multi_gbk_object[ii].features[jj].qualifiers["db_xref"] = [ "ITEP:%s" %(indexed_array[querytup][1]) ]
            else:
                sys.stderr.write("BAD - the gene in the following location did not have a match in the table file!\n")
                sys.stderr.write(str(location) + "\n")
                sys.stderr.write(str(seq) + "\n\n")
                if querytup in indexed_array:
                    sys.stderr.write("did not match sequence in table file in the same location:\n")
                    sys.stderr.write(str(indexed_array[querytup][0]) + "\n\n")
            pass
        pass
    return multi_gbk_object, newToOriginalName

def replaceTemporaryIdsWithOriginalIds(genbank_fid, newToOriginalName, output_fid):
    '''
    (Intended for internal use)
    
    Given file objects for a genbank file and an output file, replaces the temporary IDs
    with the original ones...

    This is a workaround for the biopython limit of 16 characters for writing contig IDs
    in Genbank files.
    '''

    for line in genbank_fid:
        if line.startswith("LOCUS"):
            for new in newToOriginalName:
                if new in line:
                    line = line.replace(new, newToOriginalName[new])
                    break
                pass
            pass
        output_fid.write(line)

    return
