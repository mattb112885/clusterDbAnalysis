/* This script generates two tables
 The first one is "contigs" which contains the genome sequence for every
 contig (scaffolds). It is used to generate TBLASTN results against a 
 specified list of organisms on the fly.

 The second one is a table meant for holding TBLASTN hits that are inconsistent with
 or missing from the current annotation. */ 
CREATE TABLE contigs(
       "contig_mod" VARCHAR(256),
       "seq" VARCHAR(8000000),
       FOREIGN KEY(contig_mod) REFERENCES processed(contig_mod)
);

CREATE TABLE tblastn(
       "queryid" VARCHAR(256),
       "targetcontig" VARCHAR(256),
       "tstart" INT,
       "tend" INT,
       "evalue" FLOAT,
       "bitscore" FLOAT,
       "strand" INT,
       "strandconsistency" VARCHAR(32),
       "targetgene" VARCHAR(32),
       "targetgenestart" INT,
       "targetgeneend" INT,
       "geneoverlap" INT,
       "geneoverlappct" FLOAT,
       "geneoverlapannote" VARCHAR(2048),
       FOREIGN KEY(targetcontig) REFERENCES processed(contig_mod),
       FOREIGN KEY(queryid) REFERENCES processed(geneid),
       FOREIGN KEY(targetgene) REFERENCES processed(geneid)
);

.separator "\t"
.import db/contigs contigs