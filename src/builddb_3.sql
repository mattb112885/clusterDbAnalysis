/* This script generates two tables
 The first one is "contigs" which contains the genome sequence for every
 contig (scaffolds). It is used to generate TBLASTN results against a 
 specified list of organisms on the fly.

 The second one is a table meant for holding TBLASTN hits that are inconsistent with
 or missing from the current annotation. */ 

DROP TABLE IF EXISTS contigs;
DROP TABLE IF EXISTS tblastn;

CREATE TABLE contigs(
       "contig_mod" VARCHAR(256),
       "seq" VARCHAR(8000000),
       "organismid" VARCHAR(256),
       FOREIGN KEY(organismid) REFERENCES organisms(organismid)
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

CREATE INDEX contigcontigs ON contigs(contig_mod);
