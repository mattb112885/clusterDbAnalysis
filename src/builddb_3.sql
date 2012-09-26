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
       "querylen" INT,
       "targetcontig" VARCHAR(256),
       "targetorganism" VARCHAR(256),
       "tblaststart" INT,
       "tblastend" INT,
       "tblastlen" INT,
       "queryoverlappct" FLOAT,
       "evalue" FLOAT,
       "bitscore FLOAT,
       "hitframe" INT,
       "strandedstring" VARCHAR(16),
       "targetgeneid" VARCHAR(32),
       "targetannotation" VARCHAR(2048),
       "targetgenelen" INT,
       "targetoverlappct" FLOAT,
       FOREIGN KEY(targetorganism) REFERENCES organisms(organism),
       FOREIGN KEY(queryid) REFERENCES processed(geneid),
       FOREIGN KEY(targetcontig) REFERENCES contigs(contig_mod)
);

.separator "\t"
.import db/contigs contigs

CREATE INDEX contigcontigs ON contigs(contig_mod);
