/* This script generates two tables
 The first one is "contigs" which contains the genome sequence for every
 contig (scaffolds). It is used to generate TBLASTN results against a 
 specified list of organisms on the fly.

 The second one is a table meant for holding TBLASTN hits that are inconsistent with
 or missing from the current annotation. */ 

DROP TABLE IF EXISTS contigs;
DROP TABLE IF EXISTS tblastn;

/* You will need to make the seq bigger if you have very large genomes.
I chose 10,000,000 because this is the size of the largetst Bacterial genome
in the CBS Genome Atlas except for 5. I will throw an error if your contigs are
larger than this. If they are much smaller you can save space by reducing 
this number. */
CREATE TABLE contigs(
       "contig_mod" TEXT,
       "seq" TEXT,
       "organismid" TEXT,
       FOREIGN KEY(organismid) REFERENCES organisms(organismid)
);

.separator "\t"
.import db/contigs contigs

/* If this one isn't unique we're in trouble */
CREATE UNIQUE INDEX contigcontigs ON contigs(contig_mod);
CREATE INDEX contigorgs ON contigs(organismid);

CREATE TABLE tblastn(
       "queryid" TEXT,
       "querylen" INTEGER,
       "targetcontig" TEXT,
       "targetorganism" TEXT,
       "tblaststart" INTEGER,
       "tblastend" INTEGER,
       "tblastlen" INTEGER,
       "queryoverlappct" REAL,
       "evalue" REAL,
       "bitscore" REAL,
       "hitframe" INTEGER,
       "strandedstring" TEXT,
       "targetgeneid" TEXT,
       "targetannotation" TEXT,
       "targetgenelen" INTEGER,
       "targetoverlappct" REAL,
       "tblastn_id" TEXT,
       FOREIGN KEY(targetorganism) REFERENCES organisms(organism),
       FOREIGN KEY(queryid) REFERENCES processed(geneid),
       FOREIGN KEY(targetcontig) REFERENCES contigs(contig_mod)
);

ANALYZE;