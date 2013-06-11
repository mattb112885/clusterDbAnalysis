
/* Contains the RAW tables */

CREATE TABLE rawdata(
       "contig" TEXT,
       "geneid" TEXT PRIMARY KEY,
       "ftype" TEXT,
       "location" TEXT,
       "genestart" INTEGER,
       "geneend" INTEGER,
       "strand" TEXT,
       "annotation" TEXT,
       "aliases" TEXT,
       "figfam" TEXT,
       "evidence" TEXT,
       "nucseq" TEXT,
       "aaseq" TEXT
       );

/* BLASTP results */
CREATE TABLE blastresults(
       "querygene" TEXT,
       "targetgene" TEXT,
       "pctid" REAL,
       "alnlen" INTEGER,
       "mismatches" INTEGER,
       "gapopens" INTEGER,
       "querystart" INTEGER,
       "queryend" INTEGER,
       "substart" INTEGER,
       "subend" INTEGER,
       "evalue" REAL,
       "bitscore" REAL,
       FOREIGN KEY(querygene) REFERENCES rawdata(geneid),
       FOREIGN KEY(targetgene) REFERENCES rawdata(geneid)
       );

CREATE TABLE blastn_results(
       "querygene" TEXT,
       "targetgene" TEXT,
       "pctid" REAL,
       "alnlen" INTEGER,
       "mismatches" INTEGER,
       "gapopens" INTEGER,
       "querystart" INTEGER,
       "queryend" INTEGER,
       "substart" INTEGER,
       "subend" INTEGER,
       "evalue" REAL,
       "bitscore" REAL,
       FOREIGN KEY(querygene) REFERENCES rawdata(geneid),
       FOREIGN KEY(targetgene) REFERENCES rawdata(geneid)
       );

CREATE TABLE organisms(
       "organism" TEXT PRIMARY KEY,
       "organismid" TEXT UNIQUE
       );

/* I decided to reduce our memory footprint by just creating a view for
   the organism - gene links */
CREATE TABLE geneinfo(
       "geneid" TEXT,
       "organismid" TEXT,
       "organism" TEXT,
       "organismabbrev" TEXT,
       "contig_mod" TEXT,
       "strandsign" INTEGER,
       "aalen" INTEGER,
       "nuclen" INTEGER,
       FOREIGN KEY(geneid) REFERENCES rawdata(geneid),
       FOREIGN KEY(organismid) REFERENCES organisms(organismid),
       FOREIGN KEY(organism) REFERENCES organisms(organism),
       );

/* Note - distance is number of genes */
CREATE TABLE neighborhoods(
       "centergene" TEXT,
       "neighborgene" TEXT,
       "distance" TEXT,
       "contig_mod" TEXT,
       "startloc" INTEGER,
       "stoploc" INTEGER,
       "strand" TEXT,
       "annotation" TEXT,
       FOREIGN KEY(centergene) REFERENCES rawdata(geneid),
       FOREIGN KEY(neighborgene) REFERENCES rawdata(geneid)
       );
       
.separator "\t"
.import organisms organisms
.import db/raw_cat rawdata
.import db/mod_cat geneinfo
.import db/blastres_cat blastresults
.import db/blastnres_cat blastn_results
.import db/neighborhoods neighborhoods

/* Removing these might be slower on large data sets. Needs testing.

It was ~10% faster on a set of 4 genomes.

CREATE INDEX blastqueryidx ON blastresults (querygene);
CREATE INDEX blasttargetidx ON blastresults (targetgene);

CREATE INDEX blastnqueryidx ON blastn_results (querygene);
CREATE INDEX blastntargetidx ON blastn_results (targetgene);

*/

CREATE INDEX blastidx ON blastresults(querygene, targetgene);
CREATE INDEX blastnidx ON blastn_results(querygene, targetgene);

CREATE INDEX neighborhoodcenteridx ON neighborhoods (centergene);

/* Add self-bit score to blastp results table */

CREATE TABLE blast_self AS 
       SELECT * FROM blastresults 
       INNER JOIN (
       	     SELECT querygene, targetgene, max(bitscore) AS bitscore FROM blastresults
	     WHERE querygene = targetgene
	     GROUP BY querygene
	     ) AS s
	     ON s.querygene = blastresults.querygene AND s.bitscore = blastresults.bitscore
	     WHERE blastresults.querygene = s.querygene
	     AND blastresults.targetgene = s.targetgene
	     AND blastresults.querygene = blastresults.targetgene;


CREATE TABLE blastn_self AS 
       SELECT * FROM blastn_results 
       INNER JOIN (
       	     SELECT querygene, targetgene, max(bitscore) AS bitscore FROM blastn_results 
	     WHERE querygene = targetgene
	     GROUP BY querygene
	     ) AS s
	     ON s.querygene = blastn_results.querygene
	     WHERE blastn_results.querygene = s.querygene
	     AND blastn_results.targetgene = s.targetgene
	     AND blastn_results.querygene = blastn_results.targetgene
	     AND s.bitscore = blastn_results.bitscore;

/* We need to index the self-bit scores.
   These SHOULD be unique ... */
CREATE UNIQUE INDEX selfqueryidx ON blast_self(querygene);
CREATE UNIQUE INDEX blastnselfqueryidx ON blastn_self(querygene);

CREATE VIEW s AS
       SELECT blastresults.*, blast_self.bitscore AS queryselfbit
       FROM blastresults
       INNER JOIN blast_self ON blast_self.querygene = blastresults.querygene;

CREATE VIEW sn AS
       SELECT blastn_results.*, blastn_self.bitscore AS queryselfbit
       FROM blastn_results
       INNER JOIN blastn_self ON blastn_self.querygene = blastn_results.querygene;

CREATE TABLE blastres_selfbit AS
       SELECT s.*, blast_self.bitscore AS targetselfbit FROM s
       INNER JOIN blast_self ON blast_self.querygene = s.targetgene;


CREATE TABLE blastnres_selfbit AS
       SELECT sn.*, blastn_self.bitscore AS targetselfbit FROM sn
       INNER JOIN blastn_self ON blastn_self.querygene = sn.targetgene;

CREATE INDEX selfbitqueryidx ON blastres_selfbit(querygene);
CREATE INDEX selfbittargetidx ON blastres_selfbit(targetgene);

CREATE INDEX selfbitquerytarget ON blastres_selfbit(querygene, targetgene);

CREATE INDEX blastn_selfbitqueryidx ON blastnres_selfbit(querygene);
CREATE INDEX blastn_selfbittargetidx ON blastnres_selfbit(targetgene);

DROP VIEW s;
DROP VIEW sn;

/* Processed table (this is the table you should generally query against - it results in the 'geneinfo' tables):
   Gene ID | organism | organismID | organism abbreviation | contig ID | gene start | gene end | strand | annotation | nucleotide sequence | amino acid sequence 
 */

CREATE TABLE processed AS
       SELECT geneinfo.geneid, organisms.organism, organisms.organismid, organisms.organismid AS placeholder, geneinfo.contig_mod,  
       	      rawdata.genestart, rawdata.geneend, rawdata.strand, geneinfo.strandsign, 
              rawdata.annotation, rawdata.nucseq, rawdata.aaseq
       FROM organisms
       INNER JOIN geneinfo ON geneinfo.organismid = organisms.organismid
       INNER JOIN rawdata ON geneinfo.geneid = rawdata.geneid;

CREATE UNIQUE INDEX processedgeneids ON processed(geneid);
CREATE INDEX processedcontigs ON processed(contig_mod);
CREATE INDEX processedorganismids ON processed(organismid);

/* These tables are no longer needed. Drop them to reduce the memory footprint and reduce confusion. */
DROP TABLE blastresults;
DROP TABLE blastn_results;
DROP TABLE geneinfo;

/* This command is needed to actually get the DB to shrink upon command.
It takes a while to run and requires about twice the size of the database while it's running,
but after it's done it can save hundreds of gigabytes.

The alternative is to do a PRAGMA auto_vacuum FULL; ... or just to keep the tables around. */
VACUUM;

/* This should make some queries faster and doesn't take much time. */
ANALYZE;