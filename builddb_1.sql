
/* I found that the maximum length of proteins was 12,549 so 15,000 should be good.
   Command:
   cat db/raw_cat | python src/raw2processed.py organisms | cut -f 8 | sort -n | less */

CREATE TABLE rawdata(
       "contig" VARCHAR(32),
       "geneid" VARCHAR(32) PRIMARY KEY,
       "ftype" VARCHAR(8),
       "location" VARCHAR(32),
       "genestart" INT,
       "geneend" INT,
       "strand" VARCHAR(4),
       "annotation" VARCHAR(2048),
       "aliases" VARCHAR(2048),
       "figfam" VARCHAR(128),
       "evidence" VARCHAR(128),
       "nucseq" VARCHAR(45000),
       "aaseq" VARCHAR(15000)
       );

CREATE TABLE blastresults(
       "querygene" VARCHAR(128),
       "targetgene" VARCHAR(128),
       "pctid" FLOAT,
       "alnlen" INT,
       "mismatches" INT,
       "gapopens" INT,
       "querystart" INT,
       "queryend" INT,
       "substart" INT,
       "subend" INT,
       "evalue" FLOAT,
       "bitscore" FLOAT,
       FOREIGN KEY(querygene) REFERENCES rawdata(geneid),
       FOREIGN KEY(targetgene) REFERENCES rawdata(geneid)
       );

/* abbreviation must be less than 6 characters for PHYLIP format to work correctly... */
CREATE TABLE organisms(
       "organism" varchar(128) PRIMARY KEY,
       "organismabbrev" varchar(6) UNIQUE,
       "organismid" varchar(128) UNIQUE
       );

/* I decided to reduce our memory footprint by just creating a view for
   the organism - gene links */
CREATE TABLE geneinfo(
       "geneid" VARCHAR(32),
       "organismid" VARCHAR(128),
       "organism" VARCHAR(128),
       "organismabbrev" varchar(6),
       "contig_mod" VARCHAR(256),
       "strandsign" INT,
       "aalen" INT,
       "nuclen" INT,
       FOREIGN KEY(geneid) REFERENCES rawdata(geneid),
       FOREIGN KEY(organismid) REFERENCES organisms(organismid),
       FOREIGN KEY(organism) REFERENCES organisms(organism),
       FOREIGN KEY(organismabbrev) REFERENCES organisms(organismabbrev),
       CHECK (aalen < 15000)
       );

.separator "\t"
.import organisms organisms
.import db/raw_cat rawdata
.import db/mod_cat geneinfo
.import db/blastres_cat blastresults

CREATE INDEX blastqueryidx ON blastresults (querygene);
CREATE INDEX blasttargetidx ON blastresults (targetgene);

/* Add self-bit score */
CREATE TEMPORARY TABLE blast_self AS
       SELECT * FROM blastresults 
       WHERE blastresults.querygene = blastresults.targetgene;

CREATE INDEX selfqueryidx ON blast_self(querygene);

CREATE VIEW s AS
       SELECT blastresults.*, blast_self.bitscore
       FROM blastresults
       INNER JOIN blast_self ON blast_self.querygene = blastresults.querygene;

CREATE TABLE blastres_selfbit AS
       SELECT s.*, blast_self.bitscore FROM s
       INNER JOIN blast_self ON blast_self.querygene = s.targetgene;

CREATE INDEX selfbitqueryidx ON blastres_selfbit(querygene);
CREATE INDEX selfbittargetidx ON blastres_selfbit(targetgene);

/* Resulting table:
   Gene ID | organism | organismID | organism abbreviation | contig ID | gene start | gene end | strand | annotation | nucleotide sequence | amino acid sequence 
 */

CREATE TABLE processed AS
       SELECT geneinfo.geneid, organisms.organism, organisms.organismid, organisms.organismabbrev, geneinfo.contig_mod,  
       	      rawdata.genestart, rawdata.geneend, rawdata.strand, geneinfo.strandsign, 
              rawdata.annotation, rawdata.nucseq, rawdata.aaseq
       FROM organisms
       INNER JOIN geneinfo ON geneinfo.organismid = organisms.organismid
       INNER JOIN rawdata ON geneinfo.geneid = rawdata.geneid;

CREATE INDEX processedgeneids ON processed(geneid);
CREATE INDEX processedcontigs ON processed(contig_mod);