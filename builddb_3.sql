/* Import tblastn information */

CREATE TABLE tblastres(
       querygene VARCHAR(128),
       targetcontig VARCHAR(256),
       pctid FLOAT,
       alnlen INT,
       mismatches INT,
       gapopens INT,
       qstart INT,
       qend INT,
       targetstart INT,
       targetend INT,
       evalue FLOAT,
       bitscore FLOAT,
       FOREIGN KEY (querygene) REFERENCES processed(geneid),
       FOREIGN KEY(targetcontig) REFERENCES processed(contig_mod)
);

.separator "\t"
.import db/tblastn_cat tblastres

CREATE INDEX tblast_query ON tblastres(querygene);
CREATE INDEX tblast_contig ON tblastres(targetcontig);

/* Convert the TBLASTN results from each query [tblastres.querygene] 
   into a list of overlapping proteins [processed.geneid AS targetgene]

   If a TBLASTN hit overlaps with multiple protiens there will be multiple processed.geneid's for
   a given tblastres.querygene / tblastres.targetconfig pair.

   Note - all properties in this table except the query gene are properites of the TARGET gene - 
   the TBLASTN's estimation of the gene location, and the actual location of the containing gene (according to RAST)

   Note - this command takes a very long time (~1.5 hours on Cocoa)

   The purpose of processed.strandsign is to avoid missing hits that happen to be on the - strand...
   So I am pretty sure I was missing over half of the hits before thanks to that.

   This isn't strictly going to exclude things on different strands from corresponding with each other. However, it's likely
   that the same protein won't match a sequence on overlapping regions on opposite strands simultaneously...
   If that happens I'll be very surprised and kind of pissed off.

   Note - This is INDEPENDENT of the clustering results. It is critical that this be the case so that we can update the clusters
   without recalculating this table.

   Note - the < 9 is to enforce that the target hit can only be considered part of the query gene if the overlap is at least 9 base pairs
   (RAST sometimes annotates overlapping genes with less than 3 base pairs overlapping, and we don't want to pull all of those in)

*/

CREATE TABLE tblast_converted AS
         SELECT tblastres.querygene, 
              processed.geneid AS targetgene, processed.contig_mod, processed.organism AS targetorganism, processed.annotation AS targetannotation,
	      processed.genestart AS actual_targetstart, processed.geneend AS actual_targetend,
              tblastres.targetstart AS tblast_targetstart, tblastres.targetend AS tblast_targetend,
              tblastres.pctid, tblastres.evalue
         FROM tblastres 
         INNER JOIN processed
         ON NOT processed.strandsign*tblastres.targetend - processed.strandsign*processed.genestart < 9
         AND NOT processed.strandsign*processed.geneend - processed.strandsign*tblastres.targetstart < 9 
         AND tblastres.targetcontig = processed.contig_mod;

CREATE INDEX totest_targetgenes ON tblast_converted(targetgene);
CREATE INDEX totest_querygenes ON tblast_converted(querygene);