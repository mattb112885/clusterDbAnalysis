/* This table is intended to contain information about
user-added genes, clusters, etc.

*/

DROP TABLE IF EXISTS user_genes;

/*

Note - There are a couple unique things about this table.

1: A lot of the fields can be empty. For example,
   we could have no contig ID (e.g. if we are looking at adding
   genes from raw reads or biochemical data), no sequence,
   no location, or no run/clusterID association.

2: This table is the responsibility of the user to generate (unlike
   all the others which are at least partly automatically-generated)
   and then import.

*/
CREATE TABLE user_genes (
       user_geneid TEXT PRIMARY KEY,
       organismid TEXT NOT NULL,
       genetype TEXT NOT NULL,
       contigid TEXT,
       startloc INTEGER,
       stoploc INTEGER,
       runid TEXT,
       clusterid INTEGER,
       seq TEXT,
       annotation TEXT,
       FOREIGN KEY (organismid) REFERENCES organisms(organismid)
);

CREATE INDEX user_gene_clusterrun ON user_genes(runid, clusterid);
/* CREATE INDEX user_gene_organisms ON user_genes(organismid); */

.separator "\t"
.import userdata/user_genes user_genes

UPDATE user_genes SET contigid=NULL WHERE contigid="";
UPDATE user_genes SET startloc=NULL WHERE startloc="";
UPDATE user_genes SET stoploc=NULL WHERE stoploc="";
UPDATE user_genes SET runid=NULL WHERE runid="";
UPDATE user_genes SET clusterid=NULL WHERE clusterid="";
UPDATE user_genes SET seq=NULL WHERE seq="";

ANALYZE;