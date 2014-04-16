DROP TABLE IF EXISTS clusters;
DROP VIEW IF EXISTS clusterorgs;
DROP TABLE IF EXISTS distinctorgs;

/* This is intended to be run after builddb_1.sql
   The two scripts are separated in order to allow us to generate clusters
   based on querying the blast results from the database. */
CREATE TABLE clusters(
       "runid" TEXT,
       "clusterid" INTEGER,
       "geneid" TEXT,
       FOREIGN KEY(geneid) REFERENCES rawdata(geneid)
       );

.separator "\t"
.import db/flat_clusters clusters

CREATE INDEX clusteridx ON clusters (geneid);
CREATE INDEX clusterrunidx ON clusters (runid) ;
CREATE INDEX clusterididx ON clusters (clusterid);

CREATE INDEX clusteridclusterrunidx ON clusters(runid, clusterid);

/* Clusters adjacent to what organism they belong to */
CREATE VIEW clusterorgs AS
       SELECT clusters.*, processed.organism FROM clusters
       INNER JOIN processed ON processed.geneid = clusters.geneid
       ORDER BY clusters.clusterid DESC;

/* Run IDs and what organisms are part of those run IDs */
CREATE TABLE distinctorgs AS SELECT DISTINCT runid, organism FROM clusterorgs;

ANALYZE;