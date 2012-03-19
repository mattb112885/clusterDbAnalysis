DROP TABLE clusters;
DROP VIEW clusterorgs;
DROP VIEW distinctorgs;

/* This is intended to be run after builddb_1.sql
   The two scripts are separated in order to allow us to generate clusters
   based on querying the blast results from the database. */
CREATE TABLE clusters(
       "runid" VARCHAR(128),
       "clusterid" INT,
       "geneid" VARCHAR(32),
       FOREIGN KEY(geneid) REFERENCES rawdata(geneid)
       );

.separator "\t"
.import db/flat_clusters clusters

CREATE INDEX clusteridx ON clusters (geneid);

/* Clusters adjacent to what organism they belong to */
CREATE VIEW clusterorgs AS
       SELECT clusters.*, processed.organism FROM clusters
       INNER JOIN processed ON processed.geneid = clusters.geneid
       ORDER BY clusters.clusterid DESC;

/* Cluster IDs and what organisms are part of those cluster IDs */
CREATE VIEW distinctorgs AS SELECT DISTINCT clusterorgs.runid, clusterorgs.organism FROM clusterorgs;