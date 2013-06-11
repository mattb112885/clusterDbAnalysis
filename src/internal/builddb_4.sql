/* RPSBLAST-related db functions */

DROP TABLE IF EXISTS external_clusters;
DROP TABLE IF EXISTS rpsblast_results;

/* This table is cddid.tbl from NCBI */
CREATE TABLE external_clusters(
       "cdd_id" INTEGER,
       "external_clusterid" TEXT,
       "clustername" TEXT,
       "description" TEXT,
       "profilelength" INTEGER
);

CREATE TABLE rpsblast_results(
       "querygene" TEXT,
       "cdd_id" INTEGER,
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
       FOREIGN KEY(cdd_id) REFERENCES external_clusters(cdd_id)
);

.separator "\t"

.import db/cddid.tbl external_clusters
.import db/external_CDD rpsblast_results

CREATE UNIQUE INDEX externalclusterididx ON external_clusters (external_clusterid);
CREATE UNIQUE INDEX externalclustercddid ON external_clusters (cdd_id);
CREATE INDEX rpsblastquery ON rpsblast_results (querygene);
CREATE INDEX rpsblastcddidx ON rpsblast_results (cdd_id);

ANALYZE;