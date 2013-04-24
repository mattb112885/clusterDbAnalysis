/* RPSBLAST-related db functions */

DROP TABLE IF EXISTS external_clusters;
DROP TABLE IF EXISTS rpsblast_results;

/* This table is cddid.tbl from NCBI */
CREATE TABLE external_clusters(
       "counter" INTEGER,
       "external_clusterid" TEXT,
       "clustername" TEXT,
       "description" TEXT,
       "profilelength" INTEGER
);

CREATE TABLE rpsblast_results(
       "querygene" TEXT,
       "external_clusterid" TEXT,
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
       FOREIGN KEY(external_clusterid) REFERENCES external_clusters(external_clusterid)
);

.separator "\t"

.import db/cddid.tbl external_clusters
.import db/external_CDD rpsblast_results

CREATE INDEX externalclusterididx ON external_clusters (external_clusterid);
CREATE INDEX rpsblastquery ON rpsblast_results (querygene);
CREATE INDEX rpsblastexternalclusterididx ON rpsblast_results (external_clusterid);