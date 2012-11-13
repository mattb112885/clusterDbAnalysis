/* RPSBLAST-related db functions */

DROP TABLE IF EXISTS external_clusters;
DROP TABLE IF EXISTS rpsblast_results;

/* This table is cddid.tbl from NCBI */
CREATE TABLE external_clusters(
       "counter" INT,
       "external_clusterid" VARCHAR(36),
       "clustername" VARCHAR(36),
       "description" VARCHAR(1028),
       "profilelength" INT
);

CREATE TABLE rpsblast_results(
       "querygene" VARCHAR(128),
       "external_clusterid" VARCHAR(128),
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
       FOREIGN KEY(external_clusterid) REFERENCES external_clusters(external_clusterid)
);

.separator "\t"

.import db/cddid.tbl external_clusters
.import db/external_CDD rpsblast_results

CREATE INDEX externalclusterididx ON external_clusters (external_clusterid);
CREATE INDEX rpsblastquery ON rpsblast_results (querygene);
CREATE INDEX rpsblastexternalclusterididx ON rpsblast_results (external_clusterid);