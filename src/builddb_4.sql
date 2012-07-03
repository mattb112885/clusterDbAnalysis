CREATE TABLE IF NOT EXISTS clusterlists (
        "alignmentid" INTEGER,
	"runid" VARCHAR(256),
	"clusterid" INTEGER,
	FOREIGN KEY(runid,clusterid) REFERENCES clusters(runid,clusterid),
	FOREIGN KEY(alignmentid) REFERENCES alignments(alignmentid)
	);
	
CREATE TABLE IF NOT EXISTS alignments (
	"alignmentid" INTEGER PRIMARY KEY,
	"alnmethod" VARCHAR(36),
	"trimmethod" VARCHAR(36)
	);

CREATE TABLE IF NOT EXISTS alignedseqs (
	"alignmentid" INTEGER,
	"geneid" INTEGER,
	"alignedaaseq" VARCHAR(60000),
	FOREIGN KEY(alignmentid) REFERENCES alignments(alignmentid),
	FOREIGN KEY(geneid) REFERENCES rawdata(geneid)
	);