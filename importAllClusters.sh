#!/bin/bash

########
# importAllClusters.sh
# This is a utility that imports all of the clusters by first flattening any clusters
# that exist in the "clusters" directory and then running the necessary SQL to create and
# populate the cluster data tables in the SQLite database.
#
# You can run this script to avoid having to run MCL and just import clusterings that
# are already generated from another source (see scripts/importExternalClustering.py
# for a convenient interface to this function)
#
# This function MUST be run as ./importAllClusters.sh (running it from another directory will fail).
########

# Modify the format of cluster files to be SQL-friendly.
# Each separate clustering is given a random (very highly likely unique) identifier
# which is subsequently used to identify which clustering solution we want from SQL
# given a set of organism names
cd clusters;
for file in *; do
    cat ${file} | flattenClusterFile.py -n ${file} > ../flatclusters/${file};
done
cd ..;

cat flatclusters/* > db/flat_clusters;

# Import clusters into the db
echo "Importing cluster information into database...";
sqlite3 db/DATABASE.sqlite < src/internal/builddb_2.sql;

# Once clusters are loaded, make a pre-built presence\absence table.
# This is a python script because it involves some operations that would be impossible or hard in pure SQL.
echo "Making a pre-built presence\absence table in the database..."
db_loadPresenceAbsence.py;

rm db/flat_clusters;