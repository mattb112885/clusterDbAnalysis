#!/bin/bash

# Load user-specifc data into the database

if [ ! -f userdata/user_genes ]; then
    echo "setup_step5.sh is intended for loading the user_genes file (located at userdata/user_genes)"
    echo "No user_genes file was found so this command is a no-op."
    echo ""
    exit 1;
fi

# Make sure we are pointing at the right version of the repo
source SourceMe.sh

# Load user_genes into the database
echo "Loading user gene data..."
sqlite3 db/DATABASE.sqlite < src/internal/builddb_5.sql;

# Re-load the presence-absence table with user-defined genes
echo "Re-generating presence-absence table..."
db_loadPresenceAbsence.py;

