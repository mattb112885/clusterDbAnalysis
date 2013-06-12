#!/bin/sh

# Usage: ./main5.sh
#
# Load user-specifc data into the database

if [ ! -f /userfiles/user_genes ]; then
    echo "Main5.sh is intended for loading the user_genes file (located in userfules/)"
    echo "No user_genes file was found so this command is a no-op."
    echo ""
    exit 1;
fi

# Load user_genes into the database
sqlite3 db/DATABASE.sqlite < src/internal/builddb_5.sql;

# Re-load the presence-absence table with user-defined genes
db_loadPresenceAbsence.py;

