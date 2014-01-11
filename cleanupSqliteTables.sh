#!/bin/bash

if [ $# -lt 1 ]; then
    echo "usage: ./cleanupSqliteTables.sh TRUE"
    echo ""
    echo "Description: When run (with TRUE as an argument),"
    echo "performs a VACUUM on the sqlite database to reclaim the"
    echo "hard drive space taken by dropped tables."
    echo ""
    echo "WARNING: Running this command requires you to have up to twice the disk"
    echo "space of the original database free on the partition containing /tmp/."
    echo "See https://www.sqlite.org/lang_vacuum.html for details."
    echo ""
    echo "It is recommended to run this command in a UNIX screen, as it can take a long time."
    echo ""
    exit 1
fi

if [ "$1" = "TRUE" ]; then
    echo "VACUUM;" | sqlite3 db/DATABASE.sqlite
fi