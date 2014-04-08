#!/bin/sh
## clean_fasta_header.sh 
##Remove descriptions from header

sed 's/\(>\w*\)\s*.*/\1/' 