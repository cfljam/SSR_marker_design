#!/bin/sh
#design primer sets from MISA output using Primer 3 and MISA helper scripts
#USAGE sh design_MISA.sh <misa_file> <source_fasta_file> <output_file>

#get directory
SCRIPT=`readlink -f $0`
SCRIPTPATH=`dirname $SCRIPT`


perl $SCRIPTPATH/p3_in.pl  $1 $2  temp.p3in

cat  temp.p3in | primer3_core --io_version=3 >  temp.p3out


perl $SCRIPTPATH/p3_out.pl  temp.p3out $1  $3

rm -f temp.*

