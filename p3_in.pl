#!/usr/bin/perl -w
# Author: Thomas Thiel
# Program name: primer3_in.pl
# Description: creates a PRIMER3 input file based on SSR search results

open (IN,"<$ARGV[0]") || die ("\nError: Couldn't open misa.pl results file (*.misa) !\n\n");

#my $filename = $ARGV[0];
#$filename =~ s/\.misa//;
open (SRC,"<$ARGV[1]") || die ("\nError: Couldn't open source file containing original FASTA sequences !\n\n");
open (OUT,">$ARGV[2]");

undef $/;
$in = <IN>;
study $in;

$/= ">";

#my $count;
while (<SRC>)
  {
  next unless (my ($id,$seq) = /(.*?)\n(.*)/s);
  $seq =~ s/[\d\s>]//g;#remove digits, spaces, line breaks,...
  while ($in =~ /$id\t(\d+)\t\S+\t\S+\t(\d+)\t(\d+)/g)
    {
    my ($ssr_nr,$size,$start) = ($1,$2,$3);
    #$count++;
    print OUT "PRIMER_SEQUENCE_ID=$id"."_$ssr_nr\nSEQUENCE=$seq\n";
    print OUT "PRIMER_PRODUCT_SIZE_RANGE=100-280\n";
    print OUT "TARGET=",$start-3,",",$size+6,"\n";
    print OUT "PRIMER_MAX_END_STABILITY=250\n=\n"
    };
  };
#print "\n$count records created.\n";
