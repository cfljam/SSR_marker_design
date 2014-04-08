#!/usr/bin/perl -w
# Author: Thomas Thiel
# Program name: misa.pl

###_______________________________________________________________________________
###
###Program name:p3_ misa_parameter.pl
###Author:       Thomas Thiel
###Release date: 14/12/01 (version 1.0)
###
###_______________________________________________________________________________
###
## _______________________________________________________________________________
##
## DESCRIPTION: Tool for the identification and localization of
##              (I)  perfect microsatellites as well as
##              (II) compound microsatellites (two individual microsatellites,
##                   disrupted by a certain number of bases)
##
## SYNTAX:   misa.pl <FASTA file>
##
##    <FASTAfile>    Single file in FASTA format containing the sequence(s).
##
##    In order to specify the search criteria, an additional file containing
##    the microsatellite search parameters is required named "misa.ini", which
##    has the following structure:
##      (a) Following a text string beginning with 'def', pairs of numbers are
##          expected, whereas the first number defines the unit size and the
##          second number the lower threshold of repeats for that specific unit.
##      (b) Following a text string beginning with 'int' a single number defines
##          the maximal number of bases between two adjacent microsatellites in
##          order to specify the compound microsatellite type.
##    Example:
##      definition(unit_size,min_repeats):          1-10 2-6 3-5 4-5 5-5 6-5
##      interruptions(max_difference_for_2_SSRs):   100
##
## EXAMPLE: misa.pl seqs.fasta
## Modified by Leshi Chen for primer design
## _______________________________________________________________________________
##


#§§§§§ DECLARATION §§§§§#

# Check for arguments. If none display syntax #


if (@ARGV == 0)
  {
  open (IN,"<$0");
  while (<IN>) {if (/^\#\# (.*)/) {$message .= "$1\n"}};
  close (IN);
  die $message;
  };

# Check if help is required #

if ($ARGV[0] =~ /-help/i)
  {
  open (IN,"<$0");
  while (<IN>) {if (/^\#\#\#(.*)/) {$message .= "$1\n"}};
  close (IN);
  die $message;
  };

# Open FASTA file #

open (IN,"<$ARGV[0]") || die ("\nError: FASTA file doesn't exist !\n\n");
#open (OUT,">$ARGV[0].misa"); updated by Leshi chen for galaxy integration
open (OUT,">$ARGV[1]");
print OUT "ID\tSSR nr.\tSSR type\tSSR\tsize\tstart\tend\n";

# Reading arguments updated by Leshi chen to get local path otherwise will create error #
#use Cwd 'abs_path';
#use Cwd 'getcwd';
#print getcwd()&"misa.ini";
#print OUT abs_path($0);
#open (SPECS,"\/root\/galaxy_dist\/tools\/pfr_2010\/"."misa.ini") || die ("\nError: Specifications file doesn't exist ! \n\n misa.ini not found ! \n\n");
my $arg_def= $ARGV[2]||'';
my $arg_interuption= $ARGV[3]||'';
#my $tmb = '';
#my $_ = '';
my %typrep;
my $amb = 0;

%typrep =  $arg_def  =~/(\d+)-(\d+)/gi;
#print "1:" , $arg_def , "\n";
#print "hh: ",  %typrep , "\n";
#print $arg_def , "\n";
#print $arg_interuption ,"\n";
#print  $arg_def  =~/(\d+)/gi , "\n";
#%typrep =  $arg_def  =~/(\d+)/gi;
print %typrep , "\n";
$amb = $arg_interuption;
print $amb , "\n";
#while (<SPECS>)#
  # {#
  # %typrep = $1 =~ /(\d+)/gi if (/^def\S*\s+(.*)/i);#
  # if (/^int\S*\s+(\d+)/i) {$amb = $1}#
  # };#
my @typ = sort { $a <=> $b } keys %typrep;
print @typ . "\n";
#die (%typrep , "--" , @typ , "--" , $amb);
#§§§§§ CORE §§§§§#

$/ = ">";
my $max_repeats = 1; #count repeats
my $min_repeats = 1000; #count repeats
my (%count_motif,%count_class); #count
my ($number_sequences,$size_sequences,%ssr_containing_seqs); #stores number and size of all sequences examined
my $ssr_in_compound = 0;
my ($id,$seq);
while (<IN>)
  {
  next unless (($id,$seq) = /(.*?)\n(.*)/s);
  my ($nr,%start,@order,%end,%motif,%repeats); # store info of all SSRs from each sequence
  $seq =~ s/[\d\s>]//g; #remove digits, spaces, line breaks,...
  $id =~ s/^\s*//g; $id =~ s/\s*$//g;$id =~ s/\s/_/g; #replace whitespace with "_"
  $number_sequences++;
  $size_sequences += length $seq;
  for ($i=0; $i < scalar(@typ); $i++) #check each motif class
    {
    my $motiflen = $typ[$i];
    my $minreps = $typrep{$typ[$i]} - 1;
    if ($min_repeats > $typrep{$typ[$i]}) {$min_repeats = $typrep{$typ[$i]}}; #count repeats
    my $search = "(([acgt]{$motiflen})\\2{$minreps,})";
    while ( $seq =~ /$search/ig ) #scan whole sequence for that class
      {
      my $motif = uc $2;
      my $redundant; #reject false type motifs [e.g. (TT)6 or (ACAC)5]
      for ($j = $motiflen - 1; $j > 0; $j--)
        {
        my $redmotif = "([ACGT]{$j})\\1{".($motiflen/$j-1)."}";
        $redundant = 1 if ( $motif =~ /$redmotif/ )
        };
      next if $redundant;
      $motif{++$nr} = $motif;
      my $ssr = uc $1;
      $repeats{$nr} = length($ssr) / $motiflen;
      $end{$nr} = pos($seq);
      $start{$nr} = $end{$nr} - length($ssr) + 1;
      # count repeats
      # count_motifs doesn't required as statistic has been removed - modified by leshi
      #$count_motifs{$motif{$nr}}++; #counts occurrence of individual motifs
      $motif{$nr}->{$repeats{$nr}}++; #counts occurrence of specific SSR in its appearing repeat
      $count_class{$typ[$i]}++; #counts occurrence in each motif class
      if ($max_repeats < $repeats{$nr}) {$max_repeats = $repeats{$nr}};
      };
    };
  next if (!$nr); #no SSRs
  $ssr_containing_seqs{$nr}++;
  @order = sort { $start{$a} <=> $start{$b} } keys %start; #put SSRs in right order
  $i = 0;
  my $count_seq; #counts
  my ($start,$end,$ssrseq,$ssrtype,$size);
  while ($i < $nr)
    {
    my $space = $amb + 1;
    if (!$order[$i+1]) #last or only SSR
      {
      $count_seq++;
      my $motiflen = length ($motif{$order[$i]});
      $ssrtype = "p".$motiflen;
      $ssrseq = "($motif{$order[$i]})$repeats{$order[$i]}";
      $start = $start{$order[$i]}; $end = $end{$order[$i++]};
      next
      };
    if (($start{$order[$i+1]} - $end{$order[$i]}) > $space)
      {
      $count_seq++;
      my $motiflen = length ($motif{$order[$i]});
      $ssrtype = "p".$motiflen;
      $ssrseq = "($motif{$order[$i]})$repeats{$order[$i]}";
      $start = $start{$order[$i]}; $end = $end{$order[$i++]};
      next
      };
    my ($interssr);
    if (($start{$order[$i+1]} - $end{$order[$i]}) < 1)
      {
      $count_seq++; $ssr_in_compound++;
      $ssrtype = 'c*';
      $ssrseq = "($motif{$order[$i]})$repeats{$order[$i]}($motif{$order[$i+1]})$repeats{$order[$i+1]}*";
      $start = $start{$order[$i]}; $end = $end{$order[$i+1]}
      }
    else
      {
      $count_seq++; $ssr_in_compound++;
      $interssr = lc substr($seq,$end{$order[$i]},($start{$order[$i+1]} - $end{$order[$i]}) - 1);
      $ssrtype = 'c';
      $ssrseq = "($motif{$order[$i]})$repeats{$order[$i]}$interssr($motif{$order[$i+1]})$repeats{$order[$i+1]}";
      $start = $start{$order[$i]};  $end = $end{$order[$i+1]};
      #$space -= length $interssr
      };
    while ($order[++$i + 1] and (($start{$order[$i+1]} - $end{$order[$i]}) <= $space))
      {
      if (($start{$order[$i+1]} - $end{$order[$i]}) < 1)
        {
        $ssr_in_compound++;
        $ssrseq .= "($motif{$order[$i+1]})$repeats{$order[$i+1]}*";
        $ssrtype = 'c*';
        $end = $end{$order[$i+1]}
        }
      else
        {
        $ssr_in_compound++;
        $interssr = lc substr($seq,$end{$order[$i]},($start{$order[$i+1]} - $end{$order[$i]}) - 1);
        $ssrseq .= "$interssr($motif{$order[$i+1]})$repeats{$order[$i+1]}";
        $end = $end{$order[$i+1]};
        #$space -= length $interssr
        }
      };
    $i++;
    }
  continue
    {
    print OUT "$id\t$count_seq\t$ssrtype\t$ssrseq\t",($end - $start + 1),"\t$start\t$end\n"
    };
  };

close (OUT);
#open (OUT,">$ARGV[0].statistics"); updated by Leshi chen for galaxy integration
# the statistics part has been removed as we only need misa for primer

