#!/usr/bin/perl


#
#This program converts the input from junction.pl into a GFF file format
#
#By: Brian Cho, BCGSC
#Supervisor: Malachi Griffith
#

use strict;
use Getopt::Long;
use Term::ANSIColor qw(:constants);

my $start;
my $end;
my $count;
my $donor;
my $acceptor;
my $groupID;
my $startless20;
my $endmore20;
my @count_array;
my $n = 0;

my $min_count = '';
my $track_name = '';
my $track_color = '';
my $strand = '';
my $chr = '';

GetOptions( 'strand=s' => \$strand,
            'chr=i' => \$chr,
            'min_count=i' => \$min_count,
            'track_name=s' => \$track_name,
            'track_color=s' => \$track_color );

unless ($strand && $chr && $min_count && $track_name && $track_color)
{
  print MAGENTA, "Please specify the strandedness (of the gene!), chromosome number, minimum count, library ID, and track_color, using the following format:\n";
  print "--strand='-' --chr=10 --min_count=10 --track_name=RA0001 --track_color=255,0,0\n";
  print "makeGFF.pl terminated\n\n", RESET;
  exit();
}

my %lines;
my $lc = 0;
while (my $line = <STDIN>)
{
  
  chomp $line;
  
  if ($line =~ /^(\d+)\:(\d+)\t\d+\t\d+\t(\d+)\t(\d+)\t(\d+)/)
  {
    $start = $1;
    $end = $2;
    $count = $3;
    $donor = $4;
    $acceptor = $5;
  }
  
  #Filtering only for canonical donor and canonical acceptor; as well as filtering for reads with counts greater than or equal to $min_count
  if ($donor == 1 && $acceptor == 1 && $count >= $min_count)
  {
    $lc++;
    $count_array[$lc-1] = $count;

    $startless20 = $start - 20;
    $endmore20 = $end + 20;

    $lines{$lc}{record1}="chr$chr\t\tjunction\t$startless20\t$start\t.\t$strand\t.";
    $lines{$lc}{record2}="chr$chr\t\tjunction\t$end\t$endmore20\t.\t$strand\t.";

    $lines{$lc}{start} = $startless20;
    $lines{$lc}{end} = $endmore20;
  }

  

}

$groupID = $lc;

my $grand_start = $lines{1}{start};
my $grand_end = $lines{$lc}{end};

print "browser position chr$chr:$grand_start-$grand_end\n";
print "browser hide all\n";
print "browser full knownGene\n";
print "track name=$track_name description=\"Junctions of $track_name\" visibility=2 color=$track_color\n";

foreach my $lc (sort {$a <=> $b} keys %lines)
{
  print $lines{$lc}{record1}; print "\t$groupID\_$count_array[$n]\n";
  print $lines{$lc}{record2}; print "\t$groupID\_$count_array[$n]\n";
  $groupID--; $n++;
}

exit();

