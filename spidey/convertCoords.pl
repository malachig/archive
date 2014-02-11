#!/usr/bin/perl


#
#This program modifies the coordinates to be relative to the start/end of a chromosome
#
#By: Brian Cho, BCGSC
#Supervisor: Malachi Griffith
#

use strict;
use Getopt::Long;
use Term::ANSIColor qw(:constants);

my $strand = '';
my $start = '';
my $end = '';

GetOptions( 'strand=s' => \$strand,
            'chr_start=i' => \$start,
            'chr_end=i' => \$end );

unless ($strand && $start && $end)
{
  print MAGENTA, "Please specify the strandedness, start coordinates, and end coordinates (of the gene!), using the following format:\n";
  print "--strand='+' --chr_start=14606372 --chr_end=14694051\n";
  print "convertCoords.pl terminated\n\n", RESET;
  exit();
}

while (my $line = <STDIN>)
{
  chomp $line;
  if ($line =~ /^(\d+)\:(\d+)\t\d+\t\d+\t(\d+)\t(\d+)\t(\d+)/)
  {
    #Convert by adding to start coord of the gene relative to chromosome
    #print $1+$start,"\:",$2+$start,"\t",$1+$start,"\t",$2+$start,"\t",$3,"\t",$4,"\t",$5,"\n";

    #Same as above, except subtract by 1 to fix off-by-1 error
    print $1+$start-1,"\:",$2+$start-1,"\t",$1+$start-1,"\t",$2+$start-1,"\t",$3,"\t",$4,"\t",$5,"\n";
    
    #Convert by subtracting from end coord of the gene relative to chromosome
    #print $end-$1,"\:",$end-$2,"\t",$end-$1,"\t",$end-$2,"\t",$3,"\t",$4,"\t",$5,"\n";

    #Same as above, except add/subtract by 1 to fix the off-by-1 error
  }
  

}

exit();
