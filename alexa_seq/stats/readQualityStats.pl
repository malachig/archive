#!/usr/bin/perl
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#DESCRIPTION: takes in  a *.fastq file STREAM and summarizes the qualities observed
#This fastq file must be the type generated from PRB and SEQ files by MAQ!!

#Build a cumulative distribution of all qualities and do the same thing by read position

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);

my %quality_counts;
my %pos_qualities;

my $total = 0;
my $counter = 0;
my $progress_counter = 0;
my $grand_quality_count = 0;

while (<>){
  $counter++;

  #Qualities are on every 4th line of a fastq file
  if ($counter == 4){
    $counter = 0;
    chomp($_);

    $progress_counter++;
    if ($progress_counter == 100000){
      $| = 1; print STDERR "."; $| = 0;
      $progress_counter = 0;
    }

    my $qual_string = $_;
    my @qual_asci = split (//, $qual_string);

    my $pos = 0;

    foreach my $qual (@qual_asci){
      $pos++;

      my $quality = ord($qual)-33;
      $grand_quality_count++;

      if ($quality_counts{$quality}){
	$quality_counts{$quality}{qual}++;
      }else{
	$quality_counts{$quality}{qual} = 1;
      }

      if ($pos_qualities{$pos}){
	$pos_qualities{$pos}{cumulative_quality} += $quality;
	$pos_qualities{$pos}{quality_count}++;
      }else{
	$pos_qualities{$pos}{cumulative_quality} = $quality;
	$pos_qualities{$pos}{quality_count} = 1;
      }
    }
  }
}

print STDERR "\n";

print "#Total Quality counts:\n";
print "#Quality\tCount\tPercent Of Total\tCumulative Percent (Forward)\tCumulative Percent (Reverse)\n";

my $qual_percent_cum_forward = 0;
foreach my $quality (sort {$a <=> $b} keys %quality_counts){
  my $quality_percent = ($quality_counts{$quality}{qual}/$grand_quality_count)*100;
  $qual_percent_cum_forward += $quality_percent;
  $quality_counts{$quality}{qual_percent} = $quality_percent;
  $quality_counts{$quality}{qual_percent_cum_forward} = $qual_percent_cum_forward;
}

my $qual_percent_cum_reverse = 0;
foreach my $quality (sort {$b <=> $a} keys %quality_counts){
  my $quality_percent = ($quality_counts{$quality}{qual}/$grand_quality_count)*100;
  $qual_percent_cum_reverse += $quality_percent;
  $quality_counts{$quality}{qual_percent_cum_reverse} = $qual_percent_cum_reverse;
}

foreach my $quality (sort {$a <=> $b} keys %quality_counts){
  print "$quality\t$quality_counts{$quality}{qual}\t$quality_counts{$quality}{qual_percent}\t$quality_counts{$quality}{qual_percent_cum_forward}\t$quality_counts{$quality}{qual_percent_cum_reverse}\n";
}

print "#\n#Average Quality at each read position:\n";
print "#POS\tQuality Count\tCumulative Quality\tAverage Quality\n";
foreach my $pos (sort {$a <=> $b} keys %pos_qualities){

  my $average_pos_quality = $pos_qualities{$pos}{cumulative_quality} / $pos_qualities{$pos}{quality_count};

  print "$pos\t$pos_qualities{$pos}{quality_count}\t$pos_qualities{$pos}{cumulative_quality}\t$average_pos_quality\n";
}

print "\n\n";

exit();
