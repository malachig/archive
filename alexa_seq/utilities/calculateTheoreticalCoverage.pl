#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);

my $er_size = 62;
my $j_size = 62;
my $read_size = 42;
my $transcript_size = 200;

#A.) Imagine a theoretical 'Exon Region' of 62-bases within a theoretical transcript of 200 bases (say from position 69 to 131)
print BLUE, "\nSimulate an exon region within a transcript", RESET;
my %coverage;

#Initialize the transcript
for (my $i = 1; $i <= $transcript_size; $i++){
  $coverage{$i} = 0;
}

#Throw every possible read of the specified length onto this transcript
my $end = 0;
my $start = 0;
while($end <= $transcript_size){
  $start++;
  $end = ($start+$read_size)-1;

  for (my $i = $start; $i <= $end; $i++){
    $coverage{$i}++;
  }
}

#Print out the coverage object
foreach my $position (sort {$a <=> $b} keys %coverage){
  print YELLOW "\n$position = $coverage{$position}", RESET;
}

#What is the cumulative coverage of a theoretical exon region with the coordinates 69-131?
my $cumulative_coverage = 0;
for (my $i = 69; $i <= 131-1; $i++){
  $cumulative_coverage += $coverage{$i};
}

print BLUE, "\n\nCoverage potential of an exon region centred within this theoretical transcript is: $cumulative_coverage", RESET;

#B.) Now imagine a theoretical exon junction sequence of size 62
print BLUE, "\n\n\nSimulate an exon junction within a transcript", RESET;

#Initialize the junction
my %coverage_j;
for (my $i = 1; $i <= $j_size; $i++){
  $coverage_j{$i} = 0;
}

my $end = 0;
my $start = 0;
while($end <= $j_size){
  $start++;
  $end = ($start+$read_size)-1;
  for (my $i = $start; $i <= $end; $i++){
    $coverage_j{$i}++;
  }
}

#Print out the coverage object
foreach my $position (sort {$a <=> $b} keys %coverage_j){
  print YELLOW "\n$position = $coverage_j{$position}", RESET;
}

#What is the cumulative coverage of this theoretical junction?
my $cumulative_coverage_j = 0;
for (my $i = 1; $i <= $j_size; $i++){
  $cumulative_coverage_j += $coverage_j{$i};
}

print BLUE, "\n\nCoverage potential of an exon junction of $j_size bases is: $cumulative_coverage_j" RESET;






