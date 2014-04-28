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
use File::Basename;

#Load the ALEXA modules
BEGIN {
  use Cwd 'abs_path';
  if (abs_path($0) =~ /(.*)\/.*\/.*\.pl/){
    push (@INC, $1);
  }
}
use utilities::utility qw(:all);

my $lib_names_file = '';
my $average_length_file = '';

GetOptions ('lib_names_file=s'=>\$lib_names_file, 'average_length_file=s'=>\$average_length_file);

print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify a library names file using:  --lib_names_file", RESET;
print GREEN, "\n\tSpecify an average read length file using: --average_length_file", RESET;
print GREEN, "\n\tCounts from this file will be merged into the library names file", RESET;
print GREEN, "\n\nExample:  mergeQualityReadCounts.pl  --lib_names_file=/projects/malachig/alexa_seq/batch_jobs/LnCAP_AR_KnockIn/LnCAP_AR_KnockIn_Lib_Names.txt  --average_length_file=/projects/malachig/alexa_seq/batch_jobs/LnCAP_AR_KnockIn/HS1441/averageReadLength.txt\n\n", RESET;

unless ($lib_names_file && $average_length_file){
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  exit();
}
unless (-e $lib_names_file && -e $average_length_file){
  print RED, "\n\nOne of the input files was not found\n\n", RESET;
  exit();
}

#Get values from the average lengths file
my %lengths;
open(LENGTH, "$average_length_file") || die "\n\nCould not open average_length file: $average_length_file\n\n";
while(<LENGTH>){
  chomp($_);
  my @line = split("\t", $_);
  $lengths{$line[0]}{average_length} = $line[1];
}
close(LENGTH);
#print Dumper %lengths;

#Get entries from the library data file
my %libs;
open(LIB, "$lib_names_file") || die "\n\nCould not open lib names file: $lib_names_file";
my $order = 0;
while(<LIB>){
  $order++;
  chomp($_);
  my @line = split("\t", $_);
  $libs{$line[0]}{line} = \@line;
  $libs{$line[0]}{order} = $order;
}
close(LIB);
#print Dumper %libs;

#Replace NA value with actual counts in the lib data entries
foreach my $lib (sort keys %lengths){
  if ($libs{$lib}){
    my $length = $lengths{$lib}{average_length};
    my @line = @{$libs{$lib}{line}};
    $line[4] = $length;
    $libs{$lib}{line} = \@line;
  }
}

#Print out the updated lane data records to a temp file
my $temp_file = "$lib_names_file".".tmp";
open (TEMP, ">$temp_file") || die "\n\nCould not open temp output file: $temp_file\n\n";
foreach my $lib (sort {$libs{$a}->{order} <=> $libs{$b}->{order}} keys %libs){
  my @line = @{$libs{$lib}{line}};
  print TEMP "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\n";
  #print YELLOW, "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\n", RESET;
}
close(TEMP);

#Overwrite the original lane data file with the temp one
my $cmd = "mv -f $temp_file $lib_names_file";
print BLUE, "\n\n$cmd\n\n", RESET;
system($cmd);

exit();

