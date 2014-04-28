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

my $lib_data_file = '';
my $quality_counts_file = '';

GetOptions ('lib_data_file=s'=>\$lib_data_file, 'quality_counts_file=s'=>\$quality_counts_file);

print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify a library data file using:  --lib_data_file", RESET;
print GREEN, "\n\tSpecify a quality counts file using: --quality_counts_file", RESET;
print GREEN, "\n\tCounts from this file will be merged into the library data file", RESET;
print GREEN, "\n\nExample:  mergeQualityReadCounts.pl  --lib_data_file=/projects/malachig/alexa_seq/batch_jobs/LnCAP_AR_KnockIn/LnCAP_AR_KnockIn_Lib_Data.txt  --quality_counts_file=/projects/malachig/alexa_seq/batch_jobs/LnCAP_AR_KnockIn/HS1441/qualityReadCounts.txt\n\n", RESET;

unless ($lib_data_file && $quality_counts_file){
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  exit();
}

unless (-e $lib_data_file && -e $quality_counts_file){
  print RED, "\n\nOne of the input files was not found\n\n", RESET;
  exit();
}

#Get values from the quality counts file
my %qualities;
open(QUAL, "$quality_counts_file") || die "\n\nCould not open quality counts file: $quality_counts_file\n\n";
while(<QUAL>){
  chomp($_);
  my @line = split("\t", $_);
  if ($line[0] =~ /(.*)\.txt\.gz/){
    my $id = $1;
    my $count = $line[1];
    $qualities{$id}{count} = $count;
  }
}
close(QUAL);
#print Dumper %qualities;

#Get entries from the library data file
my %lanes;
open(DATA, "$lib_data_file") || die "\n\nCould not open lib data file: $lib_data_file";
my $order = 0;
while(<DATA>){
  $order++;
  chomp($_);
  my @line = split("\t", $_);
  $lanes{$line[2]}{line} = \@line;
  $lanes{$line[2]}{order} = $order;
}
close(DATA);
#print Dumper %lanes;

#Replace NA value with actual counts in the lib data entries
foreach my $flowcell_lane (sort {$lanes{$a}->{order} <=> $lanes{$b}->{order}} keys %lanes){
  if ($qualities{$flowcell_lane}){
    my $count = $qualities{$flowcell_lane}{count};
    my @line = @{$lanes{$flowcell_lane}{line}};
    $line[7] = $count;
    $lanes{$flowcell_lane}{line} = \@line;
  }
}

#Print out the updated lane data records to a temp file
my $temp_file = "$lib_data_file".".tmp";
open (TEMP, ">$temp_file") || die "\n\nCould not open temp output file: $temp_file\n\n";
foreach my $flowcell_lane (sort {$lanes{$a}->{order} <=> $lanes{$b}->{order}} keys %lanes){
  my @line = @{$lanes{$flowcell_lane}{line}};
  print TEMP "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$line[6]\t$line[7]\t$line[8]\n";

}
close(TEMP);

#Overwrite the original lane data file with the temp one
my $cmd = "mv -f $temp_file $lib_data_file";
print BLUE, "\n\n$cmd\n\n", RESET;
system($cmd);

exit();

