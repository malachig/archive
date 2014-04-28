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

my $read_record_dir = '';
my $result_file = '';
my $library_id = '';

GetOptions ('read_record_dir=s'=>\$read_record_dir, 'result_file=s'=>\$result_file, 'library_id=s'=>\$library_id);

print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the library Id for this library using:  --library_id", RESET;
print GREEN, "\n\tSpecify a root read record directory containing files for each flowcell-lane using:  --read_record_dir", RESET;
print GREEN, "\n\tSpecify an output file to store the result using: --result_file", RESET;
print GREEN, "\n\tThe average read length for the entire library will be determined", RESET;
print GREEN, "\n\nExample:  getAverageReadLength.pl  --library_id=HS1441  --read_record_dir=/projects/malachig/alexa_seq/read_records/HS1441/  --result_file=/projects/malachig/alexa_seq/batch_jobs/LnCAP_AR_KnockIn/HS1441/averageReadLength.txt\n\n", RESET;

unless ($library_id && $read_record_dir && $result_file){
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  exit();
}

$read_record_dir = &checkDir('-dir'=>$read_record_dir, '-clear'=>"no");

opendir(DIRHANDLE, "$read_record_dir") || die "\nCannot open directory: $read_record_dir\n\n";
my @files = readdir(DIRHANDLE);
closedir(DIRHANDLE);

my $cumulative_length = 0;
my $read_count = 0;
my %lanes;
my $fc = 0;
foreach my $file (@files){
  #Skip sub-directories
  if (-d $file){
    next();
  }
  #Skip non-compressed non-text files
  my $name;
  if ($file =~ /(.*)\.txt\.gz/){
    $name=$1;
  }else{
    next();
  }
  $fc++;
  $lanes{$fc}{name} = $name;
  $lanes{$fc}{read_count}=0;
  $lanes{$fc}{cumulative_length}=0;
  print BLUE "\n$file", RESET;
  open(IN, "zcat $read_record_dir$file | ") || die "\nCould not open input file\n\n";
  my $header = 1;
  while(<IN>){
    if ($header == 1){
      $header = 0;
      next();
    }
    chomp($_);
    my @line = split("\t", $_);
    my $read1 = $line[5];
    my $read2 = $line[6];
    $cumulative_length += length($read1);
    $lanes{$fc}{cumulative_length} += length($read1);
    $read_count++;
    $lanes{$fc}{read_count}++;
    $cumulative_length += length($read2);
    $lanes{$fc}{cumulative_length} += length($read2);
    $read_count++;
    $lanes{$fc}{read_count}++;
  }
  close(IN);
}

#Calculate averagess
my $average_read_length = ($cumulative_length/$read_count);
foreach my $fc (keys %lanes){
  $lanes{$fc}{average_read_length} = ($lanes{$fc}{cumulative_length}/$lanes{$fc}{read_count});
}

#Store results in output file
open(RESULT, ">$result_file") || die "\n\nCould not open result file: $result_file (getQualityReadCounts.pl)\n\n";
foreach my $fc (sort {$a <=> $b} keys %lanes){
  print RESULT "$lanes{$fc}{name}\t$lanes{$fc}{average_read_length}\n";
}
print RESULT "$library_id\t$average_read_length\n";
close(RESULT);
print "\n\n";

exit();

