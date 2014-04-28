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

GetOptions ('read_record_dir=s'=>\$read_record_dir, 'result_file=s'=>\$result_file);

print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify a root read record directory containing files for each flowcell-lane using:  --read_record_dir", RESET;
print GREEN, "\n\nExample:  getQualityReadCounts.pl  --read_record_dir=/projects/malachig/solexa/read_records/HS0499/  --result_file=/projects/malachig/alexa_seq/batch_jobs/Neuroblastoma/Neuroblastoma_Lib_Data.txt\n\n", RESET;

unless ($read_record_dir && $result_file){
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  exit();
}

$read_record_dir = &checkDir('-dir'=>$read_record_dir, '-clear'=>"no");

#Open the lib data file
my %lib_data;
open(RESULT, ">$result_file") || die "\n\nCould not open result file: $result_file (getQualityReadCounts.pl)\n\n";

opendir(DIRHANDLE, "$read_record_dir") || die "\nCannot open directory: $read_record_dir\n\n";
my @files = readdir(DIRHANDLE);
closedir(DIRHANDLE);

foreach my $file (@files){

  #Skip sub-directories
  if (-d $file){
    next();
  }
  #Skip non-compressed non-text files
  unless ($file =~ /\.txt\.gz/){
    next();
  }

  print BLUE "\n$file", RESET;
  print RESULT "$file";
  open(IN, "zcat $read_record_dir$file | ") || die "\nCould not open input file\n\n";
  my $header = 1;
  my $passread = 0;
  while(<IN>){
    if ($header == 1){
      $header = 0;
      next();
    }
    if ($_ =~ /\S+\s+\S+\s+\S+\s+(\S+)\s+(\S+)/){
      my $r1 = $1; 
      my $r2 = $2;
      unless ($r1 =~ /Low_Quality|Low_Complexity|Duplicate|Read1_Status/){$passread++}; 
      unless ($r2 =~ /Low_Quality|Low_Complexity|Duplicate|Read2_Status/){$passread++};
    }
  }
  close(IN);
  print BLUE "\t$passread", RESET;
  print RESULT "\t$passread\n";
}

print "\n\n";
close(RESULT);
exit();

