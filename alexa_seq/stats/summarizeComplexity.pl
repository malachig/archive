#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#For a particular library, grab the estimate of library complexity previously calculated
#Grab the library complexity for all other libraries available
#Create a box-plot to display the library complexity of these libraries and mark the current library on this distribution

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);

#Load the ALEXA libraries
BEGIN {
  use Cwd 'abs_path';
  if (abs_path($0) =~ /(.*)\/.*\/.*\.pl/){
    push (@INC, $1);
  }
}
use utilities::utility qw(:all);

my $analysis_dir = '';
my $library_id = '';
my $script_dir = '';
my $complexity_type = '';

GetOptions ('analysis_dir=s'=>\$analysis_dir, 'library_id=s'=>\$library_id, 'script_dir=s'=>\$script_dir, 'complexity_type=s'=>\$complexity_type);

print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the main analysis directory using: --analysis_dir", RESET;
print GREEN, "\n\tSpecify the library to summarize using: --library_id", RESET;
print GREEN, "\n\tSpecify the main script dir using: --script_dir", RESET;
print GREEN, "\n\tSpecify whether to evaluate complexity by the sequences (SEQ) or their mapping positions (MAP) using: --complexity_type", RESET;
print GREEN, "\n\nExample: summarizeComplexity.pl  --analysis_dir=/projects/malachig/alexa_seq/  --library_id=HS04391  --script_dir=/home/malachig/svn/alexa_seq/  --complexity_type=SEQ\n\n", RESET;

unless ($analysis_dir && $library_id && $script_dir && $complexity_type){
  print RED, "\nRequired input parameter(s) missing!\n\n", RESET;
  exit();
}
unless ($complexity_type =~ /SEQ|MAP/i){
  print RED, "\n--complexity type must be 'SEQ or 'MAP'\n\n", RESET;
  exit();
}

#Check input directories
$analysis_dir = &checkDir('-dir'=>$analysis_dir, '-clear'=>"no");
$script_dir = &checkDir('-dir'=>$script_dir, '-clear'=>"no");
my $r_script_dir = "$script_dir"."R_bin/";
$r_script_dir = &checkDir('-dir'=>$r_script_dir, '-clear'=>"no");
my $r_script = "$r_script_dir"."summarizeComplexity.R";
my $stats_dir = "$analysis_dir"."figures_and_stats/";
$stats_dir = &checkDir('-dir'=>$stats_dir, '-clear'=>"no");
my $outdir = "$stats_dir"."$library_id/LibraryQuality/";
$outdir = &checkDir('-dir'=>$outdir, '-clear'=>"no");
my $temp_file = "$outdir"."temp.txt" ;

#Get all the libraryComplexity file
my @file_list;
my $file_list;
if ($complexity_type =~ /SEQ/i){
  $file_list = "ls $analysis_dir"."figures_and_stats/*/LibraryQuality/LibraryComplexity_SEQ.txt";
}elsif($complexity_type =~ /MAP/i){
  $file_list = "ls $analysis_dir"."figures_and_stats/*/LibraryQuality/LibraryComplexity_MAP.txt";
}
my @temp = `$file_list`;

push(@file_list, @temp);

my %files;
my $order = 0;
foreach my $path (@file_list){
  chomp($path);
  if ($path =~ /figures_and_stats\/(.*)\/LibraryQuality/){
    $order++;
    $files{$1}{path}=$path;
    $files{$1}{order}=$order;
    if ($1 eq "$library_id"){
      $files{$1}{order}=0;
    }
  }
}

my $lib_count = keys %files;
if ($lib_count <= 5){
  print RED, "\n\nFound 5 or fewer libraries processed so far - rerun this job later when you have processed more libraries", RESET;
  exit();
}else{
  print BLUE, "\n\nFound $lib_count library complexity values to summarize", RESET;
}

#Build a temp file containing all the values for each library.  Place the current library on the first line of this file
open(TEMP, ">$temp_file") || die "\n\nCould not open temp output file: $temp_file\n\n";
print TEMP "LibraryID\tOrder\tUniqueReadsPerMillion\tPercentUniqueReadsPerMillion\tRedundantReadsPerMillion\tPercentRedundantReadsPerMillion\tDistinctRedundantReads\tPercentDistinctRedundantReads\n";
foreach my $library (sort {$files{$a}->{order} <=> $files{$b}->{order}} keys %files){
  my $order = $files{$library}{order};
  my $path = $files{$library}{path};
  open(IN, "$path") || die "\n\nCould not open input file: $path\n\n";
  my $header = 1;
  while(<IN>){
    if ($header == 1){
      $header = 0;
      next();
    }
    chomp($_);
    my @line = split("\t", $_);

    if ($line[2] =~ /(.*)\((.*)\%\)/){
      print TEMP "$library\t$order\t$1\t$2";
    }else{
      print TEMP "$library\t$order\tNA\tNA";
    }

    if ($line[3] =~ /(.*)\((.*)\%\)/){
      print TEMP "\t$1\t$2";
    }else{
      print TEMP "\tNA\tNA";
    }

    if ($line[4] =~  /(.*)\((.*)\%\)/){
      print TEMP "\t$1\t$2\n";
    }else{
      print TEMP "\tNA\tNA\n";
    }
  }
  close(IN);
}
close(TEMP);

#Execute the R script
my $type_uc = uc($complexity_type);
my $file_base = "LibraryComplexity_"."$type_uc";
my $r_cmd = "$r_script $outdir $file_base";
print BLUE, "\n\n$r_cmd", RESET;
system($r_cmd);

#Delete the temp file
my $rm_cmd = "rm -f $temp_file";
print BLUE, "\n\n$rm_cmd", RESET;
system($rm_cmd);

print "\n\n";

exit();

