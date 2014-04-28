#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#The purpose of this script is to check for completion of BLAST parsing jobs for a particular sequence database
#Get a list of all libraries and lanes for a particular project from an input file
#For each lane of each library, check the corresponding LOG file for completion status
#Also check for uncompressed results files that might indicate a failed job

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

my $project_name = '';
my $analysis_dir = '';
my $ensembl_version = '';

GetOptions ('project_name=s'=>\$project_name, 'analysis_dir=s'=>\$analysis_dir, 'ensembl_version=i'=>\$ensembl_version);

print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the main analysis directory using: --analysis_dir", RESET;
print GREEN, "\n\tSpecify the project name for this project using:  --project_name", RESET;
print GREEN, "\n\tSpecify the ensembl version used for this analysis using: --ensembl_version", RESET;
print GREEN, "\n\nExample:  checkParsingJobs.pl  --project_name=FL_Trans  --analysis_dir=/projects/malachig/alexa_seq/  --ensembl_version=53\n\n", RESET;

unless ($project_name && $analysis_dir && $ensembl_version){
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  exit();
}
$analysis_dir = &checkDir('-dir'=>$analysis_dir, '-clear'=>"no");

#Get the libraries and lanes from a the lib_data file
my $lib_lanes_file = "$analysis_dir"."batch_jobs/$project_name/$project_name"."_Lib_Data.txt";
print BLUE, "\n\nGetting list of libraries and lanes for this project from: $lib_lanes_file", RESET;
open(LANES, "$lib_lanes_file") || die "\n\nCould not open lib data file: $lib_lanes_file\n\n";
my %lanes;
my $c = 0;
while(<LANES>){
  $c++;
  chomp($_);
  my @line = split("\t",$_);
  $lanes{$c}{lib} = $line[0];
  $lanes{$c}{lane} = $line[2];
}
close(LANES);

my $total_lane_count = keys %lanes;

#Define parsing types
my %types;
$types{'1'}{name} = "Repeat";
$types{'1'}{dir} = "Repeats";
$types{'1'}{file} = "Repeats";
$types{'2'}{name} = "GeneExon";
$types{'2'}{dir} = "ENST_v"."$ensembl_version";
$types{'2'}{file} = "ENST_v"."$ensembl_version";
$types{'3'}{name} = "Junction";
$types{'3'}{dir} = "Junctions_v"."$ensembl_version";
$types{'3'}{file} = "Junctions_v"."$ensembl_version";
$types{'4'}{name} = "Boundary";
$types{'4'}{dir} = "Boundaries_v"."$ensembl_version";
$types{'4'}{file} = "Boundaries_v"."$ensembl_version";
$types{'5'}{name} = "Intron";
$types{'5'}{dir} = "Introns_v"."$ensembl_version";
$types{'5'}{file} = "Introns_v"."$ensembl_version";
$types{'6'}{name} = "Intergenic";
$types{'6'}{dir} = "Intergenics_v"."$ensembl_version";
$types{'6'}{file} = "Intergenics_v"."$ensembl_version";

#First check the read records file for each lib/lane to see if there are any uncompressed files
#If present, these indicate that jobs either in progress or have failed at some point...
print BLUE, "\n\nChecking for presence of read record files:", RESET;
my $rr_count_gz = 0;
my $rr_count_txt = 0;
foreach my $c (sort {$a <=> $b} keys %lanes){
  my $lib = $lanes{$c}{lib};
  my $lane = $lanes{$c}{lane};
  my $rr_file_gz = "$analysis_dir"."read_records/$lib"."/$lane".".txt.gz";
  my $rr_file_txt = "$analysis_dir"."read_records/$lib"."/$lane".".txt";

  if (-e $rr_file_gz){
    $rr_count_gz++;
  }
  if (-e $rr_file_txt){
    $rr_count_txt++;
  }
}
if ($rr_count_txt){
  print YELLOW, "\n\tFound $rr_count_gz Read Record files and $rr_count_txt uncompressed Read Record files (in process or represent errors?)", RESET;
}elsif($rr_count_gz != $total_lane_count){
  print YELLOW, "\n\tFound $rr_count_gz Read Record files BUT expected $total_lane_count", RESET;
}else{
  print BLUE, "\n\tFound $rr_count_gz Read Record files (none in process)", RESET;
}

#Now go through each sequence type and check for uncompressed mapping result files
#Also check for completion status in the log files for each job
foreach my $t (sort {$a <=> $b} keys %types){
  my $type_name = $types{$t}{name};
  my $type_dir = $types{$t}{dir};
  my $type_file = $types{$t}{file};
  print MAGENTA, "\n\nType: $type_name", RESET;

  my $s_count_gz = 0;
  my $s_count_txt = 0;
  my $log_count = 0;
  my $completed_job_count = 0;

  foreach my $c (sort {$a <=> $b} keys %lanes){
    my $lib = $lanes{$c}{lib};
    my $lane = $lanes{$c}{lane};

    #First check for expected mapping summary files
    my $s_file_gz = "$analysis_dir"."read_records/$lib"."/$type_dir/$lane"."_"."$type_file".".txt.gz";
    my $s_file_txt = "$analysis_dir"."read_records/$lib"."/$type_dir/$lane"."_"."$type_file".".txt";
    if (-e $s_file_gz){
      $s_count_gz++;
    }
    if (-e $s_file_txt){
      $s_count_txt++;
    }

    #Now check log file for completion status
    my $log_file = "$analysis_dir"."logs/$lib"."/$lane/parseBlastResults_"."$type_name"."_LOG.txt";
    if (-e $log_file){
      $log_count++;

      my $result = `tail $log_file | grep "SCRIPT COMPLETE"`;
      chomp($result);
      if ($result){
        $completed_job_count++;
      }
    }
  }

  my $test1 = 0;
  my $test2 = 0;
  if ($s_count_txt){
    print YELLOW, "\n\tFound $s_count_gz Summary files and $s_count_txt uncompressed Summary files (in process or represent errors?)", RESET;
  }elsif($s_count_gz != $total_lane_count){
    print YELLOW, "\n\tFound $s_count_gz Summary files (expecting: $total_lane_count)", RESET;
  }else{
    print BLUE, "\n\tFound $s_count_gz Summary files (none in process)", RESET;
    $test1 = 1;
  }

  my $message = "Found $log_count log files and $completed_job_count of these indicated that the jobs were completed successfully";
  if ($log_count == $total_lane_count && $completed_job_count == $total_lane_count){
    print BLUE, "\n\t$message", RESET;
    $test2 = 1;
  }else{
    print YELLOW, "\n\t$message", RESET;
  }
  if ($test1 && $test2){
    print MAGENTA, "\n\tDONE", RESET;
  }
}

print "\n\n";

exit();


