#!/usr/bin/perl -w

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use Benchmark;

my $usage=<<INFO;

./createLaneDataLines.pl 

This script starts with a fastq input data dir and a list of library subdirs.  

It searches through these subdirs and identifies lanes or sublanes of data.

For each lane/sublane it prints out a LANE definition line to be used in an ALEXA-seq configuration file

e.g.

./createLaneDataLines.pl  --input_data_dir=/gscmnt/gc2142/techd/analysis/alexa_seq/input_data/  --library_list='M_BW-5289-Wt_M_lot_1-cDNA-1-lib1,M_BW-5290-Wt_M_lot_2-cDNA-1-lib1,M_BW-5291-Wt_M_lot_3-cDNA-1-lib1,M_BW-5292-Wt_F_lot_1-cDNA-1-lib1,M_BW-5293-Wt_F_lot_2-cDNA-1-lib1,M_BW-5294-Wt_F_lot_3-cDNA-1-lib1,M_BW-5295-Null_M_lot_1-cDNA-1-lib1,M_BW-5296-Null_M_lot_2-cDNA-1-lib1,M_BW-5297-Null_M_lot_3-cDNA-1-lib1,M_BW-5298-Null_F_lot_1-cDNA-1-lib1,M_BW-5299-Null_F_lot_2-cDNA-1-lib1,M_BW-5300-Null_F_lot_3-cDNA-1-lib1'

INFO

my $input_data_dir = '';
my $library_list = '';

GetOptions ('input_data_dir=s'=>\$input_data_dir, 'library_list=s'=>\$library_list);

unless ($input_data_dir && $library_list){
  print RED, "\n\nParameters missing\n\n", RESET;
  print GREEN, "$usage", RESET;
  exit();
}

unless (-e $input_data_dir && -d $input_data_dir){
  print RED, "\n\nInput data dir does not appear valid: $input_data_dir\n\n", RESET;
  exit();
}
unless ($input_data_dir =~ /\/$/){
  $input_data_dir .= "/";
}

#Goal is to find all lanes/sublanes of data and create an entry like the following for each of them:
#LANE  MM5289  20BP9ABXX  1  /gscmnt/gc2142/techd/analysis/alexa_seq/input_data/M_BW-5289-Wt_M_lot_1-cDNA-1-lib1/20BP9ABXX_1/  fastq  80  0  2  1  sanger_phred

#Strip white space from string
$library_list =~ s/\s+//g;
my @lib_list = split(",", $library_list);

my %lane_records;
my $lane_record_count = 0;
foreach my $lib_name (sort @lib_list){

  #Extract the library number from the full library name
  my $species_code;
  my $lib_number;
  if ($lib_name =~ /^(\w+)\_\w+\-(\d+)\-/){
    $species_code = $1;
    $lib_number = $2;
  }elsif($lib_name =~ /^(\w+)\_.*\-(\d{7})\-/){
    $species_code = $1;
    $lib_number = $2;
  }elsif($lib_name =~ /^(\w+)\_.*\-(\d{6})\_/){
    $species_code = $1;
    $lib_number = $2;
  }else{
    print RED, "\n\nCould not identify library number from full library name: $lib_name\n\n", RESET;
    exit();
  }

  #Convert species code
  my $species;
  if ($species_code eq "M"){
    $species = "MM";
  }elsif($species_code eq "H"){
    $species = "HS";
  }else{
    print RED, "\n\nUnrecognized species code: $species_code\n\n", RESET;
    exit();
  }

  #Get lane dirs for this library
  my $lib_dir = "$input_data_dir"."$lib_name/";
  opendir(DIR, "$lib_dir") || die "\n\nCould not open lib dir: $lib_dir\n\n";
  my @subdirs = readdir(DIR);
  closedir(DIR);
  print BLUE, "\n$lib_dir", RESET;

  foreach my $subdir (@subdirs){
    if ($subdir =~ /^\./){
      next();
    }
    my $lane;
    my $flowcell;
    if ($subdir =~ /(\w+)\_(\d+)/){
      $flowcell = $1;
      $lane = $2;
    }else{
      print YELLOW, "\n\tLane dir: $subdir not understood - skipping", RESET;
      next();
    }

    #Get lane files for this lane dir
    my $lib_lane_dir = "$lib_dir"."$subdir/";
    print BLUE, "\n\t$flowcell $lane $lib_lane_dir", RESET;
    opendir(DIR, "$lib_lane_dir") || die "\n\nCould not open lib lane dir: $lib_lane_dir\n\n";
    my @lanefiles = readdir(DIR);
    closedir(DIR);

    #If there are sublane files, use these, otherwise use the original files

    my %lanes;
    my $sublanes = 0;
    my $sublane_test = "$lib_lane_dir"."s_"."$lane"."0_1_sequence.txt.bz2";
    if (-e $sublane_test){
      foreach my $file (@lanefiles){
        if ($file =~ /s\_($lane\d+)\_\d+\_sequence\.txt\.bz2/){
          $lanes{$1}=1;
        }
      }
    }else{
      foreach my $file (@lanefiles){
        if ($file =~ /s\_(\d+)\_\d+\_sequence\.txt\.bz2/){
          $lanes{$1}=1;
        }
      }
    }

    foreach my $lane (sort {$a cmp $b} keys %lanes){

      #Get the read length from the top of the file
      my $read_length;
      my $file_path = "$lib_lane_dir"."s_"."$lane"."_1_sequence.txt.bz2";
      my $wc_cmd = "bzcat $file_path | head -n 2 | tail -n 1 | wc";
      my $result = `$wc_cmd`;
      chomp($result);
      if ($result =~ /\d+\s+\d+\s+(\d+)/){
        $read_length = $1 - 1;
      }else{
        print RED, "\n\nCould not determine read length from string: $result ($wc_cmd)", RESET;
        exit();
      }
      #Make sure a positive integer was returned
      unless ($read_length > 0){
        print RED, "\n\nInvalid read length\n\n", RESET;
        exit();
      }


      my $read_file_type = "fastq";
      my $read_trim = 0;
      my $max_n = 2;
      my $min_phred = 1;
      my $qual_type = "sanger_phred";
      my $lane_file = "";
      my $lane_record = "LANE  $species$lib_number  $flowcell  $lane  $lib_lane_dir  $read_file_type  $read_length  $read_trim  $max_n  $min_phred  $qual_type";
      print BLUE, "\n\t\t$lane_record", RESET;
      $lane_record_count++;
      $lane_records{$lane_record_count} = $lane_record;
    }
  }
}

#Print out the final list of lane records
print "\n\n";
foreach my $l (sort {$a <=> $b} keys %lane_records){
  print "\n$lane_records{$l}";
}
print "\n\n";


exit();


