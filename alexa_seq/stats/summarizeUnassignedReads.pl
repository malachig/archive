#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#Go through all read record files in the specified root directory
#Build a list of read ids where the status is current 'Unassigned'.  Initialize the value for this read ID as '0'
#Store the grand total read count
#Next go through all mapping results files and if the read ID is found, update the read ids list so that its value is '1'
#Go through the list and determine the number of Unassigned reads that did and did not appear in the mapping results

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;
use Benchmark;

#Load the ALEXA modules
BEGIN {
  use Cwd 'abs_path';
  if (abs_path($0) =~ /(.*)\/.*\/.*\.pl/){
    push (@INC, $1);
  }
}
use utilities::utility qw(:all);

#Initialize command line options
my $read_records_dir = '';

GetOptions ('read_records_dir=s'=>\$read_records_dir);

#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the root directory to read records for a single library using: --read_record_dir", RESET;
print GREEN, "\n\nExample: summarizeUnassignedReads.pl  --read_records_dir=/projects/malachig/solexa/read_records/HS04391/\n\n", RESET;

#Make sure all options were specified
unless ($read_records_dir){
  print RED, "\nRequired input parameter(s) missing!\n\n", RESET;
  exit();
}

$read_records_dir = &checkDir('-dir'=>$read_records_dir, '-clear'=>"no");

#Get read records files
my @read_files = `ls $read_records_dir*.txt.gz`;
my %read_files;
my $count = 0;
foreach my $file (@read_files){
  $count++;
  chomp($file);
  $read_files{$count}{file} = $file;
}

#Get mapping results files
my @map_files = `ls $read_records_dir*/*.txt.gz`;
my %map_files;
my %type_counts;
$count = 0;
foreach my $file (@map_files){
  $count++;
  chomp($file);
  $map_files{$count}{file} = $file;

  #Store the file type - to keep track of where the similarities are coming from:
  my $type;
  if ($file =~ /Lane\d+\_(\w+)\_v53\.txt\.gz/){
    $type = $1;
  }elsif($file =~ /Lane\d+\_(\w+)\.txt\.gz/){
    $type = $1;
  }else{
    print RED, "\n\nFile name not understood: $file\n\n",RESET;
    exit();
  }
  $map_files{$count}{type} = $type;
  $type_counts{$type} = 0;
}

#Build a list of all unassigned reads
my %ureads;
print BLUE, "\n\nGetting list of unassigned quality reads", RESET;
my $grand_read_count = 0;
foreach my $c (sort {$read_files{$a}{file} cmp $read_files{$b}{file}} keys %read_files){

  my $file = $read_files{$c}{file};
  print BLUE, "\nProcessing: $file", RESET;
  open (READS, "zcat $file |") || die "\nCould not open file: $file\n\n";

  my $header = 1;
  while(<READS>){
    if ($header == 1){
      $header = 0;
      next();
    }
    chomp($_);
    my @line = split("\t", $_);

    my $r1_id = $line[1];
    my $r2_id = $line[2];
    my $r1_status = $line[3];
    my $r2_status = $line[4];

    #Count quality reads toward grand total
    unless ($r1_status =~ /Duplicate|Low_Complexity|Low_Quality/){
      $grand_read_count++; 
    }
    unless ($r2_status =~ /Duplicate|Low_Complexity|Low_Quality/){
      $grand_read_count++; 
    }

    #If a read is still unassigned, store it
    if ($r1_status eq "Unassigned"){
      $ureads{$r1_id} = 1;
    }
    if ($r2_status eq "Unassigned"){
      $ureads{$r2_id} = 1;
    }
  }
  close(READS);
  my $current_uread_count = keys %ureads;
  print BLUE, "\n\tGrand reads so far = $grand_read_count\tGrand unassigned reads so far = $current_uread_count", RESET;
}

my $grand_uread_count = keys %ureads;

#Now go through all the map results and mark reads that are found in these
print BLUE, "\n\nChecking for all unassigned reads in the map files", RESET;

foreach my $c (sort {$map_files{$a}{file} cmp $map_files{$b}{file}} keys %map_files){

  my $file = $map_files{$c}{file};
  my $type = $map_files{$c}{type};
  print BLUE, "\nProcessing: $file", RESET;
  open (MAP, "zcat $file |") || die "\nCould not open file: $file\n\n";

  my $header = 1;
  my %columns;
  while(<MAP>){
    chomp($_);
    my @line = split("\t", $_);

    if ($header == 1){
      my $column_count = 0;
      foreach my $column (@line){
        $columns{$column}{column_pos} = $column_count;
        $column_count++;
      }
      $header = 0;
      next();
    }

    #Get the read id(s) for this line
    my @read_ids;
    if ($columns{'R1_ID'} && $columns{'R2_ID'}){
      my $r1_id = $line[$columns{'R1_ID'}{column_pos}];
      my $r2_id = $line[$columns{'R2_ID'}{column_pos}];
      push(@read_ids, $r1_id);
      push(@read_ids, $r2_id);

    }elsif($columns{'Read_ID'}){
      my $r_id = $line[$columns{'Read_ID'}{column_pos}];
      push(@read_ids, $r_id);

    }else{
      print RED, "\n\nRead ID column names not understood for file: $file\n\n", RESET;
      exit();
    }
    #Check each read ID from this line against the master list
    foreach my $rid (@read_ids){
      if (undef($ureads{$rid})){
        next();
      }else{
        delete($ureads{$rid});
        $type_counts{$type}++;
      }
    }
  }
}

#Now count the final result
my $grand_mapped_count = keys %ureads;

print BLUE, "\n\nFINAL RESULTS:", RESET;
print BLUE, "\nGrand total quality reads = $grand_read_count", RESET;
print BLUE, "\nGrand total unassigned quality reads = $grand_uread_count", RESET;
print BLUE, "\nGrand total of these with some similarity to database sequences = $grand_mapped_count", RESET;

print BLUE, "\n\nSimilarities found (by type):", RESET;
foreach my $type (sort {$a cmp $b} keys %type_counts){
  print BLUE, "\n\t$type: $type_counts{$type}", RESET;
}
print BLUE, "\n\n", RESET;




exit();












