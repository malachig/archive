#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#This script joins read record files that were generated on a lane by lane basis into one master read record file with a single header line
#Note that a read record file is assumed to have the format:
#Read_ID Read1_ID        Read2_ID        Read1_Status    Read2_Status    Read1_Seq       Read2_seq       Read1_Ambig_Count       Read2_Ambig_Count     Read1_mdust_bases        Read2_mdust_bases

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);

my $read_record_dir = '';
my $joined_file = '';
my $expected_read_count = '';

GetOptions ('read_record_dir=s'=>\$read_record_dir, 'joined_file=s'=>\$joined_file, 'expected_read_count=i'=>\$expected_read_count);

print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify a directory containing compressed read record files to be joined using: --read_record_dir", RESET;
print GREEN, "\n\tSpecify the name of the output joined file using: --joined_file", RESET;
print GREEN, "\n\tSpecify the expected read count for all files using: --expected_read_count", RESET;
print GREEN, "\n\nExample: joinReadRecordFiles.pl  --read_record_dir=/projects/malachig/solexa/read_records/HS04391/  --joined_file=/projects/malachig/solexa/read_records/HS04391/HS04391_Lanes1-8.txt  --expected_read_count=22370010\n\n", RESET;

unless ($read_record_dir && $joined_file && ($expected_read_count =~ /\d+/)){
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  exit();
}

#Check input dir
unless ($read_record_dir =~ /.*\/$/){
  $read_record_dir = "$read_record_dir"."/";
}
unless (-d $read_record_dir && -e $read_record_dir){
  print RED, "\nInput directory does not appear to be a valid directory:\n\tread_record_dir = $read_record_dir\n\n", RESET;
  exit();
}


opendir(DIRHANDLE, "$read_record_dir") || die "\nCannot open directory: $read_record_dir\n\n";
my @files = readdir(DIRHANDLE);
closedir(DIRHANDLE);
my %read_files;

my $file_count = 0;

foreach my $file (@files){

  #Skip sub-directories
  if (-d $file){
    next();
  }

  #Skip all files except .txt.gz files with the correct name format
  my $flowcell;
  my $lane;
  if ($file =~ /([a-zA-Z0-9]+)\_Lane(\d+)\.txt\.gz$/){
    $flowcell = $1;
    $lane = $2;
  }elsif($file =~ /([a-zA-Z0-9]+)\_Lane(\d+)\_.*\.txt\.gz$/){
    $flowcell = $1;
    $lane = $2;
  }else{
    print YELLOW, "\nFile: $file does not appear to be a read record file (or file name format is not understood ...)\n\n", RESET;
    next();
  }

  my $file_path = "$read_record_dir"."$file";

  #Open each file, read the first line and make sure it has a valid header
  open(READS, "zcat $file_path |") || die "\nCould not open file: $file_path\n\n";

  while(<READS>){
    chomp($_);
    my @line = split ("\t", $_);

    if ($line[0] =~ /Read_ID/){
      $file_count++;
      $read_files{$file_count}{file_name} = $file;
      $read_files{$file_count}{flowcell} = $flowcell;
      $read_files{$file_count}{lane} = $lane;
      $read_files{$file_count}{file_path} = $file_path;
      $read_files{$file_count}{flowcell_lane} = "$flowcell"."_"."$lane";

    }else{
      print YELLOW, "\nFile: $file does not appear to be a read record file (or header format is not understood ...)\n\n", RESET;
    }

    last(); #read only the first line
  }
  close(READS);
}

#Print out the list of files to be joined
my $files_found = keys %read_files;
unless ($files_found > 0){
  print YELLOW, "\nDid not find any files to be processed - aborting\n\n", RESET;
  exit();
}

print BLUE, "\n\nFound the following files to be joined:\n", RESET;
foreach my $file_count (sort {$read_files{$a}->{flowcell_lane} cmp $read_files{$b}->{flowcell_lane}} keys %read_files){
  print BLUE, "\n\t$read_files{$file_count}{file_name}", RESET;
}
print "\n\n";
print BLUE, "\n\nDoes this seem correct (y/n)? ", RESET;
my $answer = <>;
chomp($answer);
unless ($answer =~ /y|yes/i){
  print RED, "\nAborting then ...\n\n", RESET;
  exit();
}

#Now actually join the list of files, including the header from only the first file
print BLUE, "\n\nBegin joining files... to $joined_file\n", RESET;
my $first_file = 1;
open (JOINED, ">$joined_file") || die "\nCould not open output file: $joined_file\n\n";

my $total_reads_printed = 0;
foreach my $file_count (sort {$read_files{$a}->{flowcell_lane} cmp $read_files{$b}->{flowcell_lane}} keys %read_files){
  print BLUE, "\n\tProcessing: $read_files{$file_count}{file_name}", RESET;

  open (READS, "zcat $read_files{$file_count}{file_path} |") || die "\nCould not open read file: $read_files{$file_count}{file_path}\n\n";
  while(<READS>){
    if (($_ =~ /^Read_ID/) && ($first_file == 1)){
      print JOINED "$_";
      $first_file = 0;
      next();
    }elsif(($_ =~ /^Read_ID/) && ($first_file == 0)){
      next();
    }
    print JOINED "$_";
    $total_reads_printed++;
  }
  close (READS);
}
close (JOINED);

if ($total_reads_printed == $expected_read_count){
  print BLUE, "\nTotal reads printed ($total_reads_printed) matches the expected number supplied by the user ($expected_read_count)\n\n", RESET;
}else{
  print RED, "\nTotal reads printed ($total_reads_printed) does not match the expected number supplied by the user ($expected_read_count)\n\n", RESET;
}


exit();
