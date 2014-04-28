#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#The purpose of this script is to create a fasta file from a read record file containing Solexa paired end reads
#Before using this script first create the read record file for your library from raw solexa seq files
#   i.e. /home/malachig/svn/solexa_analysis/processRawSolexaReads.pl
#The user will specify filtering criteria to exlude reads as a space separated list.
#Any value in this list which is in the status field of a read record will result in that read being skipped
#The user can also specify whether they want to allow one read of a read pair to pass through
# - For certain applications (e.g. assembly with Velvet using paired option) all reads must be present as pairs

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;

my $read_file = '';
my $fasta_file = '';
my $require_paired_reads = '';
my $exclusion_classes = '';
my $log_file = '';

GetOptions ('read_file=s'=>\$read_file, 'fasta_file=s'=>\$fasta_file,
	    'require_paired_reads=s'=>\$require_paired_reads, 'exclusion_classes=s'=>\$exclusion_classes,
	    'log_file=s'=>\$log_file);

print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the input read_file using: --read_file", RESET;
print GREEN, "\n\tSpecify the output fasta file using: --fasta_file", RESET;
print GREEN, "\n\tIf you wish to require that every read have a pair, use: --require_paired_reads=yes", RESET;
print GREEN, "\n\tSpecify the read status values that you wish exclude from the output fasta as a space separated list using: --exclusion_classes", RESET;
print GREEN, "\n\tSpecify the name of a log file using: --log_file", RESET;
print GREEN, "\n\nExample: createSolexaReadFasta.pl  --read_file=/projects/malachig/solexa/read_records/13288AAXX_Lanes1-2_HS03271_ReadRecords.txt  --require_paired_reads=yes  --exclusion_classes=\"Low_Quality Low_Complexity Duplicate\" --fasta_file=/projects/malachig/solexa/fasta_seq_data/13288AAXX_Lanes1-2_QualityFiltered.fa  --log_file=/projects/malachig/solexa/logs/HS03271/createSolexaReadFasta_FILTERED_LOG.txt\n\n", RESET;

unless ($read_file && $fasta_file && ($require_paired_reads =~ /^yes|^y|^n|^no/i) && $exclusion_classes && $log_file){
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  exit();
}

open (LOG, ">$log_file") || die "\nCould not open log file: $log_file\n\n";
print LOG "The is a log for the script: createSolexaReadFasta.pl\n\n";
print LOG "User specified the options: \nread_file = $read_file\nfasta_file = $fasta_file\nrequire_paired_reads = $require_paired_reads\nlog_file = $log_file\n\n";

#Variables to store the number of reads that are filtered at the user specified cutoffs
my $total_reads = 0;
my $total_read1 = 0;
my $total_read2 = 0;
my $read1_passes = 0;
my $read2_passes = 0;
my $read12_passes = 0;

my $counter = 0;

#Deal with the exclusion classes specified by the user
my @exclusion_classes = split(" ", $exclusion_classes);
my %exclusion_classes;
foreach my $c (@exclusion_classes){
  $exclusion_classes{$c}{exclusions} = 0;
}

unless ($read_file =~ /(.*)\.gz$/){
  print RED, "\nInput read file is not compressed!\n\n", RESET;
  exit();
}

print BLUE, "\nParsing read file: $read_file\nWriting fasta file: $fasta_file\n\n", RESET;
print LOG "\nParsing read file: $read_file\n\n";
open (READS, "zcat $read_file |") || die "\nCould not open read file: $read_file\n\n";
open (FASTA, ">$fasta_file") || die "\nCould not open fasta file: $fasta_file\n\n";

if ($require_paired_reads =~ /^yes|^y/i){
  print YELLOW, "\nOnly allowing reads where both reads of the pair pass\n\n", RESET;
}else{
  print YELLOW, "\nAllowing individual reads where the other read of a pair fails\n\n", RESET;
}

my $header = 1;
while(<READS>){
  $counter++;

  if ($counter == 10000){
    $| = 1;
    print BLUE, ".", RESET;
    $| = 0;
    $counter = 0;
  }

  #Skip header line
  if ($header == 1){
    $header = 0;
    next();
  }
  $total_reads+=2;
  $total_read1++;
  $total_read2++;
  chomp($_);
  my @line = split ("\t", $_);

  #Make sure all values are defined!!  Missing values might indicate file corruption at some point
  unless ($line[0] && $line[1] && $line[2] && $line[3] && $line[4] && $line[5] && $line[6] && ($line[7] =~ /\d+/) && ($line[8] =~ /\d+/) && ($line[9] =~ /\d+/) && ($line[10] =~ /\d+/)){
    print RED, "\n\nFound an undefined value in the input read records file!!!\n\n", RESET;
    print RED, "RECORD: $_\n\n", RESET;
    exit();
  }

  my $read1_id = $line[1];
  my $read2_id = $line[2];
  my $read1_status = $line[3];
  my $read2_status = $line[4];
  my $read1_seq = $line[5];
  my $read2_seq = $line[6];

  my $read1_test = 1;
  my $read2_test = 1;

  if ($exclusion_classes{$read1_status}){
    $read1_test = 0;
    $exclusion_classes{$read1_status}{exclusions}++;
  }
  if ($exclusion_classes{$read2_status}){
    $read2_test = 0;
    $exclusion_classes{$read2_status}{exclusions}++;
  }

  if ($require_paired_reads =~ /^yes|^y/i){

    if ($read1_test == 1 && $read2_test == 1){
      $read12_passes++;
      print FASTA ">$read1_id\n$read1_seq\n";
      print FASTA ">$read2_id\n$read2_seq\n";
    }
  }else{

    #If each read passed the tests, print it to the fasta file
    if ($read1_test == 1){
      $read1_passes++;
      print FASTA ">$read1_id\n$read1_seq\n";
    }
    if ($read2_test == 1){
      $read2_passes++;
      print FASTA ">$read2_id\n$read2_seq\n";
    }
    if ($read1_test == 1 && $read2_test == 1){
      $read12_passes++;
    }
  }
}

my $percent_pass = (($read1_passes + $read2_passes)/$total_reads)*100;
my $percent_pass_f = sprintf("%.2f", $percent_pass);

print BLUE, "\n\nThe results of filtering can be summarized as follows:", RESET;
print BLUE, "\n\tTotal reads parsed = $total_reads", RESET;
print BLUE, "\n\tTotal read count:\tRead1 = $total_read1\tRead2 = $total_read2", RESET;
print BLUE, "\n\tTotal passing read count:\tRead1 = $read1_passes\tRead2 = $read2_passes", RESET;
print BLUE, "\n\tTotal percent of all reads (R1 + R2) passing = $percent_pass_f%", RESET;
print BLUE, "\n\tTotal records where both reads passed = $read12_passes\n", RESET;

print LOG "\n\nThe results of filtering can be summarized as follows:";
print LOG "\n\tTotal reads parsed = $total_reads";
print LOG "\n\tTotal read count:\tRead1 = $total_read1\tRead2 = $total_read2";
print LOG "\n\tTotal passing read count:\tRead1 = $read1_passes\tRead2 = $read2_passes";
print LOG "\n\tTotal percent of all reads (R1 + R2) passing = $percent_pass_f%";
print LOG "\n\tTotal records where both reads passed = $read12_passes\n";

print BLUE, "\nExclusion class (# individual reads being filtered out)", RESET;
print LOG "\nExclusion class (# individual reads being filtered out)";
foreach my $ec (sort {$exclusion_classes{$b}->{exclusions} <=> $exclusion_classes{$a}->{exclusions}} keys %exclusion_classes){
  print BLUE, "\n\t$ec ($exclusion_classes{$ec}{exclusions})", RESET;
  print LOG "\n\t$ec ($exclusion_classes{$ec}{exclusions})";
}

close (FASTA);
close (READS);
close (LOG);

#Compress the fasta file
my $compress_cmd = "gzip -f $fasta_file";
print BLUE, "\n\nCompressing the resulting fasta file:  $compress_cmd", RESET;
system ($compress_cmd);

print "\n\nSCRIPT COMPLETE\n\n";

exit();

