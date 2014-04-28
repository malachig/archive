#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#Simple script to collapse duplicate reads
#This is a hack to reduce the dataset size for input to the Velvet assembler which used obscene amounts of memory and disk space when processing large numbers of reads

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;

use lib '/home/malachig/perl/BerkeleyDB_x64';
use BerkeleyDB;

#Initialize command line options
my $input_fasta = '';
my $output_fasta = '';
my $filter_n_reads = '';
my $min_length = '';

GetOptions ('input_fasta=s'=>\$input_fasta, 'output_fasta=s'=>\$output_fasta, 'filter_n_reads=i'=>\$filter_n_reads, 'min_length=i'=>\$min_length);

#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the path to an input fasta (compressed) using:  --input_fasta", RESET;
print GREEN, "\n\tSpecify the path to an output fasta (will be compressed) using:  --output_fasta", RESET;
print GREEN, "\n\tTo filter reads containing N's use: --filter_n_reads=1", RESET;
print GREEN, "\n\tTo filter reads less than a particular length use: --min_length=42", RESET;
print GREEN, "\n\nExample: collapseDuplicateReads.pl  --input_fasta=/projects/malachig/solexa/fasta_seq_data/HS04391/HS04391_Lanes1-23_QualityFiltered_Unpaired_1.fa.gz  --output_fasta=/projects/malachig/solexa/velvet_analysis/SecondAttempt/HS04391/unpaired/HS04391_Lanes1-23_QualityFiltered_Unpaired_Collapsed.fa  --filter_n_reads=1  --min_length=42\n\n", RESET;

unless ($input_fasta && $output_fasta){
  print RED, "\nRequired input parameter(s) missing!\n\n", RESET;
  exit();
}

my $pid = $$;

my %h;

my $block_size = 1000000;
my $n_reads = 0;
my $short_reads = 0;

print BLUE, "\n\nBegin processing input, storing unique read sequences\n\n", RESET;
open(IN, "zcat $input_fasta |") || die "\nCould not open input file: $input_fasta\n\n";
my $c = 0;
my $gc = 0;
my $last = 1;
while(<IN>){
  chomp($_);
  if ($_ =~ /^\>/){
    next();
  }

  if ($filter_n_reads){
    if ($_ =~ /N/i){
      $n_reads++;
      next();
    }
  }
  if ($min_length){
    unless (length($_) >= $min_length){
      $short_reads++;
      next();
    }
  }

  $c++;
  $gc++;
  if($c == $block_size){
    my $message = &getMemUsage;
    $c = 0;
    my $unique_reads = keys %h;
    my $diff = ($unique_reads-$last)/$block_size;
    $last = $unique_reads;
    my $percent_diff = sprintf("%.2f", ($diff)*100);
    print BLUE, "\n\t$gc reads processed, $unique_reads reads stored, $percent_diff% new ($message)", RESET;
  }

  if ($h{$_}){
    $h{$_}++;
  }else{
    $h{$_}=1;
  }
}
close(IN);

print BLUE, "\n\nRemoved $n_reads reads containing N's and $short_reads reads that were too short", RESET;

print BLUE, "\n\nBegin printing output to: $output_fasta\n\n", RESET;
open (OUT, ">$output_fasta") || die "\n\nCould not open output file: $output_fasta for writing\n\n";
$c = 0;
$gc = 0;
my $rid = 0;
foreach my $seq (keys %h){
  $c++;
  $gc++;
  $rid++;
  my $count = $h{$seq};
  my $read_name = "$rid"."_"."$count";
  print OUT ">$read_name\n$seq\n";

  if($c == $block_size){
    my $message = &getMemUsage;
    $c = 0;
    print BLUE, "\n\t$gc reads printed", RESET;
  }
}
close(OUT);

#Compress the output file
print BLUE, "\n\nCompressing fasta file\n\n", RESET;
my $cmd = "gzip -f $output_fasta";
system($cmd);

exit();


###########################################################################################################################################
#getMemoryUsage                                                                                                                           #
###########################################################################################################################################
sub getMemUsage{
  my $ps_query = `ps -p $pid -o pmem,rss`;
  my @process_info = split ("\n", $ps_query);
  my $memory_usage = '';
  my $memory_usage_p = '';
  if ($process_info[1] =~ /(\S+)\s+(\S+)/){
    $memory_usage_p = $1;
    $memory_usage = $2;
  }
  my $memory_usage_m = sprintf("%.1f", ($memory_usage/1024));
  my $message = "Memory usage: $memory_usage_m Mb ($memory_usage_p%)"; 

  return($message)
}
