#!/usr/bin/perl -w
#Written by Obi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#The purpose of this script is to create dummy R2 fastq files
#This is an imperfect hack to allow processing of single-end data

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;

my $file = ''; #fastq file
my $dir = ''; #working dir

GetOptions ('file=s'=>\$file, 'dir=s'=>\$dir);

print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the fastq file to create dummy R2 file for:  --file", RESET;
print GREEN, "\n\nExample:  createDummyR2_fastq.pl  --dir=/scratch/obig/data/pancreas/620AGAAXX --file=s_5_1_sequence.txt.gz\n\n", RESET;

unless ($file && $dir){
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  exit();
}

chomp $file;
chomp $dir;

#Create name for new file
my $dummy_filename;
if ($file=~/s_(\d+)_(\d+)_sequence\.txt/){
  $dummy_filename = "s_"."$1"."_2_"."sequence.txt";
}else{
  print RED, "File name format not recognized: $file\n", RESET;
}

print BLUE, "processing $file\n", RESET;

#If file zipped, unzip it
if ($file =~/\.gz/){
  my $zcat_cmd="zcat $dir/$file > $dir/tempfile";
  print BLUE, "$zcat_cmd\n", RESET;
  system($zcat_cmd);
}else{
  my $cp_cmd="cp $file $dir/tempfile";
  print BLUE, "$cp_cmd\n", RESET;
  system($cp_cmd);
}

#Go through file and create mirroring dummy file
open (R1FILE, "$dir/tempfile") or die "can't open $dir/tempfile\n";
open (R2FILE, ">$dir/$dummy_filename") or die "can't open $dir/$dummy_filename for write\n";

print BLUE, "converting tempfile to $dummy_filename\n", RESET;

#Grab 4 lines at a time as this is the (insane) format of fastq
my @lines;
my $i=0;
while (<R1FILE>){
  push @lines, $_;
  next if (!eof(R1FILE)) and scalar(@lines) <4;
  #process @lines
  $i++;
  my $seq_header = $lines[0];
  my $seq = $lines[1];
  my $q_header = $lines[2];
  my $q = $lines[3];
  chomp ($seq_header, $seq, $q_header, $q);
  my $seq_length=length($seq);

  my $dummy_seq_header = $seq_header;
  $dummy_seq_header=~s/\#0\/1$/\#0\/2/;
  my $dummy_q_header = $q_header;
  $dummy_q_header=~s/\#0\/1$/\#0\/2/;

  my $dummy_seq = '.' x $seq_length;
  my $dummy_q = 'B' x $seq_length;
  print R2FILE "$dummy_seq_header\n$dummy_seq\n$dummy_q_header\n$dummy_q\n";
  #print "$seq_header\n$dummy_seq\n$q_header\n$dummy_q\n";
  #print Dumper(@lines);
  @lines = (); #clear
  if ($i==1000){print "."; $i=0;} #Print out dot every 1000 sequences to track progress
}

close (R1FILE);
close (R2FILE);

#gzip resulting file
my $gzip_cmd = "gzip $dir/$dummy_filename";
print BLUE, "\n$gzip_cmd\n", RESET;
system($gzip_cmd);

#Delete remaining temp file
my $rm_cmd = "rm $dir/tempfile";
print BLUE, "$rm_cmd\n\n", RESET;
system ($rm_cmd);

exit();
