#!/usr/bin/perl -w
#Written by Obi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#The purpose of this script is to create dummy R2 qseq files
#This is an imperfect hack to allow processing of single-end data

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;

my $file_dir = '';   #Directory of qseq files

GetOptions ('file_dir=s'=>\$file_dir);

print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the dir containing qseq files to create dummy R2 files for:  --file_dir", RESET;
print GREEN, "\n\nExample:  createDummyR2_qseq.pl  --file_dir=/scratch/obig/data/CPTAC_ovary/\n\n", RESET;

unless ($file_dir){
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  exit();
}

#Go through the input dir and get a list of qseq files
print BLUE, "\nGetting list of qseq files from: $file_dir\n\n", RESET;
my %jobs;

my @filelist = `ls $file_dir`;

#print Dumper (@filelist);

foreach my $file (@filelist) {
  chomp $file;
  #Use a regular expression to ignore files beginning with a period
  next if ($file =~ m/^\./);
  next unless ($file =~ /_qseq/);

  #Create name for new file
  my $dummy_filename;
  if ($file=~/s_(\d+)_(\d+)_(\d+)_qseq\.txt/){
    $dummy_filename = "s_"."$1"."_2_"."$3"."_qseq.txt";
  }else{
     print BLUE, "File name format not recognized: $file\n", RESET;
  }

  my $filepath = "$file_dir/"."$file";
  print "processing $filepath\n";

  #If file zipped, unzip it
  if ($file =~/\.gz/){
    my $zcat_cmd="zcat $filepath > $file_dir/tempfile";
    print BLUE, "$zcat_cmd\n", RESET;
    system($zcat_cmd);
  }else{
    my $cp_cmd="cp $filepath $file_dir/tempfile";
    print BLUE, "$cp_cmd\n", RESET;
    system($cp_cmd);
  }

  #Go through file and create mirroring dummy file
  open (R1FILE, "$file_dir/tempfile") or die "can't open $file_dir/tempfile\n";
  open (R2FILE, ">$file_dir/$dummy_filename") or die "can't open $file_dir/$dummy_filename for write\n";
    while (<R1FILE>){
      my @data = split ("\t", $_);
      my $seq = $data[8];
      my $q = $data[9];
      my $seq_length=length($seq);
      my $dummy_seq = '.' x $seq_length;
      my $dummy_q = 'B' x $seq_length;
      print R2FILE "$data[0]\t$data[1]\t$data[2]\t$data[3]\t$data[4]\t$data[5]\t$data[6]\t2\t$dummy_seq\t$dummy_q\t$data[10]"
    }
  close (R1FILE);
  close (R2FILE);
  
  #gzip resulting file
  my $gzip_cmd = "gzip $file_dir/$dummy_filename";
  print BLUE, "$gzip_cmd\n\n", RESET;
  system($gzip_cmd);
}

#Delete remaining temp file
my $rm_cmd = "rm $file_dir/tempfile";
system ($rm_cmd);

exit();
