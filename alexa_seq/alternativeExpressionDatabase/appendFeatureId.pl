#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#This script assigns helps to assign a unique feature ID to a file of features
#Since many feature are extracted by jobs that run on the cluster, this is difficult to do at the time of feature creation
#Each feature ID will consist of a prefix code followed by a number (from 1 to the total number of features)
#The user MUST specify a prefix that is unique to each feature type:
#For example:
# G = genes
# T = transcripts
# ER = exon regions
# EJ = exon junctions
# EB = exon boundaries
# IN = introns
# AIN = active intron region
# SIN = silent intron region
# IG = intergenic region
# AIG = active intergenic region
# SIG = silent intergenic region

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;

#Initialize command line options
my $infile = '';
my $outfile = '';
my $prefix = '';

GetOptions ('infile=s'=>\$infile,'outfile=s'=>\$outfile, 'prefix=s'=>\$prefix);

print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify a compressed feature database file as input using: --infile", RESET;
print GREEN, "\n\tSpecify a name for the updated output database file using:  --outfile", RESET;
print GREEN, "\n\tSpecify the prefix to use for the feature ID using:  --prefix", RESET;
print GREEN, "\n\nExample: appendFeatureId.pl   --infile=genes/genes_annotated.txt.gz  --outfile=genes/tmp.txt --prefix='G'\n\n", RESET;

#Make sure all options were specified
unless ($infile && $outfile && $prefix){
  print RED, "\nRequired input parameter(s) missing!\n\n", RESET;
  exit();
}

my $count = 0;

#Check input file
unless ((-e $infile) && ($infile =~ /\.gz$/)){
  print RED, "\n\nInput file was not found or is not compressed!!\n\n", RESET;
  exit();
}

print BLUE, "\n\nProcessing: $infile, appended FID values (with prefix: $prefix) and printing to $outfile\n\n", RESET;

#Open input and output files
open (IN, "zcat $infile |") || die "\n\nCould not open input file: $infile\n\n";
open (OUT, ">$outfile") || die "\n\nCould not open output file: $outfile\n\n";

my $header = 1;
my $c = 0;
while(<IN>){
  chomp($_);
  my $line = $_;
  if ($header){
    $header = 0;
    if ($line =~ /FID/){
      print RED, "\nFile ($infile) already contains an FID column!\n\n", RESET;
      close(IN);
      close(OUT);
      system("rm -f $outfile");
      exit();
    }
    print OUT "$line\tFID\n";
    next();
  }
  $c++;
  my $fid = "$prefix$c";
  print OUT "$line\t$fid\n";
}
close(IN);
close(OUT);

print BLUE, "\n\nAppended FID values to $c records\n\n", RESET;

#Compress output file
my $cmd1 = "gzip -f $outfile";
print BLUE, "\n\nCompressing:\n$cmd1\n\n", RESET;
system($cmd1);


#Overwrite original file with new file
my $new_outfile = "$outfile".".gz";
my $cmd2 = "mv $new_outfile $infile";
print BLUE, "\n\nReplacing old file with new appended file:\n$cmd2\n\n", RESET;
system($cmd2);

exit();

