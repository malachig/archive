#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#Take UCSC mRNA/EST/xmRNA/xEST files and partition into pieces by chromosome

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;
use IO::File;

#Load the ALEXA modules
BEGIN {
  use Cwd 'abs_path';
  if (abs_path($0) =~ /(.*)\/.*\/.*\.pl/){
    push (@INC, $1);
  }
}
use utilities::utility qw(:all);

my $ucsc_align_dir = '';
my $out_dir = '';

GetOptions ('ucsc_align_dir=s'=>\$ucsc_align_dir,'out_dir=s'=>\$out_dir);

print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify a directory containing UCSC mRNA/EST/xmRNA and xEST files using: --ucsc_align_dir", RESET;
print GREEN, "\n\tSpecify an output directory for partitioned files using:  --out_dir", RESET;
print GREEN, "\n\nExample: partitionUcscAlignmentFiles.pl  --ucsc_align_dir=/projects/malachig/sequence_databases/hs_53_36o/mrna_est/   --out_dir=/projects/malachig/sequence_databases/hs_53_36o/mrna_est/partitions/\n\n", RESET;

unless ($ucsc_align_dir && $out_dir){
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  exit();
}

#Check the temp dir before proceeding.  Make sure it is empty
$ucsc_align_dir = &checkDir('-dir'=>$ucsc_align_dir, '-clear'=>"no");
$out_dir = &checkDir('-dir'=>$out_dir, '-clear'=>"yes");

#Get the UCSC alignment files (one each for target species mRNA and EST and one each for all other species mRNA and EST)
my $ucsc_mrna_table = "$ucsc_align_dir"."all_mrna.txt.gz";
my $ucsc_est_table = "$ucsc_align_dir"."all_est.txt.gz";
my $ucsc_xmrna_table = "$ucsc_align_dir"."xenoMrna.txt.gz";
my $ucsc_xest_table = "$ucsc_align_dir"."xenoEst.txt.gz";

unless (-e $ucsc_mrna_table || -e $ucsc_est_table || -e $ucsc_xmrna_table || -e $ucsc_xest_table){
  print RED, "\n\nCould not find any of the following UCSC alignment files needed (looking in: $ucsc_align_dir)\n\t$ucsc_mrna_table\n\t$ucsc_est_table\n\t$ucsc_xmrna_table\n\t$ucsc_xest_table\n\n", RESET;
  exit();
}

if (-e $ucsc_mrna_table){
  &partitionUcscAlignments('-align_file'=>$ucsc_mrna_table, '-seq_type'=>"mrna", '-out_dir'=>$out_dir);
}else{
  print YELLOW, "\n\nDid not find: $ucsc_mrna_table in $ucsc_align_dir - skipping\n\n", RESET;
}

if (-e $ucsc_est_table){
  &partitionUcscAlignments('-align_file'=>$ucsc_est_table, '-seq_type'=>"est", '-out_dir'=>$out_dir);
}else{
  print YELLOW, "\n\nDid not find: $ucsc_est_table in $ucsc_align_dir - skipping\n\n", RESET;
}

if (-e $ucsc_xmrna_table){
  &partitionUcscAlignments('-align_file'=>$ucsc_xmrna_table, '-seq_type'=>"xmrna", '-out_dir'=>$out_dir);
}else{
  print YELLOW, "\n\nDid not find: $ucsc_xmrna_table in $ucsc_align_dir - skipping\n\n", RESET;
}

if (-e $ucsc_xest_table){
  &partitionUcscAlignments('-align_file'=>$ucsc_xest_table, '-seq_type'=>"xest", '-out_dir'=>$out_dir);
}else{
  print YELLOW, "\n\nDid not find: $ucsc_xest_table in $ucsc_align_dir - skipping\n\n", RESET;
}


exit();


##############################################################################################################################################
#Partition alignments                                                                                                                        #
##############################################################################################################################################
sub partitionUcscAlignments{
  my %args = @_;
  my $file = $args{'-align_file'};
  my $type = $args{'-seq_type'};
  my $out_dir = $args{'-out_dir'};

  print BLUE, "\n\nPartitioning UCSC alignments of the type: $type (from $file)", RESET;

  #Open a file foreach chromosome name encountered and write records to these files according to their chromosome
  my %chr_list;
  my $record_count = 0;
  open (ALIGN, "zcat $file |") || die "\nCould not open UCSC alignment file: $file\n\n";
  my $open_fh = 0;
  while (<ALIGN>){
    chomp($_);
    my @line = split("\t", $_);
    my $chr = $line[14];
    my $file_name = "$out_dir"."$chr"."_"."$type".".txt";
    $chr_list{$chr}{file_name} = "$file_name";

    #Unless this chr was already encountered 
    #unless ($chr_list{$chr}){
      #$open_fh++;
      #print YELLOW, "\n$open_fh: Opening new file handle for: $file_name", RESET;
      #my $fh = IO::File->new(">$file_name") || die "\nCould not open output file: $file_name\n\n";
      #$chr_list{$chr}{fh} = $fh;
    #}
    $record_count++;
    #my $fh = $chr_list{$chr}{fh};
    #print $fh "$_\n";

    #Print without opening too many file handles
    open(TEMP, ">>$file_name") || die "\nCould not open temp file: $file_name";
    print TEMP "$_\n";
    close(TEMP);
  }
  close(ALIGN);

  my $chr_count = keys %chr_list;
  print BLUE, "\n\tPrinted $record_count total records corresponding to $chr_count unique chromosome names", RESET;

  #Now close all the file handles
  #foreach my $chr (keys %chr_list){
    #my $fh = $chr_list{$chr}{fh};
    #close($fh);
  #}

  #Now compress all the resulting files
  print BLUE, "\n\tCompressing files\n", RESET;
  foreach my $chr (keys %chr_list){
    my $path = $chr_list{$chr}{file_name};
    my $cmd = "gzip $path";
    system($cmd);
  }

  return();
}







