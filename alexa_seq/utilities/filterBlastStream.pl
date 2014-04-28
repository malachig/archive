#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#Filter BLAST output on the fly via a pipe to of standard out from a BLAST job to this script
#Assumes the BLAST input stream is tab delimited BLAST output
#Filter according to bit_score:

# ** MAKE SURE YOU PICK A REASONABLE BIT SCORE CUTOFF !! **
# Remember that Bit Score depends on the size of the search space (i.e. the chosen WORD LENGTH matters and the database used may also matter):
# To chose a good Bit Score cutoff, run your command without filtering and produce a representative blast output file
# This should be done with reads of the size you want to use and with the word size and database you intend to use
# Then use the script 'summarizeBlastBitScoreProfile.pl' to summarize the Bit Scores resulting for reads of these lengths with different amount of gaps or mismatches
# Examine the output of this script and chose a reasonable cutoff

#For example, if you map 36-mer reads to the EnsEMBL transcriptome (v.49) with a word size of 11
#A BitScore cutoff of 40.0  will allow the following to pass:
#- Perfect matches must be 20 bp or greater
#- Matches with 1 mismatches must be 24 bp or greater
#- Matches with 2 mismatches must be 28 bp or greater
#- Matches with 3 mismatches must be 32 bp or greater
#- Matches with 4 mismatches must be 36 bp or greater
#- Matches with 1 gaps must be 28 bp or greater
#- Matches with 2 gaps must be 36 bp or greater

#Tabular blast output format.  Fields:
#        (0)Query id,(1)Subject id,(2)% identity,(3) alignment length,(4) mismatches,(5) gap openings,
#        (6)q. start, (7)q. end, (8)s. start, (9)s. end, (10)e-value, (11)bit score


#Example usage:
#BLAST COMMAND or cat of blast file | filterBlastStream.pl  --min_bit_score=40.0 > filtered_blast_file.txt

#Include an option that allows the user to require that the hit spans the centre of the target sequence by at least X

use strict;
use Getopt::Long;
use Term::ANSIColor qw(:constants);

my $min_bit_score = '';
my $centre_span = '';
my $target_size = '';

GetOptions ('min_bit_score=f'=>\$min_bit_score, 'target_size=i'=>\$target_size, 'centre_span=i'=>\$centre_span);

unless ($min_bit_score =~ /^\d+/){
  print RED, "\nRequired input parameter(s) missing or incorrect format!\n\n", RESET;
  print RED, "\nExample usage: BLAST COMMAND or cat of blast file | filterBlastStream.pl  --min_bit_score=40.0 > filtered_blast_file.txt\n", RESET;
  print RED, "\nExample usage: BLAST COMMAND or cat of blast file | filterBlastStream.pl  --min_bit_score=40.0  --centre_span=5  --target_size=103 > filtered_blast_file.txt\n", RESET;
  exit();
}

#Sanity check of cutoff format
chomp($min_bit_score);
unless ($min_bit_score =~ m/^\d+.\d+$/ || $min_bit_score =~ m/^\d+$/){
  print "\nfilterBlastStream.pl was supplied a --min_bit_score value that does not appear to be a valid number!\n\n";
  exit();
}


while(<STDIN>){
  my @line = split("\t", $_);
  my $bit_score = $line[11];

  #A.) If the user cares about the read hits spaning the centre of the subject sequence
  if ($centre_span && $target_size){
    my $lower = ($target_size/2)-$centre_span;
    my $upper = ($target_size/2)+$centre_span;

    my $subject_start = $line[8]; #Position on the subject sequence at which the probe alignment begins
    my $subject_end = $line[9]; #Position on the subject sequence at which the probe alignment ends

    #Swap BLAST coords if reversed
    my $original_start;
    my $hit_strand = "+";
    if ($subject_start > $subject_end){
      $original_start = $subject_start;
      $subject_start = $subject_end;
      $subject_end = $original_start;
      $hit_strand = "-";
    }
  
    #Make sure the alignment overlaps the centre of the sequence (+/- the centre_span)
    if ($bit_score >= $min_bit_score && $subject_start <= $lower && $subject_end >= $upper){
      print $_;
    }

  }else{
    #B.) Hits without consideration for spanning the centre of the target sequence
    if ($bit_score >= $min_bit_score){
      print $_;
    }
  }
}

exit();







