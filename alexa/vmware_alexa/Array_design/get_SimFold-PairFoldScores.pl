#!/usr/bin/perl -w

#Written by Malachi Griffith
#This script gets RNA/DNA pairing scores using the RNAsoft algorithms described by Andronescu et al. 2003 and 
#Andronescu et al. 2004.  Specifically it uses simfold to:
#  "predict the minimum free energy (MFE) secondary structure of a given input RNA or DNA sequence."
#It then uses pairfold to:
#  "predict the MFE secondary structure of two interacting RNA or DNA molecules, and suboptimal structures."
#In this case I am looking for the interaction between a probe and itself.  That is, SimFold finds the within probe 
#folding potential and PairFold finds the potential of two seperate copies of the probe to form a secondary structure
#Note: NimbleGen use in situ DNA synthesis to create their oligonucleotide arrays.  This means that the sequences on the chip 
#are DNA not RNA.  PairFold and SimFold calculate RNA folding energies by default - I use DNA parameters which means that
#the 'dna_or_rna' variables in the MultiRNA-Fold-1.1 code should were changed from 'RNA' to 'DNA'

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;
use Benchmark;

#ALEXA libraries
#When a script is initiated, use the full path of the script location at execution to add the perl module libraries to @INC
#This should allow this scripts to work regardless of the current working directory or the script location (where it was unpacked).
#The /utilities directory must remain in the same directory as this script but the entire code directory can be moved around
BEGIN {
  my $script_dir = &File::Basename::dirname($0);
  push (@INC, $script_dir);
}
use utilities::Probes qw(:all);

my $sequence = '';
my $infile = '';
my $pairfold_dir = '';
my $results_dir = '';
my $file_num = '';

GetOptions ('sequence=s'=>\$sequence, 'infile=s'=>\$infile, 'pairfold_dir=s'=>\$pairfold_dir, 'results_dir=s'=>\$results_dir, 'file_num=i'=>\$file_num);

#Usage instructions
print GREEN, "\nThis function calculates the within self and self-self probe folding scores", RESET;
print GREEN, "\nTo test with a single probe sequence try the following", RESET;
print GREEN, "\nTest usage: get_SimFold-PairFoldScores.pl --sequence ACCCCCTCCTTCCTTGGATCAAGGGGCTCAA --pairfold_dir=~/MultiRNAFold-1.1/\n", RESET;
print GREEN, "\n\tTo generate for a large number of probes, provide a tab-delimited input file containing probe sequences", RESET;
print GREEN, "\n\nUsage: get_SimFold-PairFoldScores.pl  --infile=/home/user/alexa/ALEXA_version/fold/probe_files/probes1  --results_dir=/home/user/alexa/ALEXA_version/fold/results/  --pairfold_dir=~/MultiRNAFold-1.1/  --file_num=1\n\n", RESET;

#Allow test with a single probe
if ($sequence && $pairfold_dir){

  unless ($pairfold_dir =~ /.*\/$/){
    $pairfold_dir = "$pairfold_dir"."/";
  }

  print BLUE, "\nTesting on a single sequence:\n\n", RESET;
  my $simfold_score = &pairFoldCalc ('-sequence1'=>$sequence, '-program'=>'simfold', '-silent'=>0, '-bin_dir'=>$pairfold_dir);
  my $pairfold_score = &pairFoldCalc ('-sequence1'=>$sequence, '-program'=>'pairfold', '-silent'=>0, '-bin_dir'=>$pairfold_dir);
  print "\n\n";
  exit();
}

unless ($infile && $pairfold_dir && $results_dir && $file_num){
  print RED, "\nPlease specify an input file, pairfold_dir, results_dir and file_num\n\n", RESET;
  exit();
}

unless ($pairfold_dir =~ /.*\/$/){
  $pairfold_dir = "$pairfold_dir"."/";
}
unless ($results_dir =~ /.*\/$/){
  $results_dir = "$results_dir"."/";
}

my $params = "$pairfold_dir"."params/"."turner_parameters_fm363_constrdangles.txt";

print BLUE, "\nBegining folding of sequences in $infile\n\n", RESET;

#Execute simfold command
my $t0 = new Benchmark;
my $simfold_result_file = "$results_dir"."results/"."simfold"."$file_num";
my $simfold_bin = "$pairfold_dir"."simfold";
my $cmd1 = "$simfold_bin"." -f $infile -o $simfold_result_file -p $params -v";
print BLUE, "\nSIMFOLD", RESET;
print BLUE, "\nExecuting: $cmd1\n", RESET;
system($cmd1);
my $t1 = new Benchmark;
my $td1 = timediff($t1, $t0);
print YELLOW, "\nThis process took: ",timestr($td1),"\n\n", RESET;

#Execute pairfold command
$t0 = new Benchmark;
my $pairfold_result_file = "$results_dir"."results/"."pairfold"."$file_num";
my $pairfold_bin = "$pairfold_dir"."pairfold";
my $cmd2 = "$pairfold_bin"." -f $infile -o $pairfold_result_file -v";
print BLUE, "\nPAIRFOLD", RESET;
print BLUE, "\nExecuting: $cmd2\n", RESET;
system($cmd2);
$t1 = new Benchmark;
$td1 = timediff($t1, $t0);
print YELLOW, "\nThis process took: ",timestr($td1),"\n\n", RESET;

#Paste input file with second column of simfold and pairfold files
my $result_file1 = "$results_dir"."results/"."pairfold_scores"."$file_num".".tmp";
my $cmd3 = "cut -f 2 $pairfold_result_file > $result_file1";
print BLUE, "\nExecuting: $cmd3\n", RESET;
system($cmd3);

my $result_file2 = "$results_dir"."results/"."probe_scores"."$file_num";
my $cmd4 = "paste $simfold_result_file $result_file1 > $result_file2";
print BLUE, "\nExecuting: $cmd4\n", RESET;
system($cmd4);

#Delete the simfold and pairfold files
my $cmd5 = "rm -f $simfold_result_file $pairfold_result_file $result_file1";
print BLUE, "\nExecuting: $cmd5\n", RESET;
system($cmd5);

exit();
