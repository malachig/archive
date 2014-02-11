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
#the 'dna_or_rna' variables in the MultiRNA-Fold-1.8 code should were changed from 'RNA' to 'DNA'
#I also modified this code to accept a file as input rather than a single sequence or pair of sequences

#I recommend using MultiRNAFold v.1.8 or later available from www.rnasoft.ca

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;

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
my $probe_column = '';
my $infile = '';
my $outfile = '';
my $pairfold_dir = '';

GetOptions ('sequence=s'=>\$sequence, 'probe_column=i'=>\$probe_column, 'infile=s'=>\$infile,
	    'outfile=s'=>\$outfile, 'pairfold_dir=s'=>\$pairfold_dir);

#Usage instructions
print GREEN, "\nThis function calculates the within self and self-self probe folding scores", RESET;
print GREEN, "\nTo test with a single probe sequence try the following", RESET;
print GREEN, "\nTest usage: get_SimFold-PairFoldScores.pl --sequence ACCCCCTCCTTCCTTGGATCAAGGGGCTCAA --pairfold_dir=~/MultiRNAFold-1.1/\n", RESET;
print GREEN, "\n\tTo generate for a large number of probes, provide a tab-delimited input file containing probe sequences", RESET;
print GREEN, "\n\tSpecify which column contains the probe sequences, the file will be regenerated with the fold scores appended", RESET;
print GREEN, "\n\nUsage: get_SimFold-PairFoldScores.pl --probe_column=3  --infile=/home/user/alexa/ALEXA_version/unfiltered_probes/exonJunctionProbes.txt  --outfile=/home/user/alexa/ALEXA_version/unfiltered_probes/exonJunctionProbes_foldScores.txt  --pairfold_dir=~/MultiRNAFold-1.1/\n\n", RESET;

#Allow test with a single probe
if ($sequence && $pairfold_dir){
  print BLUE, "\nTesting on a single sequence:\n\n", RESET;
  my $simfold_score = &pairFoldCalc ('-sequence1'=>$sequence, '-program'=>'simfold', '-silent'=>0, '-bin_dir'=>$pairfold_dir);
  my $pairfold_score = &pairFoldCalc ('-sequence1'=>$sequence, '-program'=>'pairfold', '-silent'=>0, '-bin_dir'=>$pairfold_dir);
  print "\n\n";
  exit();
}

unless ($infile && $probe_column && $pairfold_dir && $outfile){
  print RED, "\nPlease specify an input file, probe column, pairfold_dir, and output file\n\n", RESET;
  exit();
}

$probe_column = $probe_column - 1;

my %probes;
my $header = "Probe_Count\tSimFold_score\tPairFold_score\n";

#Get all data from the probe input file
my $header_line = 1;

print BLUE, "\nProcessing file: $infile\n\n", RESET;
open (INFILE, "$infile") || die "\nCould not open $infile";
while (<INFILE>){

  #Deal with the header line
  if ($header_line == 1){
    $header_line = 0;
    next();
  }

  my @arr = split("\t", $_);

  my $probe_id = $arr[0];
  my $probe_seq = $arr[$probe_column];
  $probes{$probe_id}{sequence} = $probe_seq;
}
close (INFILE);

#Location of binaries on each cluster node (the 'local' location of each binary)
#This folder was copied to every node on the cluster by systems
#my $bin_dir = "/tmp/malachi_temp/";
#my $bin_dir = "/opt/MultiRNAFold-1.1/";
my $number_probes = keys %probes;
print BLUE, "\nGetting simfold and pairfold scores for $number_probes probes\n\n", RESET;
foreach my $probe_id (sort {$a <=> $b} keys %probes){

  #Get the simfold and pair fold scores for each probe sequence
  my $seq = $probes{$probe_id}{sequence};
  my $simfold_score = &pairFoldCalc ('-sequence1'=>$seq, '-program'=>'simfold', '-silent'=>1, '-bin_dir'=>$pairfold_dir);
  my $pairfold_score = &pairFoldCalc ('-sequence1'=>$seq, '-program'=>'pairfold', '-silent'=>1, '-bin_dir'=>$pairfold_dir);

  $probes{$probe_id}{simfold_score} = $simfold_score;
  $probes{$probe_id}{pairfold_score} = $pairfold_score;

}

#Print results to output file
print BLUE, "\nPrinting results to outfile: $outfile\n\n", RESET;
open (OUTFILE, ">$outfile") || die "\nCould not open $outfile";

#Print out the probe_ids and each score to the specified output file
print OUTFILE "$header";
foreach my $probe_id (sort {$a <=> $b} keys %probes){
  print OUTFILE "$probe_id\t$probes{$probe_id}{simfold_score}\t$probes{$probe_id}{pairfold_score}\n";
}
close (OUTFILE);

exit();


