#!/usr/bin/perl -w

#Written by Malachi Griffith
#This script simply looks in a directory for blast files, parses through each of them and summarizes the number of
#hits between probes and and sequences.
#It will be used to determine the hits of negative control probes to mRNA,EST or Ensembl Sequences.
#It will then summarize the number of hits per target and output the results as an updated probe file
#Since negative probes do not have a target, the number of hits will be stored as 'non-target-hits'

#Use to summarize EST, mRNA, and Ensembl BLAST results
use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);

#Initialize command line options
my $probe_file = '';
my $blast_file = '';
my $blast_type = '';
my $outfile = '';

GetOptions('probe_file=s'=>\$probe_file, 'blast_file=s'=>\$blast_file, 'blast_type=s'=>\$blast_type, 'outfile=s'=>\$outfile);

#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify a probe file containing negative control probes using: --probe_file", RESET;
print GREEN, "\n\tSpecify a blast result file using: --blast_file", RESET;
print GREEN, "\n\tSpecify a blast type (EST, mRNA, or Ensembl) using: --blast_type", RESET;
print GREEN, "\n\tSpecify an output probe file using: --outfile", RESET;
print GREEN, "\n\nExample: summarizeBlastResults_NegativeProbes.pl --probe_file ../probe_dev/negativeControlProbes.txt --blast_file ../probe_dev/blast/blast_results/blast1.txt --blast_type mRNA --outfile negativeControlProbes_spec.txt\n\n", RESET;

#Specify the blast results file.  User must also specify the source file containing the probes that were BLASTed
unless($probe_file && $blast_file && $blast_type && $outfile){
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  exit();
}

#Determine the type of blast results to summarize
unless ($blast_type eq "EST" || $blast_type eq "mRNA" || $blast_type eq "enst" || $blast_type eq "probe" || $blast_type eq "genomic"){
  print RED, "\nBlast type not understood!\n\n", RESET;
  exit();
}

my $non_target_hits_col = "$blast_type"."_Non-TargetHits";
my $target_hits_col = "$blast_type"."_TargetHits";
my $non_target_hits_length_col = "$blast_type"."_largestNon-TargetHitLength";
my $target_hits_length_col = "$blast_type"."_largestTargetHitLength";

#Use a hash to store the number of hits for each probe and the longest hit observed
my %probe_hits;

#Parse all the blast results files to find the number of probe hits

print BLUE, "\nProcessing $blast_file", RESET;

open (BLAST_FILE, "$blast_file") || die "\nCould not open blast file: $blast_file\n\n";

while (<BLAST_FILE>){

  my @line = split ("\t", $_);

  my $probe_id = $line[0];
  my $alignment_length = $line[3];

  #See if any hits have already been recorded for this probe
  if ($probe_hits{$probe_id}){
    $probe_hits{$probe_id}{hit_count}++;
    if ($alignment_length > $probe_hits{$probe_id}{max_align_length}){
      $probe_hits{$probe_id}{max_align_length} = $alignment_length;
    }
  }else{
    $probe_hits{$probe_id}{hit_count} = 1;
    $probe_hits{$probe_id}{max_align_length} = $alignment_length;
  }
}
close (BLAST_FILE);

#Go through each line in the probe file and append results based on the parsing of the blast files
open (PROBES, "$probe_file") || die "\nCould not open input probe file: $probe_file\n\n";
open (OUTFILE, ">$outfile") || die "\nCould not open output file: $outfile\n\n";

my $header_line = 1;
while (<PROBES>){
  chomp($_);

  my @line = split("\t", $_);
  my $probe_id = $line[0];

  if ($header_line == 1){
    $header_line = 0;
    print OUTFILE "Probe_Count\t$non_target_hits_col\t$target_hits_col\t$non_target_hits_length_col\t$target_hits_length_col\n";
    next();
  }

  #Check if there are any hits of this probe to a sequence

  #If there are no hits, all values are written as 0
  unless ($probe_hits{$probe_id}){
    print OUTFILE "$probe_id\t0\t0\t0\t0\n";
    next();
  }

  #If there was at least one hit, report the hits as non-target hits (since a negative control probe has no target)
  print OUTFILE "$probe_id\t$probe_hits{$probe_id}{hit_count}\t0\t$probe_hits{$probe_id}{max_align_length}\t0\n";
}

close (PROBES);
close (OUTFILE);

exit();
