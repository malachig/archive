#!/usr/bin/perl -w

#Written by Malachi Griffith
#This script parse through probe files, gets probe sequences and calculates the cycles required for synthesis

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;

my $probe_file = '';
my $out_file = '';

GetOptions ('probe_file=s'=>\$probe_file, 'out_file=s'=>\$out_file);

#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify a probe files using: --probe_file", RESET;
print GREEN, "\n\tSpecify a new name for the probe file after values are appended using:  --out_file", RESET;
print GREEN, "\n\nExample: calculateNimbleGenSynthesisCycles.pl  --probe_file=/projects/malachig/GTAC_Chip_Design/unfiltered_probes/negativeControlProbes.txt  --out_file=/projects/malachig/GTAC_Chip_Design/unfiltered_probes/negativeControlProbes_cycles.txt\n\n", RESET;


#Make sure all options were specified
unless ($probe_file && $out_file){
  print RED, "\nRequired input parameter(s) missing!\n\n", RESET;
  exit();
}


open (INFILE, "$probe_file") || die "\nCould not open input probe file $probe_file\n\n";
open (OUT, ">$out_file") || die "\nCould not open output probe file $out_file\n\n";

my $first_line = 1;
my $probe_seq_column;

while (<INFILE>){
  chomp($_);
  my $line = $_;
  my @line = split("\t", $line);

  #Process the header line and determine the column which contains the probe_tm values
  if ($first_line == 1){

    print OUT "$line\tNimbleGenCycles\n";

    my $column_found = 0;
    my $column_count = 0;
    foreach my $header (@line){
      if ($header eq "Sequence"){
	$column_found = 1;
	$probe_seq_column = $column_count;
      }
      $column_count++;
    }

    #If a probe_tm column was not found, skip this file
    unless ($column_found == 1){
      print RED, "\nFile: $probe_file does not contain a sequence column!\n\n", RESET;
      exit();
    }
    $first_line = 0;
    next();
  }

  my $cycles = &calculate_nimblegen_cycles('-sequence'=>$line[$probe_seq_column]);

  print OUT "$line\t$cycles\n";

}

close(OUT);
close(INFILE);

exit();


####################################################################################################
#The following subroutine was provided by NimbleGen                                                #
#Copyright 2001-2002, NimbleGen Systems Inc. All rights reserved.                                  #
####################################################################################################
sub calculate_nimblegen_cycles {
  my %args = @_;
  my $probe = $args{'-sequence'};

  my $synthesis_order = "ACGT";
  my @order = split("",$synthesis_order);

  # grab the probe sequence from the subroutine call and clean it up
  $probe = uc $probe;
  $probe =~ tr/AGCT//c;
  my @probe = split("",$probe);

  my $cycles = 0;
  # we know we have to do for every bp in the oligo
  for my $position (0 .. $#probe) {
    # an 'N' means skip 4 cycles.
    if($probe[-1] eq 'N') {
      pop @probe;
      $cycles+=4;
      if(scalar(@probe) == 0) {
	return $cycles;
      }
      next;
    }
    foreach my $base (@order) {
      $cycles++;
      # examine the 3' bp. If it equals the bp we're currently
      # making, pop it off the stack. When we get to the end
      # return the current cycle count
      if ($probe[-1] eq $base) { pop @probe };
      if (scalar(@probe) == 0) {
	return ($cycles);
      }
    }
  }
}
