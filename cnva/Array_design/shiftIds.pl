#!/usr/bin/perl -w
#Written by Malachi Griffith

#The purpose of this script is to shift Probe and ProbeSet IDs so that they start at the specified number


use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;

#Initialize command line options
my $probe_file = '';
my $first_probe_id = '';
my $first_probeset_id = '';
my $out_file = '';

GetOptions ('probe_file=s'=>\$probe_file, 'out_file=s'=>\$out_file,
	    'first_probe_id=i'=>\$first_probe_id, 'first_probeset_id=i'=>\$first_probeset_id);

print GREEN, "\n\nExample: shiftIds.pl  --probe_file=randomControlProbes.txt  --first_probe_id=15  --first_probeset_id=5  --out_file=randomControlProbes_shifted.txt\n\n", RESET;

unless ($probe_file && ($first_probe_id =~ /\d+/) && ($first_probeset_id =~ /\d+/) && $out_file){
  print RED, "\nRequired input parameter(s) missing!\n\n", RESET;
  exit();
}


open (PROBE, "$probe_file") || die "\nCould not open input file: $probe_file\n\n", RESET;
open (OUT, ">$out_file") || die "\nCould not open output file: $out_file\n\n", RESET;

my $current_probe_id = 0;
my $current_probeset_id = 0;
$first_probeset_id--;

while(<PROBE>){
  my $line = $_;
  chomp($line);
  my @line = split("\t", $line);

  #watch for header line
  unless ($line[0] =~ /^\d+/){
    print OUT "Probe_Count\tProbeSet_ID\t$line\n";
    next();
  }

  my $probe_id = $line[0];
  my $probeset_id = $line[1];

  if ($probeset_id == $current_probeset_id){

  }else{
   $current_probeset_id = $probeset_id;
   $first_probeset_id++;
  }

  print OUT "$first_probe_id\t$first_probeset_id\t$line\n";
  $first_probe_id++;

}



close (PROBE);
close (OUT);

exit();

