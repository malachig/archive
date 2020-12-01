#!/usr/bin/perl -w

#Written by Malachi Griffith
#The purpose of this script is to randomly generate a series of negative control probes to uniformly represent a particular Tm range
#The user will specify the Target Tm, the range around this Tm, the number of probes to keep within each 0.1 degree 'bin' within this range and an output file
#The user must also specify the probe length
#The user will specify the starting probe_id and probeset_id.  In this case a probeset ID corresponds to a single probe_length/tm_bin combination

#Before running this script, determine the probe_length variability and Tm range of all experimental probes (exon-exon, exon-intron, exon, intron)

#In general the neccessary code can be summarized as follows:

#(1.) Get input from the user
#(2.) Randomly generate a probe sequence, calculate it's Tm, place it in the appropriate bin (unless the bin is already) full
#(3.) Repeat step 2 until each bin has the required number of probes
#(4.) Print out a probe file with the usual probe file format (most values will be 'na')
#(5.) Repeat steps 2,3,4 for every probe length in the specified range

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

#Initialize command line options
my $target_length = '';
my $max_length_variance = '';
my $target_tm = '';
my $tm_range = '';
my $bin_size = '';
my $verbose = '';
my $force = '';
my $probe_dir = '';
my $outfile = '';
my $logfile = '';

GetOptions('target_length=i'=>\$target_length, 'max_length_variance=i'=>\$max_length_variance, 'target_tm=f'=>\$target_tm, 'tm_range=f'=>\$tm_range,
	   'bin_size=i'=>\$bin_size, 'verbose=s'=>\$verbose, 'force=s'=>\$force, 'probe_dir=s'=>\$probe_dir, 'outfile=s'=>\$outfile, 
	   'logfile=s'=>\$logfile);

#(1.) Provide instruction to the user
print GREEN, "\n\nThis script randomly generates probes of the specified length and attempts to create a set that uniformly represent a specified Tm range", RESET;
print GREEN, "\nNeccessary options:", RESET;
print GREEN, "\n\tSpecify the desired probe length using: --probe_length (say 36)", RESET;
print GREEN, "\n\tSpecify the probe length range using: --max_length_variance (say 10?)", RESET;
print GREEN, "\n\tSpecify the target Tm (median of probe set) using: --target_tm (say 67.0)", RESET;
print GREEN, "\n\tSpecify the range around this target Tm to allow using: --tm_range (say 8.0?)", RESET;
print GREEN, "\n\tSpecify the number of probes to select within each 0.1 degree bin within the specified range using: --bin_size (say 1000)", RESET;
print GREEN, "\n\tIf you want more detailed output use: --verbose=yes", RESET;
print GREEN, "\n\tTo force selection of current max probe and probeset IDs without prompting use: --force=yes", RESET;
print GREEN, "\n\tSpecify the path to the directory of current probe files (used to determine starting probe ID value) using: --probe_dir", RESET;
print GREEN, "\n\tSpecify the output file to write the resulting probes to using: --output_file (say negativeControlProbes.txt)", RESET;
print GREEN, "\n\tAlso specify the name of a log file to be used to store a record of the design process using: --logfile", RESET;
print GREEN, "\n\n Usage: generate_NegativeControlProbes.pl  --target_length=36  --max_length_variance=10  --target_tm=67.0  --tm_range=3.5  --bin_size=1000  --probe_dir=/home/user/alexa/ALEXA_version/unfiltered_probes/  --outfile=/home/user/alexa/ALEXA_version/unfiltered_probes/negativeControlProbes.txt  --logfile=/home/user/alexa/ALEXA_version/logs/generate_NegativeControlProbes_LOG.txt\n\n", RESET;

#Make sure all options were specified
unless ($target_length && ($max_length_variance || $max_length_variance == 0) && $target_tm && $tm_range && $bin_size && $probe_dir && $outfile && $logfile){
  print RED "\nRequired input parameter(s) missing!\n\n", RESET;
  exit();
}

#Create a series of nucleotide populations to draw random sequences from
#These will range from high GC to high AT.  Hopefully this will help derive sequence near the extremes of the Tm spectrum for each probe length
my $nucleotide_dist_ref = &getNucleotideDistributions();
my $nucleotide_dist_count = keys %{$nucleotide_dist_ref};
my $current_nucleotide_dist = 1;

#Go through each possible length as specified by the user ($target_length +/- $max_length_variance) and generate random probes across the specified Tm range
my $probe_count = 0;
my $total_probe_count = 0;

#Get the starting probe_id and probeset_id by examining the specified input probe file directory.  If it is empty start with 1
my @current_ids;
my $current_probe_id;
my $current_probeset_id;

if ($force){
  if ($force eq "yes"){
    @current_ids = &getCurrentProbeIds('-probe_dir'=>$probe_dir, '-force'=>"yes");
    $current_probe_id = $current_ids[0];
    $current_probeset_id = $current_ids[1];
  }else{
    print RED, "\nForce option not understood (set to 'yes' if you wish to use this option)\n\n", RESET;
    exit();
  }
}else{
  @current_ids = &getCurrentProbeIds('-probe_dir'=>$probe_dir);
  $current_probe_id = $current_ids[0];          #Unique probe ID used for each successful probe
  $current_probeset_id = $current_ids[1];       #Unique count of exon-exon junctions with successful probes
}

open (OUTFILE, ">$outfile") || die "\nCould not open output file: $outfile\n\n";
open (LOG, ">$logfile") || die "\nCould not open log file: $logfile\n\n";

#Print out the parameters supplied by the user to the logfile for future reference
print LOG "\nUser Specified the following options:\ntarget_length = $target_length\nmax_length_variance = $max_length_variance\ntarget_tm = $target_tm\ntm_range = $tm_range\nbin_size = $bin_size\nprobe_dir = $probe_dir\noutfile = $outfile\nlogfile = $logfile\n\n";

#Print the header line
print OUTFILE "Probe_Count\tProbeSet_ID\tGene_ID\tSequence\tProbe_length\tProbe_Tm\tProbe_Type\tExon1_IDs\tUnit1_start\tUnit1_end\tExon2_IDs\tUnit2_start\tUnit2_end\tmasked_bases\n";

print BLUE, "\nPrinting standard output to LogFile: $logfile\n\n", RESET;
print LOG "\nPrinting standard output to LogFile: $logfile\n\n";

srand();

my $progress_count = 0;
for (my $probe_length = $target_length-$max_length_variance; $probe_length <= $target_length+$max_length_variance; $probe_length++){
  my $t0 = new Benchmark;

  #(2.) Randomly generate a probe sequence, calculate it's Tm, place it in the appropriate bin (unless the bin is already) full
  print BLUE, "\nSearching for probes of length $probe_length that evenly cover all Tm bins for $target_tm +/- $tm_range", RESET;
  print LOG "\nSearching for probes of length $probe_length that evenly cover all Tm bins for $target_tm +/- $tm_range";

  #Create bins for the Tm range provided
  my %bins_full;
  my $bins_ref = &generateBins('-target_tm'=>$target_tm, '-tm_range'=>$tm_range);
  my $bin_count = keys %{$bins_ref};

  #Calculate the total number of probes needed to fill all bins to the required level
  my $target_probe_count = ($bin_count * $bin_size);

  print BLUE, "\nFound $bin_count bins.  Searching for $bin_count x $bin_size = $target_probe_count probes\n", RESET;
  print LOG "\nFound $bin_count bins.  Searching for $bin_count x $bin_size = $target_probe_count probes";

  while ($probe_count < $target_probe_count){
    &generateRandomProbe('-probe_length'=>$probe_length, '-target_tm'=>$target_tm, '-tm_range'=>$tm_range, '-bin_size'=>$bin_size,
			 '-bins_full_ref'=>\%bins_full, '-bins_ref'=>$bins_ref);
  }

  $probe_count = 0;

  #(4.)  Print out a probe file with the usual probe file format (most values will be 'na')
  print BLUE, "\nFound $total_probe_count probes so far, printing probes for this target length", RESET;
  print LOG "\nFound $total_probe_count probes so far, printing probes for this target length";

  &printProbeInfo('-probe_object'=>$bins_ref, '-outfile'=>$outfile, '-current_length'=>$probe_length);

  my $t1 = new Benchmark;
  my $td1 = timediff($t1, $t0);
  print YELLOW, "\nThis probe length took: ",timestr($td1),"\n\n", RESET;
  print LOG "\nThis probe length took: ",timestr($td1),"\n\n";

}

print BLUE, "\n\nFound $total_probe_count probes for lengths of $target_length +/- $max_length_variance and Tms of $target_tm +/- $tm_range\n\n", RESET;
print LOG "\n\nFound $total_probe_count probes for lengths of $target_length +/- $max_length_variance and Tms of $target_tm +/- $tm_range\n\n";

close (LOG);
close (OUTFILE);

exit();


###################################################################################################################################################
#Create a series of nucleotide populations to draw random sequences from                                                                          #
#These will range from high GC to high AT.  Hopefully this will help derive sequence near the extremes of the Tm spectrum for each probe length   #
###################################################################################################################################################
sub getNucleotideDistributions{
  my %nucleotide_distributions;

  my $dist_count = 0;
  my $temp_count;

  #High GC population (60% GC, 40% AT)
  my @gc60;
  for (my $x = 0; $x < 30; $x++) {
    push(@gc60, "G");
    push(@gc60, "C");
  }
  for (my $x = 0; $x < 20; $x++) {
    push(@gc60, "A");
    push(@gc60, "T");
  }
  $dist_count++;
  $temp_count = @gc60;
  $nucleotide_distributions{$dist_count}{size} = $temp_count;
  $nucleotide_distributions{$dist_count}{dist} = \@gc60;

  #Uniform population (50% GC, 50% AT)
  my @uni;
  for (my $x = 0; $x < 25; $x++) {
    push(@uni, "G");
    push(@uni, "C");
    push(@uni, "A");
    push(@uni, "T");
  }
  $dist_count++;
  $temp_count = @uni;
  $nucleotide_distributions{$dist_count}{size} = $temp_count;
  $nucleotide_distributions{$dist_count}{dist} = \@uni;

  #Low GC population (40% GC, 60% AT)
  my @gc40;
  for (my $x = 0; $x < 20; $x++) {
    push(@gc40, "G");
    push(@gc40, "C");
  }
  for (my $x = 0; $x < 30; $x++) {
    push(@gc40, "A");
    push(@gc40, "T");
  }
  $dist_count++;
  $temp_count = @gc40;
  $nucleotide_distributions{$dist_count}{size} = $temp_count;
  $nucleotide_distributions{$dist_count}{dist} = \@gc40;

  return(\%nucleotide_distributions);
}


################################################################################################
#Create bins for the Tm range provided                                                         #
################################################################################################
sub generateBins{
  my %args = @_;
  my $tm = $args{'-target_tm'};
  my $range = $args{'-tm_range'};

  my $lower_tm = ($tm - $range);
  my $upper_tm = ($tm + $range);
  my %bins_hash;

  my $i;
  for ($i = $lower_tm; $i <= $upper_tm; $i += 0.1){
    my $x = sprintf("%.1f", $i);
    my %probes;
    $bins_hash{$x}{probes} = \%probes;
  }
  return(\%bins_hash);
}


################################################################################################
#Generate random probe                                                                         #
################################################################################################
sub generateRandomProbe{
  my %args = @_;
  my $length = $args{'-probe_length'};
  my $tm = $args{'-target_tm'};
  my $range = $args{'-tm_range'};
  my $max_bin_size = $args{'-bin_size'};
  my $bins_full_ref = $args{'-bins_full_ref'};
  my $bins_ref = $args{'-bins_ref'};

  #Get the current nucleotide population
  my @nucleotides = @{$nucleotide_dist_ref->{$current_nucleotide_dist}->{dist}};

  $current_nucleotide_dist++;
  if ($current_nucleotide_dist > $nucleotide_dist_count){
    $current_nucleotide_dist = 1;
  }

  #Generate a random probe of the specified length
  my $probe_seq = '';
  for (1..$length){
    my $base = $nucleotides[rand(@nucleotides)];
    $probe_seq = "$probe_seq"."$base";
  }

  #Calculate the Tm of this probe
  my $temp_k = tmCalc ('-sequence'=>$probe_seq, '-silent'=>1);
  my $Tm_celsius = tmConverter ('-tm'=>$temp_k, '-scale'=>'Kelvin');

  #Check if this probe is within the required range
  my $lower_tm = ($tm - $range);
  my $upper_tm = ($tm + $range);
  unless ($Tm_celsius >= $lower_tm && $Tm_celsius <= $upper_tm){
    return();
  }

  my $tm_bin = sprintf("%.1f", $Tm_celsius);

  #Get the current hash of probes for this bin
  my $probes_ref = $bins_ref->{$tm_bin}->{probes};

  #Make sure the bin is not already full;
  my $current_bin_size = keys %{$probes_ref};
  if ($current_bin_size >= $max_bin_size){
    $bins_full_ref->{$tm_bin}->{tmp} = '';
    return();
  }

  #Update the probes hash for this bin
  $probes_ref->{$probe_seq}->{tmp} = '';

  $current_bin_size = keys %{$probes_ref};
  my $full_bins = keys %{$bins_full_ref};

  #Increment the total number of probes selected
  $probe_count++;
  $total_probe_count++;
  $progress_count++;

  if ($verbose eq "yes"){
    print CYAN, "\n\tProbe_count = $probe_count ($full_bins bins full)\tBin: $tm_bin has $current_bin_size probes.  Adding: $probe_seq", RESET;
    #print LOGFILE "\n\tProbe_count = $probe_count ($full_bins bins full)\tBin: $tm_bin has $current_bin_size probes.  Adding: $probe_seq";
  }else{
    if ($progress_count == 1000){
      $progress_count = 0;
      $| = 1;
      print BLUE, ".";
      $| = 0;
    }
  }
  return()
}


################################################################################################
#Given a probe object containing probe sequences for an exon region, print to a probe file.    #
################################################################################################
sub printProbeInfo{
  my %args = @_;
  my $bins_ref = $args{'-probe_object'};
  my $outfile = $args{'-outfile'};
  my $length = $args{'-current_length'};

  foreach my $bin (sort {$a <=> $b} keys %{$bins_ref}){

    $current_probeset_id++;

    #Get the probes for this bin
    my $probes_ref = $bins_ref->{$bin}->{probes};

    foreach my $probe_seq (keys %{$probes_ref}){

      #1.) Negative Control 'Experimental' probe
      #Calculate the Tm of this probe
      my $temp_k = tmCalc ('-sequence'=>$probe_seq, '-silent'=>1);
      my $Tm_celsius = tmConverter ('-tm'=>$temp_k, '-scale'=>'Kelvin');
      my $tm_formatted = sprintf("%.3f", $Tm_celsius);

      #Print the probe info to an output file
      $current_probe_id++;

      print OUTFILE "$current_probe_id\t$current_probeset_id\tna\t$probe_seq\t$length\t$tm_formatted\tControl-Negative\tna\tna\tna\tna\tna\tna\tna\n";
    }
  }
  return();
}
