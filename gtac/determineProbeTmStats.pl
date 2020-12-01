#!/usr/bin/perl -w

#Written by Malachi Griffith
#This script parses through probe files in a user specified directory and combines probe Tm values from all files into a single array
#The mean, median and other descriptive statistics of this array of probe Tm's is then reported

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

use utilities::Descriptive;

my $probe_dir = '';
my $logfile = '';

GetOptions ('probe_dir=s'=>\$probe_dir, 'logfile=s'=>\$logfile);

#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify a directory containing probe files using: --probe_dir", RESET;
print GREEN, "\n\tAlso specify the name of a log file to be used to store a record of this script: --logfile", RESET;
print GREEN, "\n\nExample: determineMedianTm.pl  --probe_dir=/home/user/alexa/ALEXA_version/unfiltered_probes/  --logfile=/home/user/alexa/ALEXA_version/logs/determineProbeTmStats_LOG.txt\n\n", RESET;

#Make sure all options were specified
unless ($probe_dir && $logfile){
  print RED, "\nRequired input parameter(s) missing!\n\n", RESET;
  exit();
}

open (LOG, ">$logfile") || die "\nCould not open logfile: $logfile\n\n";

#Print out the parameters supplied by the user to the logfile for future reference
print LOG "\nUser Specified the following options:\nprobe_dir = $probe_dir\nlogfile = $logfile\n\n";

#First get all Tm values from the probe files in the specied directory
my $probe_tms_ref = &parseProbeFiles('-dir'=>$probe_dir);

#Go through each collection of probe Tms and generate simple stats
foreach my $group (sort {$a <=> $b} keys %{$probe_tms_ref}){

  my $probe_tm_array = $probe_tms_ref->{$group}->{tms};
  my $probe_count =  @{$probe_tm_array};

  my $stat = Statistics::Descriptive::Full->new();
  $stat->add_data(@{$probe_tm_array});
  my $mean_tm = $stat->mean();
  my $median_tm = $stat->median();
  my $min_tm = $stat->min();
  my $max_tm = $stat->max();
  my $tm_range = $stat->sample_range();

  my $name = $probe_tms_ref->{$group}->{name};
  print BLUE, "\n\nThe following values apply to the $probe_count probe_tms of: $name", RESET;
  print BLUE, "\n\tMin Tm = $min_tm\tMax Tm = $max_tm\tRange = $tm_range", RESET;
  print BLUE, "\n\tMean Tm = $mean_tm\tMedian Tm = $median_tm\n\n", RESET;

  print LOG "\nThe following values apply to the combined list of probe_tms in the input directory:";
  print LOG "\n\tMin Tm = $min_tm\tMax Tm = $max_tm\tRange = $tm_range";
  print LOG "\n\tMean Tm = $mean_tm\tMedian Tm = $median_tm\n\n";
}

close (LOG);

exit();


#################################################################################################################################
#parseProbeFiles - Get probe tms from all files and combine into a single array                                                 #
#################################################################################################################################
sub parseProbeFiles{
  my %args = @_;
  my $probe_dir = $args{'-dir'};

  my %files;
  my %probe_tms;
  my $groups_found = 0;

  #First make sure the specified base path exists and is a directory
  unless (-e $probe_dir && -d $probe_dir){
    print RED, "\nSpecified directory: $probe_dir does not appear valid!\n\n", RESET;
    exit();
  }

  unless ($probe_dir =~ /.*\/$/){
    $probe_dir = "$probe_dir"."/";
  }

  #Get all the input files in the specified directory
  print BLUE, "\nSearching $probe_dir for files\n", RESET;
  print LOG "\nSearching $probe_dir for files\n";

  opendir(DIRHANDLE, "$probe_dir") || die "\nCannot open directory: $probe_dir\n\n";
  my @test_files = readdir(DIRHANDLE);

  my $file_count = 0;
  foreach my $test_file (@test_files){
    my $file_path = "$probe_dir"."$test_file";

    #Skip directories within the specified directory
    if (-e $file_path && -d $file_path){
      print "\n\t$file_path  is a directory - skipping";
      next();
    }
    $file_count++;

    $files{$file_count}{name} = $test_file;
    $files{$file_count}{file_path} = $file_path;
  }

  my @grand_probe_tms;

 FILE:foreach my $file_count (sort {$a <=> $b} keys %files){
    my $file_path = $files{$file_count}{file_path};

    print BLUE, "\n\nProcessing: $file_path", RESET;
    print LOG "\n\nProcessing: $file_path";

    open (INFILE, "$file_path") || die "\nCould not open input probe file $file_path\n\n";

    my $first_line = 1;
    my $probe_tm_column;

    while (<INFILE>){
      chomp($_);

      #Process the header line and determine the column which contains the probe_tm values
      if ($first_line == 1){

	my @header = split("\t", $_);
	
	my $column_found = 0;
	my $column_count = 0;
	foreach my $header (@header){
	  if ($header eq "Probe_Tm"){
	    $column_found = 1;
	    $probe_tm_column = $column_count;
	  }
	  $column_count++;
	}

	#If a probe_tm column was not found, skip this file
	if ($column_found == 1){
	  $groups_found++;
	  my @tms;
	  $probe_tms{$groups_found}{name} = $files{$file_count}{name};
	  $probe_tms{$groups_found}{tms} = \@tms;

	}else{
	  print BLUE, "\nFile: $file_path does not contain a Probe_Tm column - skipping this file\n\n", RESET;
	  print LOG "\nFile: $file_path does not contain a Probe_Tm column - skipping this file\n\n";
	  next FILE;
	}
	$first_line = 0;
	next();
      }

      #Grab the probe tm value from this data line
      my @line = split("\t", $_);
      push(@{$probe_tms{$groups_found}{tms}}, $line[$probe_tm_column]);
      push(@grand_probe_tms, $line[$probe_tm_column]);

    }
    close (INFILE);
    my $current_tm_count = @grand_probe_tms;
    print BLUE, "\nHave $current_tm_count Probe Tm values thus far", RESET;
    print LOG "\nHave $current_tm_count Probe Tm values thus far";
  }

  #Add grand total collection of probe tms seperately
  $groups_found++;
  $probe_tms{$groups_found}{name} = "ALL_PROBES";
  $probe_tms{$groups_found}{tms} = \@grand_probe_tms;

  return(\%probe_tms);
}
