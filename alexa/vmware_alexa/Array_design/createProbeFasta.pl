#!/usr/bin/perl -w

#Written by Malachi Griffith
#This script takes examines a specified directory and identifies probe files containing unique probe IDs and probe sequences
#Next a fasta file containing all these probe sequences is generated
#This file can then be used to create a blast database for checking the specificity of probes
#All microarray probes should be unique, even near perfect matches should be avoided
#NOTE: Testing specificity against Negative-Control probes is not neccessary - these probes will be skipped

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);

my $probe_dir = '';
my $outfile = '';
my $logfile = '';

GetOptions ('probe_dir=s'=>\$probe_dir, 'outfile=s'=>\$outfile, 'logfile=s'=>\$logfile);

#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tThis script creates a fasta file of probe sequences probe files in a specified directory", RESET;
print GREEN, "\n\tNegative control probes will not be included!!", RESET;
print GREEN, "\n\tSpecify the full path to a directory containing probes file using: --probe_dir", RESET;
print GREEN, "\n\tSpecify the name of the resulting outfile using: --outfile", RESET;
print GREEN, "\n\tSpecify a logfile for output using: --logfile", RESET;
print GREEN, "\n\nExample: createProbeFasta.pl  --probe_dir=/home/user/alexa/ALEXA_version/unfiltered_probes/  --outfile=/home/user/alexa/ALEXA_version/data_sources/probes_blast_database/all_experimental_probes.fa  --logfile=/home/user/alexa/ALEXA_version/logs/specificity/createProbeFasta_LOG.txt\n\n", RESET;

unless ($probe_dir && $outfile && $logfile){
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  exit();
}

#Open logfile for output
open (LOG, ">$logfile") || die "\nCould not open logfile: $logfile\n\n";
print LOG "\nUser Specified the following options:\nprobe_dir = $probe_dir\noutfile = $outfile\nlogfile = $logfile\n\n";

my %files = &getProbeFiles('-dir'=>$probe_dir);

my %probes;  #Keep track of probe_ids to make sure they are all unique

open (OUTFILE, ">$outfile") || die "\nCould not open output fasta file $outfile\n\n";

#Go through each file and start printing out probe fasta records
my $probe_records_printed = 0;
foreach my $file_count (sort {$a <=> $b} keys %files){

  #Open the input file, check for valid unique IDs, and output as fasta format
  my $file_path = $files{$file_count}{file_path};
  open (INFILE, "$file_path") || die "\nCould not open input probe file $file_path\n\n";

  my $header = 1;
  while (<INFILE>){
    #Skip the header line
    if ($header == 1){
      $header = 0;
      next();
    }

    my @line = split("\t", $_);

    #Process each probe line
    if ($line[0] =~ /^(\d+)/){
      my $probe_id = $1;
      my $probeset_id = $line[$files{$file_count}{probeset_column} - 1];
      my $seq = $line[$files{$file_count}{probe_column} - 1];
      my $type = $line[$files{$file_count}{type_column} - 1];

      #Skip probes that have the type 'Control-Negative'
      if ($type eq "Control-Negative"){
	next();
      }

      #Make sure this probe_id has not been previously observed
      if ($probes{$probe_id}){
	print RED, "\nNon-unique probe id: $probe_id found\n\n", RESET;
	close (OUTFILE);
	close (LOG);
	exit();
      }else{
	$probes{$probe_id}{tmp} = '';
      }

      my $probe_name = "$probe_id"."_"."$probeset_id";

      #Print this probe to the fasta file
      chomp($seq);
      print OUTFILE ">$probe_name\n$seq\n";
      $probe_records_printed++;

    }else{
      print RED, "\nProbe ID format missing or not understood!!\n\n", RESET;
      close (OUTFILE);
      close (LOG);
      exit();
    }
  }
  close (INFILE);
}
print BLUE, "\nPrinted a total of $probe_records_printed fasta probe records to $outfile\n\n", RESET;
print LOG "\nPrinted a total of $probe_records_printed fasta probe records to $outfile\n\n";

close (OUTFILE);
close (LOG);

exit();


#################################################################################################################################
#getProbeFiles                                                                                                                  #
#################################################################################################################################
sub getProbeFiles{
  my %args = @_;
  my $probe_dir = $args{'-dir'};

  my %files;

  unless ($probe_dir =~ /.*\/$/){
    $probe_dir = "$probe_dir"."/";
  }

  #First make sure the specified base path exists and is a directory
  unless (-e $probe_dir && -d $probe_dir){
    print RED, "\nSpecified directory: $probe_dir does not appear valid!\n\n", RESET;
    close (LOG);
    exit();
  }

  #Get all the input files in the repeat masked result directory
  print BLUE, "\nSearching $probe_dir for files\n", RESET;
  print LOG "\nSearching $probe_dir for files\n";

  my %possible_files;
  my $possible_file_count = 0;
  opendir(DIRHANDLE, "$probe_dir") || die "\nCannot open directory: $probe_dir\n\n";
  my @test_files = readdir(DIRHANDLE);

  foreach my $test_file (@test_files){
    my $file_path = "$probe_dir"."$test_file";

    #Skip directories within the specified directory
    if (-e $file_path && -d $file_path){
      print BLUE, "\n\t$file_path  is a directory - skipping", RESET;
      print LOG "\n\t$file_path  is a directory - skipping";
      next();
    }
    $possible_file_count++;

    $possible_files{$possible_file_count}{file_name} = $test_file;
    $possible_files{$possible_file_count}{source_dir} = $probe_dir;
    $possible_files{$possible_file_count}{file_path} = $file_path;
  }

  my $file_num = keys %possible_files;
  print BLUE, "\n\nFound $file_num probe files in the specified directory\n", RESET;
  print LOG "\n\nFound $file_num probe files in the specified directory\n";

  #Check each file for expected columns
  FILE: foreach my $file_count (sort {$a <=> $b} keys %possible_files){
      my $probe_file = $possible_files{$file_count}{file_path};

      open (PROBE_FILE, "$probe_file") || die "\nCould not open probe file: $probe_file\n\n";

      my $probeset_column;
      my $probeset_column_found = 0;
      my $probe_column;
      my $probe_id_column_found = 0;
      my $sequence_column_found = 0;
      my $probe_type_column;
      my $probe_type_column_found = 0;

      #Only process the first line
      while (<PROBE_FILE>){
	my $line_record = $_;
	chomp ($line_record);

	my @line = split ("\t", $line_record);

	#Watch for the header line which is assumed to start with "Probe_Count" and contain "Sequence"
	my $column_count = 0;
	foreach my $value (@line){

	  $column_count++;
	  if ($value =~ /^Probe_Count$/){
	    $probe_id_column_found = 1;
	  }
	  if ($value =~ /^ProbeSet_ID$/){
	    $probeset_column_found = 1;
	    $probeset_column = $column_count;
	  }
	  if ($value =~ /^Sequence$/){
	    $sequence_column_found = 1;
	    $probe_column = $column_count;
	  }
	  if ($value =~ /^Probe_Type$/){
	    $probe_type_column_found = 1;
	    $probe_type_column = $column_count;
	  }
	}

	#If the file has the neccessary columns, add it to the list of files to be processed
	if ($probe_id_column_found == 1 && $probeset_column_found == 1 && $sequence_column_found == 1 && $probe_type_column_found == 1){
	  $files{$file_count}{file_path} = $probe_file;
	  $files{$file_count}{file_name} = $possible_files{$file_count}{file_name};
	  $files{$file_count}{source_dir} = $possible_files{$file_count}{source_dir};
	  $files{$file_count}{probeset_column} = $probeset_column;
	  $files{$file_count}{probe_column} = $probe_column;
	  $files{$file_count}{type_column} = $probe_type_column;

	}else{
	  print YELLOW, "\nFile: $possible_files{$file_count}{file_name} does not appear to be a valid probe file (expected column missing) - skipping\n\n", RESET;
	  print LOG "\nFile: $possible_files{$file_count}{file_name} does not appear to be a valid probe file (expected column missing) - skipping\n\n";
	}
	close (PROBE_FILE);
	next FILE;
      }
    }

  #Summarize the files found - list old and new new names
  print BLUE, "\n\nFile summary:", RESET;
  print LOG "\n\nFile summary:";
  foreach my $file_count (sort {$a <=> $b} keys %files){
    print BLUE, "\n\tFile: $file_count\n\t\tINFILE: $files{$file_count}{file_path}", RESET;
    print LOG "\n\tFile: $file_count\n\t\tINFILE: $files{$file_count}{file_path}";
  }
  print BLUE, "\n\n", RESET;
  print LOG "\n\n";

  return(%files);
}

