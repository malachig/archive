#!/usr/bin/perl -w

#Written by Malachi Griffith
#This script parse probe files in a specified directory, grabs probe sequences and determines if all probe sequences present are unique.
#If a probe sequence is not unique, it is noted and related info is printed to the screen

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);

my $probe_dir = '';
my $logfile = '';

GetOptions ('probe_dir=s'=>\$probe_dir, 'logfile=s'=>\$logfile);

#Usage instructions
print GREEN, "\nThis searches for non-unique probe sequences in probe files in the specified directory", RESET;
print GREEN, "\nThis script assumes that the first column will contain a unique probe id", RESET;
print GREEN, "\n\tSpecify the full path to a directory containing probe files using: --probe_dir", RESET;
print GREEN, "\n\tSpecify a logfile for output using: --logfile", RESET;
print GREEN, "\n\nUsage: checkProbeUniqueness  --probe_dir=/home/user/alexa/ALEXA_version/filtered_probes  --logfile=/home/user/alexa/ALEXA_version/logs/checkProbeUniqueness_LOG.txt\n\n";

unless ($probe_dir && $logfile){
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  exit();
}
open (LOG, ">$logfile") || die "\nCould not open logfile: $logfile\n\n";
print LOG "\nUser Specified the following options:\nprobe_dir = $probe_dir\nlogfile = $logfile\n\n";

my %sequences;

#Get the files to be processed
my %files = &getProbeFiles('-dir'=>$probe_dir);

my $probe_count = 0;
foreach my $file_count (sort {$a <=> $b} keys %files){

  #Open the input file
  my $probe_file = $files{$file_count}{file_path};
  open (PROBE_FILE, "$probe_file") || die "\nCould not open $probe_file\n\n";

  my %columns = %{$files{$file_count}{columns}};

  #Parse through a file containing probe sequences, store relevant info and note non-unique sequences
  while (<PROBE_FILE>){

    my $line_record = $_;
    chomp ($line_record);

    my @line = split("\t", $line_record);

    my $probe_id = $line[0];

    my $relevant_data = "$line[$columns{'Probe_Count'}{column_pos}]\t$line[$columns{'ProbeSet_ID'}{column_pos}]\t$line[$columns{'Gene_ID'}{column_pos}]\t$line[$columns{'Sequence'}{column_pos}]\t$line[$columns{'Probe_Type'}{column_pos}]";

    if ($probe_id =~ /\d+/){
      $probe_count++;
      my $probe_seq = $line[$columns{'Sequence'}{column_pos}];
      chomp ($probe_seq);

      if ($sequences{$probe_seq}){
	my $probe_ids_ref = $sequences{$probe_seq}{probe_ids};
	$probe_ids_ref->{$probe_id}->{lr} = $relevant_data;
      }else{
	my %probe_ids;
	$probe_ids{$probe_id}{lr} = $relevant_data;
	$sequences{$probe_seq}{probe_ids} = \%probe_ids;
      }
    }
  }
  close (PROBE_FILE);
}

#Summarize findings from parsing
my $unique_sequences = keys %sequences;
print BLUE, "\nFound $probe_count Probes and $unique_sequences Unique Probe Sequences\n\n", RESET;
print LOG "\nFound $probe_count Probes and $unique_sequences Unique Probe Sequences\n\n";

#Now actually display all the non-unique probes
my $non_unique_count = 0;
my $non_unique_cluster_count = 0;

print BLUE, "\nThe following sequences are not unique:\n\n", RESET;
print LOG "\nThe following sequences are not unique:\n\n";

foreach my $seq (keys %sequences){

  my $probe_ids_ref = $sequences{$seq}{probe_ids};
  my $count = keys %{$probe_ids_ref};

  if ($count > 1){
    $non_unique_cluster_count++;
    print BLUE, "\nNon-unique probe cluster: $non_unique_cluster_count", RESET;
    print LOG "\nNon-unique probe cluster: $non_unique_cluster_count";
    foreach my $probe_id (sort {$a <=> $b} keys %{$probe_ids_ref}){
      $non_unique_count++;
      print BLUE, "\n\t$probe_ids_ref->{$probe_id}->{lr}", RESET;
      print LOG "\n\t$probe_ids_ref->{$probe_id}->{lr}";
    }
  }
}
print BLUE, "\n\nFound a total of $non_unique_count non-unique probes from an input list of $probe_count probes\n\n", RESET;
print LOG "\n\nFound a total of $non_unique_count non-unique probes from an input list of $probe_count probes\n\n";

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

      my %columns;

      #Only process the first line
      while (<PROBE_FILE>){
	my $line_record = $_;
	chomp ($line_record);

	my @line = split ("\t", $line_record);

	#Watch for the header line which is assumed to start with "Probe_Count" and contain "Sequence"
	my $column_count = 0;
	foreach my $column (@line){

	  $columns{$column}{column_pos} = $column_count;
	  $column_count++;
	}

	#If the file has the neccessary columns, add it to the list of files to be processed
	if ($columns{'Probe_Count'} && $columns{'ProbeSet_ID'} && $columns{'Gene_ID'} && $columns{'Sequence'} && $columns{'Probe_Type'}){
	  $files{$file_count}{file_path} = $probe_file;
	  $files{$file_count}{file_name} = $possible_files{$file_count}{file_name};
	  $files{$file_count}{source_dir} = $possible_files{$file_count}{source_dir};
	  $files{$file_count}{columns} = \%columns;

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
