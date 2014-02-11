#!/usr/bin/perl -w

#Written by Malachi Griffith
#The purpose of this script is to simply look for simple repeats.  i.e. I want to identify probe sequences with low complexity
#For example a stretch of AAAAAAA or ATATATATATAT, might be considered problematic.  Some of these very simple repeats can be missed by
#RepeatMasker which looks for LINEs, SINEs, ALUs, MIRs, LTR elements, etc.
#This script is essentially a wrapper for the publicly available low complexity masker: 'mdust'.
#It can be downloaded from http://compbio.dfci.harvard.edu/tgi/software

#This simple C program takes a fasta file as input and prints the output to the screen
#The complexity of the repeat that will be masked by this algorithm can be adjusted using the -v option
#for example..  mdust test.fa -v 28   -> is the default value and will not mask many bases unless the repeat are quite large
#Values less than 10 start to mask stuff in a somewhat haphazard fashion
#A value of 11 or less seems to mask mono-nucleotide repeats of 6 or longer
#A value of 15 or less seems to mask mono-nucleotide repeats of 7 or longer
#A value of 20 or less seems to mask mono-nucleotide repeats of 8 or longer
#A value of 12 or less seems to mask di-nucleotide repeats of 5 or longer (ATATATATAT)
#A value of 17 or less seems to mask di-nucleotide repeats of 6 or longer (ATATATATATAT)
#A value of 22 or less seems to mask di-nucleotide repeats of 7 or longer (ATATATATATATAT)
#A value of 14 or less seems to mask tri-nucleotide repeats of 5 of longer (GATGATGATGATGAT)
#A value of 19 or less seems to mask tri-nucleotide repeats of 6 or longer (GATGATGATGATGATGAT)
#The number of mdust masked bases will be appended to input probe file

#Each of the probe files specified in the input directory will be processed and appended files will be written with the results
#These new files will have the same name as the orginal file except '_mdust' will be added

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
use utilities::utility qw(:all);

#Initialize command line options
my $mdust_bin = '';
my $probe_dir = '';
my $temp_dir = '';
my $logfile = '';

GetOptions ('mdust_bin=s'=>\$mdust_bin, 'probe_dir=s'=>\$probe_dir, 'temp_dir=s'=>\$temp_dir,
	    'logfile=s'=>\$logfile);

#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tThis script processes each of the probe files in a specified directory, by identifying low complexity regions of probe sequences", RESET;
print GREEN, "\n\tSpecify the full path to your mdust binary using: --mdust_bin", RESET;
print GREEN, "\n\tSpecify the full path to a directory containing input probe files to be processed using: --probe_dir", RESET;
print GREEN, "\n\t\tOutput files will be generated automatically and will have the same names with _mdust appended", RESET;
print GREEN, "\n\tSpecify the full path to a temp directory for writing temp files using: --temp_dir", RESET;
print GREEN, "\n\tAn mdust cutoff value of 11 will be used (11 is very stringent and will mask a mononucleotide repeat of 6 or longer)", RESET;
print GREEN, "\n\t\tA value of 27 is the mdust default, and will only mask much longer repetitive sequences", RESET;
print GREEN, "\n\t\tSee script comments for more detailed examples", RESET;
print GREEN, "\n\tSpecify a logfile for output using: --logfile", RESET;
print GREEN, "\n\nExample: test_ProbeComplexity.pl  --mdust_bin=/home/user/BioSw/dust/mdust/mdust  --probe_dir=/home/user/alexa/ALEXA_version/unfiltered_probes/  --temp_dir=/home/user/alexa/ALEXA_version/mdust  --logfile=/home/user/alexa/ALEXA_version/logs/test_ProbeComplexity_LOG.txt\n\n", RESET;

unless ($mdust_bin && $probe_dir && $temp_dir && $logfile){
  print RED, "\nRequired input parameter(s) missing!\n\n", RESET;
  exit();
}

my $mdust_cutoff = 11;

#Open logfile for output
open (LOG, ">$logfile") || die "\nCould not open logfile: $logfile\n\n";
print LOG "\nUser Specified the following options:\nmdust_bin = $mdust_bin\nprobe_dir = $probe_dir\nmdust_cutoff = $mdust_cutoff\ntemp_dir = $temp_dir\nlogfile = $logfile\n\n";

#Create a $temp_dir called 'mdust'
my $working_dir = &createNewDir('-path'=>$temp_dir, '-new_dir_name'=>"mdust");

my %files = &getProbeFiles('-dir'=>$probe_dir);

foreach my $file_count (sort {$a <=> $b} keys %files){

  #Specify a directory where probe files can be temporarily be created and fed into mdust
  my $mdust_tempfile = "$working_dir"."/"."temp.fa";

  #Open the input file
  my $probe_file = $files{$file_count}{file_path};
  open (PROBE_FILE, "$probe_file") || die "\nCould not open $probe_file\n\n";

  #Open the output file
  my $outfile = $files{$file_count}{new_file_path};
  open (OUT_FILE, ">$outfile") || die "\nCould not open $outfile\n\n";

  my $line_count = 0;
  my $block_size = 1000;
  my $block_count = 0;

  print BLUE, "\n\nProcesing file: $probe_file\n", RESET;
  print LOG "\n\nProcesing file: $probe_file\n";

  print BLUE, "\n\nProcesing input probe sequences with $mdust_bin (One '.' = $block_size probes tested)\n", RESET;
  print LOG "\n\nProcesing input probe sequences with $mdust_bin (One '.' = $block_size probes tested)\n";

  #Go through the probefile and get the probe sequences to be processed in blocks
  my %probes;
  while (<PROBE_FILE>){
    $line_count++;
    #Skip the header line
    if ($line_count == 1){
      chomp($_);
      print OUT_FILE "$_\tMdustBases_cutoff_$mdust_cutoff\n";
      next();
    }

    my @line = split("\t", $_);

    #Get the probe ID and sequence from each line
    my $id = $line[0];
    my $probe_sequence = $line[$files{$file_count}{probe_column} - 1];

    $probes{$id}{sequence} = $probe_sequence;
    chomp($_);
    $probes{$id}{line_info} = $_;

    #Once a full set of probes has been collected process them with mdust
    if ((keys %probes) == $block_size){
      $block_count++;
      print BLUE, ".", RESET;
      print LOG ".";

      &processProbes('-probe_hash'=>\%probes, '-mdust_tempfile'=>$mdust_tempfile);

      %probes = ();
    }
  }
  #Process the final incomplete block of probes 
  &processProbes('-probe_hash'=>\%probes, '-mdust_tempfile'=>$mdust_tempfile);

  print BLUE, "\n\n", RESET;
  print LOG "\n\n";

  close (PROBE_FILE);
  close (OUT_FILE);
}

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

      my $probe_column = 0;
      my $probe_id_column_found = 0;
      my $sequence_column_found = 0;
      my $mdust_column_found = 0;

      #Only process the first line
      while (<PROBE_FILE>){
	my $line_record = $_;
	chomp ($line_record);

	my @line = split ("\t", $line_record);

	#Watch for the header line which is assumed to start with "Probe_Count" and contain "Sequence"
	my $probe_column_count = 0;
	foreach my $value (@line){

	  $probe_column_count++;
	  if ($value =~ /^Probe_Count$/){
	    $probe_id_column_found = 1;
	  }
	  if ($value =~ /^Sequence$/){
	    $sequence_column_found = 1;
	    $probe_column = $probe_column_count;
	  }
	  if ($value =~ /Mdust/){
	    $mdust_column_found = 1;
	  }
	}
	#If the file already has mdust scores, skip it
	if ($mdust_column_found == 1){
	  print YELLOW, "\nFile: $possible_files{$file_count}{file_name} already has mdust scores - skipping\n\n", RESET;
	  print LOG "\nFile: $possible_files{$file_count}{file_name} already has mdust scores - skipping\n\n";
	  close (PROBE_FILE);
	  next FILE;
	}

	#If the file has the neccessary columns, add it to the list of files to be processed
	if ($probe_id_column_found == 1 && $sequence_column_found == 1){
	  $files{$file_count}{file_path} = $probe_file;
	  $files{$file_count}{file_name} = $possible_files{$file_count}{file_name};
	  $files{$file_count}{source_dir} = $possible_files{$file_count}{source_dir};
	  $files{$file_count}{probe_column} = $probe_column;

	  #Create a new file name and path
	  my $new_file_name = "$files{$file_count}{file_name}"."_mdust";
	  $files{$file_count}{new_file_name} = $new_file_name;
	  $files{$file_count}{new_file_path} = "$probe_dir"."$new_file_name";

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
    print BLUE, "\n\tFile: $file_count\n\t\tINFILE: $files{$file_count}{file_path}\n\t\tOUTFILE: $files{$file_count}{new_file_path}", RESET;
    print LOG "\n\tFile: $file_count\n\t\tINFILE: $files{$file_count}{file_path}\n\t\tOUTFILE: $files{$file_count}{new_file_path}";
  }
  print BLUE, "\n\n", RESET;
  print LOG "\n\n";

  return(%files);
}


###################################################################################################
#Process Probes with mdust
###################################################################################################
sub processProbes{
  my %args = @_;
  my $probe_list_ref = $args{'-probe_hash'};
  my $mdust_tempfile = $args{'-mdust_tempfile'};

  #Create a temp fasta file with the probe sequences in this block
  open (FASTA_FILE, ">$mdust_tempfile") || die "\nCould not open $mdust_tempfile\n\n";
  foreach my $probe_id (sort {$a <=> $b} keys %{$probe_list_ref}){
    print FASTA_FILE ">$probe_id\n$probe_list_ref->{$probe_id}->{sequence}\n";
  }
  close (FASTA_FILE);

  my $result = `$mdust_bin $mdust_tempfile -v $mdust_cutoff`;

  #Split the resulting output back into single probe entries
  my @dusted_results = split (">", $result);

  #Remove the first entry which is empty
  shift(@dusted_results);

  #Grab the probe ID and sequence again
  foreach my $entry (@dusted_results){
    if ($entry =~ /^(\d+)(.*)/s){
      my $probe_id = $1;
      my $seq = $2;
      chomp($seq);

      #clean-up sequence data by removing spaces, >,  and making uppercase
      $seq =~ s/[\s\>]//g;
      my $upper_seq = uc($seq);

      #Count the number of N's
      my @dusted_n_count = $upper_seq =~ m/N/g;
      my $dusted_n_count = @dusted_n_count;

      $probe_list_ref->{$probe_id}->{n_count} = $dusted_n_count;
      $probe_list_ref->{$probe_id}->{dusted_sequence} = $seq;

    }else{
      print RED, "\nCould not understand output format of entry from mdust!\n\n", RESET;
      close (LOG);
      exit();
    }
  }

  #Print out the probe info and append the mdust cutoff and number of dust masked N's (representing simple repeats) found
  foreach my $probe_id (sort {$a <=> $b} keys %{$probe_list_ref}){
    print OUT_FILE "$probe_list_ref->{$probe_id}->{line_info}\t$probe_list_ref->{$probe_id}->{n_count}\n";
  }
  return();
}
