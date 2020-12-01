#!/usr/bin/perl -w

#Written by Malachi Griffith
#Take a master file of probe info and divide into into smaller blocks and place in a temp directory in:
#  /home/user/alexa/ALEXA_version/fold/
#Then create a command for the perl script get_SimFold-PairFoldScores.pl specifying this file and the output file desired
#Basically this is just a wrapper for get_SimFold-PairFoldScores which simply uses the SimFold/PairFold binaries

#To run the cluster job, log into your cluster head node
#Make sure the batch file is okay and run the following command:
#mqsub --file $batch_file.sh --name $job_name --mkdir
#This process may differ depending on the cluster you use

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
my $master_probe_file = '';
my $expected_probe_count = '';
my $pairfold_dir = '';
my $folding_bin = '';
my $temp_dir = '';
my $batch_file = '';
my $repair = '';
my $logfile = '';

GetOptions ('master_probe_file=s'=>\$master_probe_file, 'expected_probe_count=i'=>\$expected_probe_count, 'pairfold_dir=s'=>\$pairfold_dir,
	    'folding_bin=s'=>\$folding_bin, 'temp_dir=s'=>\$temp_dir, 'batch_file=s'=>\$batch_file, 'repair=s'=>\$repair,
	    'logfile=s'=>\$logfile);

#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tThis script takes a probe file as input, breaks it into pieces and creates batch jobs for the cluster", RESET;
print GREEN, "\n\tSpecify the input probe file to be processed using: --master_probe_file", RESET;
print GREEN, "\n\tSpecify the expected number of probes using: --expected_probe_count", RESET;
print GREEN, "\n\tSpecify the full path to the directory containing pairfold executables using:  --pairfold_dir (eg. /opt/MultiRNAFold-1.1/)", RESET;
print GREEN, "\n\tSpecify the full path (including name) to the script used to generate folding scores using: --folding_bin (eg. get_SimFold-PairFoldScores.pl)", RESET;
print GREEN, "\n\tSpecify the full path to a temp working directory to store probe files and results using: --temp_dir", RESET;
print GREEN, "\n\tSpecify the name of the resulting batch file to be submitted to the cluster using: --batch_file (eg. pairfold_ExonJunction.sh)", RESET;
print GREEN, "\n\tIf some jobs have failed on the cluster and you simply want to generate a new batch file for only unsuccessful jobs, use: --repair=yes", RESET;
print GREEN, "\n\tAlso specify the name of a log file to be used to store a record of this script: --logfile", RESET;
print GREEN, "\n\nExample: createPairFoldBatch.pl  --master_probe_file=/home/user/alexa/ALEXA_version/unfiltered_probes/exonJunctionProbes.txt  --expected_probe_count=1234  --pairfold_dir=~/MultiRNAFold-1.1/  --folding_bin=~/alexa/Array_design/get_SimFold-PairFoldScores.pl  --temp_dir=/home/user/alexa/ALEXA_version/fold/  --batch_file=/home/user/alexa/ALEXA_version/batch_scripts/pairfold_ExonJunction.sh  --repair=no  --logfile=/home/user/alexa/ALEXA_version/logs/createPairFoldBatch_ExonJunction_LOG.txt\n\n", RESET;

unless ($master_probe_file && $expected_probe_count && $pairfold_dir && $folding_bin && $temp_dir && $batch_file && $logfile && $repair){
  print RED, "\nRequired input parameter(s) missing!\n\n", RESET;
  exit();
}

#Make sure the pairfold_dir actually exists
unless ($pairfold_dir =~ /.*\/$/){
  $pairfold_dir = "$pairfold_dir"."/";
}
unless (-e $pairfold_dir && -d $pairfold_dir){
  print RED, "\nSpecified pairfold directory: $pairfold_dir does not appear valid!\n", RESET;

  print RED, "\nProceed anyway? ", RESET;
  my $answer = <>;
  chomp($answer);
  unless ($answer eq "y" || $answer eq "Y"){
    print "\n\n";
    exit();
  }
}

open (LOG, ">$logfile") || die "\nCould not open logfile: $logfile\n\n";

#If the user specified the --repair option, simply generate a new batch file for those probe files that do not have a successful results file
if ($repair eq "yes"){
  print BLUE, "\nGenerating a new batch file for only those jobs that failed in an earlier submission to the cluster\n\n", RESET;
  &generateNewBatchFile('-temp_dir'=>$temp_dir);
  exit();
}

#Create a $temp_probe_dir called 'probe_files' and a $results_dir called 'results' in the specified working directory
my $temp_probe_dir = &createNewDir('-path'=>$temp_dir, '-new_dir_name'=>"probe_files");
my $results_dir = &createNewDir('-path'=>$temp_dir, '-new_dir_name'=>"results");

#Open the input probe file and output batch file
open (BATCHFILE, ">$batch_file") || die "\nCould not open batch file: $batch_file\n\n";
open (PROBEFILE, "$master_probe_file") || die "\nCould not open probe file: $master_probe_file\n\n";

my $probe_count = 0;
my $total_probe_count = 0;
my $probe_files_printed = 0;
my $block_num = 1;
my $block_size = 20000;
my $probe_info = '';
my ($probe_file,$out_file, $header_line);

#Find out the number of lines in the input file and the column containing probe sequences
my $first_line = 1;
my $probe_sequence_column;
my $num_probes = 0;
while (<PROBEFILE>){

  if ($first_line == 1){
    my @header = split("\t", $_);

    my $column_found = 0;
    my $column_count = 0;
    foreach my $header (@header){
      $column_count++;
      if ($header eq "Sequence"){
	$column_found = 1;
	$probe_sequence_column = $column_count;
      }
    }

    #If a probe_tm column was not found, skip this file
    unless ($column_found == 1){
      print BLUE, "\nFile: $master_probe_file does not contain a probe Sequence column!\n\n", RESET;
      print LOG "\nFile: $master_probe_file does not contain a probe Sequence column!\n\n";
      exit();
    }

    print BLUE, "\nFound probe sequence column in input file.  Column number is: $probe_sequence_column\n\n";
    print LOG "\nFound probe sequence column in input file.  Column number is: $probe_sequence_column\n\n", RESET;

    $first_line = 0;
    next();
  }

  if ($_ =~ /^\d+/){
    $num_probes++;
  }
}

print BLUE, "\nFound a total of $num_probes probes in the master probe file\n\n", RESET;
print LOG "\nFound a total of $num_probes probes in the master probe file\n\n";

#Make sure the right number of probes are found
unless ($num_probes == $expected_probe_count){
  print RED, "\nDid not find the expected number of probes!\n\n", RESET;
  exit();
}

close (PROBEFILE);

open (PROBEFILE, "$master_probe_file") || die "\nCould not open probe file: $master_probe_file";

#open the first output file
$probe_file = "$temp_probe_dir"."probes"."$block_num";
$out_file = "$results_dir"."probe_scores$block_num";

#get_SimFold-PairFoldScores
print BATCHFILE "$folding_bin  --infile $probe_file  --pairfold_dir $pairfold_dir  --results_dir $temp_dir  --file_num $block_num\n";
$probe_files_printed++;

open (TEMP_FILE, ">$probe_file") || die "\nCould not open temp probe file: $probe_file";

my $probes_printed = 0;
$first_line = 1;

while (<PROBEFILE>){
  chomp($_);
  my @line = split("\t", $_);

  #Grab the header line
  if ($first_line == 1){
    print TEMP_FILE "$line[0]\t$line[$probe_sequence_column-1]\n";
    $header_line = "$line[0]\t$line[$probe_sequence_column-1]\n";
    $first_line = 0;
    next();
  }

  $probe_count++;
  $total_probe_count++;

  #print "\rBlock: $block_num\tProbe: $probe_count";
  print TEMP_FILE "$line[0]\t$line[$probe_sequence_column-1]\n";
  $probes_printed++;

  #Watch for end of file
  if ($total_probe_count == $num_probes){
    next();
  }

  #Once a full block is reached start on the next file
  if($probe_count == $block_size){

    #Close the current output file
    close (TEMP_FILE);

    $block_num++;
    $probe_file = "$temp_probe_dir"."probes"."$block_num";

    #Open the new output file
    open (TEMP_FILE, ">$probe_file") || die "\nCould not open temp probe file: $probe_file";

    #Print the header line to this file
    print TEMP_FILE "$header_line";

    #Figure out the name of the new output file
    $out_file = "$results_dir"."probe_scores$block_num";

    #Provide the neccessary command for get_SimFold-PairFoldmFoldDNAScores.pl
    #the -c flag specifies the column containing the probe sequence, the -f flag specifies the input file, the -o flag specifies the output file

    #get_SimFold-PairFoldScores.pl
    print BATCHFILE "$folding_bin  --infile $probe_file  --pairfold_dir $pairfold_dir  --results_dir $temp_dir  --file_num $block_num\n";
    $probe_files_printed++;

    #Reset variables
    $probe_info = '';
    $probe_count = 0;
    $probe_file = '';
    $out_file = '';
  }
}
print BLUE, "\nFound a total of $num_probes probes in the input file", RESET;
print BLUE, "\nPrinted a total of $probes_printed probes", RESET;
print BLUE, "\nThese were printed in blocks of $block_size resulting in a total of $probe_files_printed files\n\n", RESET;

print LOG "\nFound a total of $num_probes probes in the input file";
print LOG "\nPrinted a total of $probes_printed probes";
print LOG "\nThese were printed in blocks of $block_size resulting in a total of $probe_files_printed files\n\n";

close (TEMP_FILE);
close (PROBEFILE);
close (BATCHFILE);
close (LOG);
exit();


###################################################################################################################################################
#If the user specified the --repair option, simply generate a new batch file for those probe files that do not have a successful results file     #
###################################################################################################################################################
sub generateNewBatchFile{
  my %args = @_;
  my $temp_dir = $args{'-temp_dir'};

  unless ($temp_dir =~ /.*\/$/){
    $temp_dir = "$temp_dir"."/";
  }

  my $temp_probe_dir = "$temp_dir"."probe_files"."/";
  my $results_dir = "temp_dir"."results"."/";

  my @probe_files = `dir --format=single-column $temp_probe_dir`;
  my @results_files = `dir --format=single-column $results_dir`;

  my %probe_files;
  foreach my $file (@probe_files){
    chomp($file);
    if ($file =~ /^probes(\d+)/){
      my $file_count = $1;
      $probe_files{$file_count}{file_name} = $file;
    }
  }
  my $probe_files_found = keys %probe_files;
  print BLUE, "\n\nFound $probe_files_found probe files", RESET;

  my %results_files;
  foreach my $file (@results_files){
    chomp($file);
    if ($file =~ /^probe_scores(\d+)/){
      my $file_count = $1;
      $results_files{$file_count}{file_name} = '';
    }
  }
  my $results_files_found = keys %results_files;
  print BLUE, "\nFound $results_files_found results files", RESET;

  #Create batch file entries for all probe files that do not have a corresponding results file
  my $repair_job_count = 0;
  open (BATCHFILE, ">$batch_file") || die "\nCould not open batch file: $batch_file\n\n";
  foreach my $probe_file_count (sort {$a <=> $b} keys %probe_files){

    #unless a results file was generated for this batch of probes, add this job to the new batch file
    unless ($results_files{$probe_file_count}){
      $probe_file = "$temp_probe_dir"."$probe_files{$probe_file_count}{file_name}";
      $out_file = "$results_dir"."probe_scores$probe_file_count";

      print BATCHFILE "$folding_bin  --infile $probe_file  --pairfold_dir $pairfold_dir  --results_dir $temp_dir  --file_num $block_num\n";

      $repair_job_count++;
    }
  }
  close BATCHFILE;

  print BLUE, "\nPrinted $repair_job_count jobs to the new batch file\n\n", RESET;

  return();
}
