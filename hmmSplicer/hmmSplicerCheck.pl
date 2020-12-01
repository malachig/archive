#!/usr/bin/perl -w
#Written by Malachi Griffith and Obi Griffith
#Copyright 2010 Malachi Griffith

#Load in a list of hmmSplicer jobs
#Check for completion of each and produce a new list of jobs for those that are still incomplete

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;
use IO::File;
use Benchmark;

my $jobs_file = '';
my $new_jobs_file = '';
my $clean = '';

GetOptions ('jobs_file=s'=>\$jobs_file, 'new_jobs_file=s'=>\$new_jobs_file, 'clean=s'=>\$clean);

if ($jobs_file && $new_jobs_file){
  #Parameters ok
}else{
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  print GREEN, "\nOnly failed jobs will be added to the repair file.  To include pending jobs as well use:  --clean=yes", RESET;
  print GREEN, "\n\nExample: hmmSplicerCheck.pl  --jobs_file=/projects/alexa2/hmmSplicer/jobs/hmmSplicerJobs.sh  --new_jobs_file=/projects/alexa2/hmmSplicer/jobs/hmmSplicerJobs_REPAIR.sh\n\n", RESET;
  exit();
}

open(IN, "$jobs_file") || die "\n\nCould not open input file: $jobs_file\n\n";
open(OUT, ">$new_jobs_file") || die "\n\nCould not open output file: $new_jobs_file\n\n";

print BLUE, "\n\nChecking jobs present in: $jobs_file", RESET;
my $pass_count = 0;
my $fail_count = 0;
my $pending_count = 0;
my $total = 0;
while(<IN>){
  my $line = $_;
  chomp($line);

  unless($line =~ /\w+/){
    next();
  }
  if ($line =~ /^\#/){
    next();
  }

  #Capture neccessary details from the command: --library_id=SA017  --flowcell_lane_read=30N93AAXX_4_1  --analysis_dir=/projects/alexa2/hmmSplicer/
  my $lib_id;
  my $flowcell_lane_read;
  my $analysis_dir;

  if ($line =~ /library\_id\=(\S+)/){
    $lib_id = $1;
  }else{
    print RED, "\n\nCould not identify library ID from line:\n\n$line\n\n", RESET;
    exit();
  }
  if ($line =~ /flowcell\_lane\_read\=(\S+)/){
    $flowcell_lane_read = $1;
  }else{
    print RED, "\n\nCould not identify flowcell lane read from line:\n\n$line\n\n", RESET;
    exit();
  }
  if ($line =~ /analysis\_dir\=(\S+)/){
    $analysis_dir = $1;
  }else{
    print RED, "\n\nCould not identify analysis dir from line:\n\n$line\n\n", RESET;
    exit();
  }

  my $results_dir = "$analysis_dir"."results/";
  my $lib_dir = "$results_dir"."$lib_id/";
  my $lib_fc_dir = "$lib_dir"."$flowcell_lane_read/";
  my $tmp_dir = "$lib_fc_dir"."temp/";
  my $jfb_file = "$lib_fc_dir"."junction.final.bed";
  my $jnb_file = "$lib_fc_dir"."junction.nonCanonical.bed";
  my $hmm_log_file = "$analysis_dir"."logs/$flowcell_lane_read".".log.txt";

  my $log_file_presence;
  my $result_files_presence;
  my $log_file_status = "unknown";

  print BLUE, "\n\tlib:$lib_id  lane:$flowcell_lane_read  ", RESET;
  $total++;

  #Is the log file present?
  if (-e $hmm_log_file){
    $log_file_presence = 1;
  }else{
    $log_file_presence = 0;
  }

  #Are the results files present
  if (-e $jfb_file && -e $jnb_file){
    $result_files_presence = 1;
  }else{
    $result_files_presence = 0;
  }

  #If the log file was present, what status did it report (success, failure, or neither)
  if ($log_file_presence){
    open(LOG, "$hmm_log_file") || die "\n\nCould not open log file: $hmm_log_file\n\n";
    my $log_pass = 0;
    my $log_fail = 0;
    while(<LOG>){
      chomp($_);
      if ($_ =~ /COMPLETE/){
        $log_pass = 1;
      }
      if ($_ =~ /FAILED/){
        $log_fail = 1;
      }
    }
    close(LOG);
    if ($log_fail){
      $log_file_status = "failure";
    }elsif ($log_pass){
      $log_file_status = "success";
    }else{
      $log_file_status = "neither";
    }
  }else{
    $log_file_status = "missing";
  }

  #If there is NO log file at all -> FAIL
  #If result files are present and log indicates success -> PASS
  #If result files are present and log indicated failure -> FAIL
  #If result files are present and log indicates neither -> PENDING
  #If result files are NOT present and log indicates success -> FAIL
  #If result files are NOT present and log indicates failure -> FAIL
  #If result files are NOT present and log indicates neither -> PENDING
  my $grand_status;
  if ($log_file_presence){
    if ($result_files_presence && $log_file_status eq "success"){
      $grand_status = "PASS";
    }elsif($result_files_presence && $log_file_status eq "failure"){
      $grand_status = "FAIL";
    }elsif($result_files_presence && $log_file_status eq "neither"){
      $grand_status = "PENDING";
    }elsif(!$result_files_presence && $log_file_status eq "success"){
      $grand_status = "FAIL";    
    }elsif(!$result_files_presence && $log_file_status eq "failure"){
      $grand_status = "FAIL";      
    }elsif(!$result_files_presence && $log_file_status eq "neither"){
      $grand_status = "PENDING";      
    }else{
      print RED, "\n\nCondition missing\n\n", RESET;
      exit();
    }
  }else{
    $grand_status = "FAIL";
  }

  print BLUE, "log:$log_file_status  result_files:$result_files_presence  ", RESET;

  #Summarize grand status
  if ($grand_status eq "PASS"){
    $pass_count++;
    print BLUE, "$grand_status", RESET;
  }elsif($grand_status eq "FAIL"){
    $fail_count++;
    print YELLOW, "$grand_status", RESET;
    print OUT "$line\n";
  }elsif($grand_status eq "PENDING"){
    $pending_count++;
    print CYAN, "$grand_status", RESET;
    if (defined($clean)){
      if ($clean =~ /yes|y/i){
        print OUT "$line\n";
      }
    }
  }else{
    print RED, "\n\nUndefined status\n\n", RESET;
    exit();
  }
  
}
close(IN);
close(OUT);

my $pass_p = sprintf("%.1f", (($pass_count/$total)*100));
my $fail_p = sprintf("%.1f", (($fail_count/$total)*100));
my $pending_p = sprintf("%.1f", (($pending_count/$total)*100));

print BLUE, "\n\nFound $pass_count ($pass_p%) completed jobs and $fail_count ($fail_p%) failed jobs and $pending_count ($pending_p%) pending jobs of $total total jobs\n\n", RESET;

exit();


