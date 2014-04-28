#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#The purpose of this script is to re-initialize the read records files for an entire library
#Do this if wish to redo all alignment parsing from scratch
#The status of all reads in the read records files will be reset to 'Unassigned' except for those that are assigned to 'Duplicate', 'Low_Complexity', 'Low_Quality'
#The mapping summary files for each feature type will also be deleted


use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;

#Load the ALEXA modules
BEGIN {
  use Cwd 'abs_path';
  if (abs_path($0) =~ /(.*)\/.*\/.*\.pl/){
    push (@INC, $1);
  }
}
use utilities::utility qw(:all);

my $library_id = '';
my $analysis_dir = '';
my $ensembl_version = '';

GetOptions ('library_id=s'=>\$library_id, 'analysis_dir=s'=>\$analysis_dir, 'ensembl_version=i'=>\$ensembl_version);

#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tThis script resets the read records of a library and deletes the mapping summary files so that parsing can be repeated cleanly", RESET;
print GREEN, "\n\tSpecify the library_id using: --library_id", RESET;
print GREEN, "\n\tSpecify the path to your analysis directory using: --analysis_dir", RESET;
print GREEN, "\n\tSpecify the ensembl version used for this analysis using: --ensembl_version", RESET;
print GREEN, "\n\nExample: resetReadRecords.pl  --library_id=HS0804  --analysis_dir=/projects/malachig/alexa_seq/  --ensembl_version=53\n\n", RESET;

unless ($library_id && $analysis_dir && $ensembl_version){
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  exit();
}


#Check the dir
$analysis_dir = &checkDir('-dir'=>$analysis_dir, '-clear'=>"no");

#Define the read record dir for this library
my $rr_dir = "$analysis_dir"."read_records/"."$library_id/";
$rr_dir = &checkDir('-dir'=>$rr_dir, '-clear'=>"no");

#Define the logs dir for this library
my $log_dir = "$analysis_dir"."logs/"."$library_id/";
$log_dir = &checkDir('-dir'=>$log_dir, '-clear'=>"no");


#Get the read record files
my $rr_files_ref = &getFiles('-dir'=>$rr_dir);
print BLUE, "\n\nReseting read record files:", RESET;
foreach my $c (sort {$a <=> $b} keys %{$rr_files_ref}){
  my $rr_file = $rr_files_ref->{$c}->{gz_path};
  my $rr_file_temp = $rr_files_ref->{$c}->{temp_path};
  my $rr_file_temp_gz = $rr_files_ref->{$c}->{temp_path_gz};
  print BLUE, "\n\n$c. Processing: ($rr_file -> $rr_file_temp)\n", RESET;

  &resetReadRecord('-infile'=>$rr_file, '-outfile'=>$rr_file_temp);

  #Compress the new file
  my $cmd_gz = "gzip $rr_file_temp";
  print BLUE, "\n\n$cmd_gz", RESET;
  system($cmd_gz);

  #Overwrite the old file
  my $cmd_mv = "mv -f $rr_file_temp_gz $rr_file";
  print BLUE, "\n\n$cmd_mv";
  system($cmd_mv);
}

#Clean up all the mapping summary files
print BLUE, "\n\nCleaning up all the mapping summary files for this library ID:", RESET;
my %types;
$types{'1'}{dir} = "Repeats";
$types{'1'}{file} = "Repeats";
$types{'2'}{dir} = "ENST_v"."$ensembl_version";
$types{'2'}{file} = "ENST_v"."$ensembl_version";
$types{'3'}{dir} = "Junctions_v"."$ensembl_version";
$types{'3'}{file} = "Junctions_v"."$ensembl_version";
$types{'4'}{dir} = "Boundaries_v"."$ensembl_version";
$types{'4'}{file} = "Boundaries_v"."$ensembl_version";
$types{'5'}{dir} = "Introns_v"."$ensembl_version";
$types{'5'}{file} = "Introns";
$types{'6'}{dir} = "Intergenics_v"."$ensembl_version";
$types{'6'}{file} = "Intergenics";


foreach my $t (sort {$a <=> $b} keys %types){
  my $type_dir = $types{$t}{dir};
  my $type_file = $types{$t}{file};

  foreach my $c (sort {$a <=> $b} keys %{$rr_files_ref}){
    my $rr_file_name = $rr_files_ref->{$c}->{name};

    my $file_txt = "$rr_dir"."$type_dir"."/$rr_file_name"."_"."$type_file".".txt";
    my $file_gz = "$file_txt".".gz";
    my $cmd_rm = "rm -f $file_gz $file_txt";
    print "\n$cmd_rm";
    system($cmd_rm);

  }
}

#Clean up the parsing log files:
my $cmd_rm_logs = "rm -f $log_dir"."*/parse*";
system($cmd_rm_logs);

exit();


#######################################################################################################################################################
#Get All .txt or .txt.gz files from a directory                                                                                                       #
#######################################################################################################################################################
sub getFiles{
  my %args = @_;
  my $dir = $args{'-dir'};

  #Check the dir
  $dir = &checkDir('-dir'=>$dir, '-clear'=>"no");

  #Get files from this directory
  print BLUE, "\n\nSearching $dir for files", RESET;

  my %files;
  opendir(DIRHANDLE, "$dir") || die "\nCannot open directory: $dir\n\n";
  my @test_files = readdir(DIRHANDLE);
  my $file_count = 0;

  foreach my $test_file (sort @test_files){
    my $file_path = "$dir"."$test_file";

    #Skip directories within the specified directory
    if (-e $file_path && -d $file_path){
      #print YELLOW, "\n\t$file_path  is a directory - skipping", RESET;
      next();
    }
    unless($test_file =~ /\_Lane\d+/){
      print YELLOW, "\n\t$test_file does not have the expected name format - skipping", RESET;
      next();
    }
    #If the results file is compressed uncompress it
    if ($file_path =~ /(.*)\/(.*)\.txt\.gz$/){
      $file_count++;
      $files{$file_count}{name} = "$2";
      $files{$file_count}{gz_path} = $file_path;
      $files{$file_count}{txt_path} = "$1/$2".".txt";
      $files{$file_count}{temp_path} = "$1/$2".".tmp";
      $files{$file_count}{temp_path_gz} = "$1/$2".".tmp.gz";
    }else{
      print RED, "\nFound an uncompressed file: $file_path\n\n\tRemove before proceeding\n\n", RESET;
      exit();
    }
  }
  return(\%files);
}

#######################################################################################################################################################
#Get All .txt or .txt.gz files from a directory                                                                                                       #
#######################################################################################################################################################
sub resetReadRecord{
  my %args = @_;
  my $infile = $args{'-infile'};
  my $outfile= $args{'-outfile'};

  open (OUT, ">$outfile") || die "\n\nCout not open output file: $outfile\n\n";
  open (IN, "zcat $infile |") || die "\n\nCould not open input file: $infile\n\n";
  my $header = 1;
  my $c = 0;
  while(<IN>){
    chomp($_);
    my @line = split("\t", $_);
    if($header){
      print OUT "$_\n";
      $header = 0;
      next();
    }
    $c++;

    if ($c == 1000000){
      $| = 1; print BLUE, ".", RESET;  $| = 0;
      $c = 0;
    }

    #Reset values
    unless ($line[3] =~ /Duplicate|Low_Complexity|Low_Quality/){
      $line[3] = "Unassigned";
    }
    unless ($line[4] =~ /Duplicate|Low_Complexity|Low_Quality/){
      $line[4] = "Unassigned";
    }
    print OUT "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$line[6]\t$line[7]\t$line[8]\t$line[9]\t$line[10]\n";
  }
  close(IN);
  close(OUT);

  return();
}


