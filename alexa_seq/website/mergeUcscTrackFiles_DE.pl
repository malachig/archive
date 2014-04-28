#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#This script takes custom track files that were created to display different data for the same genomic regions and merges them together
#It assumes that files to be merged are named with the format: chr2_2_NNNNNNN.txt  (i.e. chr$chromosome_$region_number_datatype.txt)

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;
use BerkeleyDB;

#Load the ALEXA modules
BEGIN {
  use Cwd 'abs_path';
  if (abs_path($0) =~ /(.*)\/.*\/.*\.pl/){
    push (@INC, $1);
  }
}
use utilities::utility qw(:all);

my $input_dir_A = '';
my $input_dir_B = '';
my $output_dir_A = '';
my $output_dir_B = '';
my $comparison_name = '';
my $short_name = '';
my $view_limit = '';
my $force = '';

GetOptions ('input_dir_A=s'=>\$input_dir_A, 'input_dir_B=s'=>\$input_dir_B, 'output_dir_A=s'=>\$output_dir_A, 'output_dir_B=s'=>\$output_dir_B, 
            'comparison_name=s'=>\$comparison_name, 'short_name=s'=>\$short_name, 'view_limit=s'=>\$view_limit, 'force=s'=>\$force);

#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tThis script merges UCSC files in a directory and places the merged versions in a new directory", RESET;
print GREEN, "\n\tSpecify the input directory for library A using: --input_dir_A", RESET;
print GREEN, "\n\tSpecify the input directory for library B using: --input_dir_B", RESET;
print GREEN, "\n\tSpecify the output directory using: --output_dir", RESET;
print GREEN, "\n\tSpecify the descriptive name of the comparison being merged using:  --comparison_name", RESET;
print GREEN, "\n\tSpecify a short name to be used in naming the UCSC track (no spaces) using:  --short_name", RESET;
print GREEN, "\n\tSpecify the desired data range to display in the UCSC DE track using:  --view_limit='lower:upper'", RESET
print GREEN, "\n\tTo force cleaning of output directories use:  --force=yes", RESET;
print GREEN, "\n\nExample: mergeUcscTrackFiles.pl  --input_dir_A=/projects/malachig/alexa_seq/temp/website/5FU/ucsc/HS04401/combined/  --input_dir_B=/projects/malachig/alexa_seq/temp/website/5FU/ucsc/HS04391/combined/  --output_dir_A=/projects/malachig/alexa_seq/temp/website/5FU/ucsc/DE/HS04401/HS04401_vs_HS04391/  --output_dir_B=/projects/malachig/alexa_seq/temp/website/5FU/ucsc/DE/HS04391/HS04401_vs_HS04391/  --comparison_name='MIP/5FU - MIP101'  --short_name='MIP_DE'  --view_limit='-2.5:2.5'\n\n", RESET;

unless ($input_dir_A && $input_dir_B && $output_dir_A && $output_dir_B && $comparison_name && $short_name && $view_limit){
  print RED, "\n\nRequired input parameter(s) missing\n\n", RESET;
  exit();
}
if ($input_dir_A eq $input_dir_B || $input_dir_A eq $output_dir_A || $input_dir_B eq $output_dir_B || $output_dir_A eq $output_dir_B){
  print RED, "\n\ninput and output A and B directories must ALL be different - output dirs will be cleared before writing!\n\n", RESET;
  exit();
}

$input_dir_A = &checkDir('-dir'=>$input_dir_A, '-clear'=>"no");
$input_dir_B = &checkDir('-dir'=>$input_dir_B, '-clear'=>"no");

if ($force){
  $output_dir_A = &checkDir('-dir'=>$output_dir_A, '-clear'=>"yes", '-force'=>$force);
  $output_dir_B = &checkDir('-dir'=>$output_dir_B, '-clear'=>"yes", '-force'=>$force);
}else{
  $output_dir_A = &checkDir('-dir'=>$output_dir_A, '-clear'=>"yes");
  $output_dir_B = &checkDir('-dir'=>$output_dir_B, '-clear'=>"yes");
}

my $files_ref_A = &getFiles('-dir'=>$input_dir_A);
my $files_ref_B = &getFiles('-dir'=>$input_dir_B);

foreach my $region (sort {$a cmp $b} keys %{$files_ref_A}){
  my $chromosome;
  if ($region =~ /(chr.*)\_\d+/){
    $chromosome = $1;
  }else{
    print RED, "\n\nCould not determine chromosome from chr_region: $region\n\n", RESET;
    exit();
  }

  print MAGENTA, "\n\nProcessing: $region\n", RESET;

  my $file_a = $files_ref_A->{$region}->{path};
  my $file_b = $files_ref_B->{$region}->{path};

  #Go through each file
  if ($file_a && $file_b){
    print BLUE, "\nCalculating DE values for:\n\t$file_a\n\t$file_b\n", RESET;
  }else{
    print YELLOW, "\nCould not find a matching file in library B for $file_a\n\n", RESET;
    next();
  }

  #The following temp files will be used to store the normalized wig track data for each library
  #These will be used to join the track data of library A and B
  #The result will be a file with: norm expression track for LibA, norm expression track for LibB, DE track
  my $temp_A = "$output_dir_A". "$region"."_DE_A_TEMP.txt";
  my $temp_B = "$output_dir_B". "$region"."_DE_B_TEMP.txt";

  #Get data for file A
  my $data_A_ref = &getTrackData('-file'=>$file_a, '-temp_file'=>$temp_A);

  #Get data for file B
  my $data_B_ref = &getTrackData('-file'=>$file_b, '-temp_file'=>$temp_B);

  #Note that in the UCSC wiggle format, data is only stored for those positions with some coverage observed
  #So the two data structures created above do not correspond to the same coordinates exactly (although substantial overlap is expected if the libraries have decent depth)
  #Need to go through each data structure and initialize the positions that were not observed in both libraries.
  print BLUE, "\n\nInitializing positions for which coverage was observed in one but not both libraries", RESET;
  foreach my $pos (keys %{$data_A_ref}){
    unless (defined($data_B_ref->{$pos})){
      $data_B_ref->{$pos} = 0;
    }
  }
  foreach my $pos (keys %{$data_B_ref}){
    unless (defined($data_A_ref->{$pos})){
      $data_A_ref->{$pos} = 0;
    }
  }
  my $data_a_size = keys %{$data_A_ref};
  my $data_b_size = keys %{$data_B_ref};
  print BLUE, "\n\tUpdated size of data objects: $data_a_size (library A), $data_b_size (library B)", RESET;

  my $tempfile = "$output_dir_A". "$region"."_DE_TEMP.txt";

  print BLUE, "\n\nPrinting temp UCSC track file containing the differential expression track: $tempfile", RESET;
  open (UCSC, ">$tempfile") || die "\nCould not open output file: $tempfile";

  #Print the browser def line
  #print UCSC "browser hide all\nbrowser full knownGene\nbrowser pack multiz28way\n\n";

  #Print the track def line
  print UCSC "track name=\"$short_name\" description=\"Log2 Differential Expression. ($comparison_name)\" type=wiggle_0 color=51,204,153 yLineMark=0.0 yLineOnOff=on visibility=full autoScale=off viewLimits=$view_limit graphType=bar smoothingWindow=off maxHeightPixels=120:80:10 priority=99\n"; 
  my $last_pos = -1000;
  foreach my $pos (sort {$a <=> $b} keys %{$data_A_ref}){

    #Add 16 to all values to prevent log of 0, and to stabilize variance
    my $dataA = ($data_A_ref->{$pos})+16;
    my $dataB = ($data_B_ref->{$pos})+16;
    my $de = ((log($dataA))/(log(2))) - ((log($dataB))/(log(2)));
    #print YELLOW, "\n\t$data2\t$data1\t$de", RESET;

    #If this is the start of a new block of contiguous positions, write out a fixedStep def line - otherwise simply write the latest data value
    if ($pos == $last_pos+1){
      print UCSC "$de\n";
    }else{
      if ($chromosome eq "chrMT"){$chromosome = "chrM";}
      print UCSC "fixedStep chrom=$chromosome start=$pos step=1\n";               
      print UCSC "$de\n";
    }
    $last_pos = $pos;
  }
  close(UCSC);

  #Join the temp DE file to the end of both input files - then delete the temp file and compress the new files
  print BLUE, "\n\nJoining new UCSC track file containing the differential expression track to each input file\n\n", RESET;
  
  my $outfile_A = "$output_dir_A". "$region"."_DE_merged.txt";
  my $outfile_B = "$output_dir_B". "$region"."_DE_merged.txt";

  my $cmd1 = "zcat $file_a | cat - $temp_B $tempfile > $outfile_A";
  my $cmd2 = "zcat $file_b | cat - $temp_A $tempfile > $outfile_B";
  my $cmd3 = "rm -f $tempfile $temp_A $temp_B";
  my $cmd4 = "gzip $outfile_A";
  my $cmd5 = "gzip $outfile_B";
  system($cmd1);
  system($cmd2);
  system($cmd3);
  system($cmd4);
  system($cmd5);
}

print "\n\nSCRIPT COMPLETE\n\n";

exit();



###############################################################################################################################
#getTrackData
###############################################################################################################################
sub getTrackData{
  my %args = @_;
  my $file = $args{'-file'};
  my $temp_file = $args{'-temp_file'};

  print BLUE, "\nStoring normalized expression data for file: $file", RESET;

  my %data;

  #Open an output file to store the DE wiggle track for this file
  open(TMP, ">$temp_file") || die "\nCould not open temp output file: $temp_file\n\n";

  #Get the normalized track data for file A
  open (FILE, "zcat $file | ") || die "\nCould not open file: $file\n\n";
  my $track_start = 0;
  my $fixed_step_start = 0;
  while(<FILE>){

    chomp($_);
    #Skip comment lines
    if ($_ =~ /^\#/){
      next();
    }
    #Strip out non-valid lines - assumes a combination of gff and fixedStep wiggle tracks only
    unless ($_ =~ /^track|^fixedStep|^\d+/){
      next();
    }

    #Watch for the beginning of the normalized track data
    if ($_ =~ /^track\s+name\=(\S+)\s+.*\s+priority\=(\d+)/){
      my $track_name = $1;
      if ($track_name =~  /\_N1/){
        print TMP "$_\n";
        $track_start = 1;
        #print YELLOW, "\nFound norm track:\n$_\n", RESET;
      }
      next();
    }

    #If we have already been processing the track of interest and a new track starts abort
    if (($track_start == 1) && ($_ =~ /^track/)){
      last();
    }

    #Unless the beginning of the norm track has been found skip line
    unless ($track_start == 1){
      next();
    }

    #Watch for start of a fixedStep block
    if ($_ =~ /^fixedStep.*start\=(\d+)/){
      #print YELLOW, "\n\tStart of fixed step: $1", RESET;

      print TMP "$_\n";

      #Initialize the fixed step start point
      $fixed_step_start = $1-1;
      next();
    }

    #At this point the line should contain data - data should be either an integer or float number
    print TMP "$_\n";

    $fixed_step_start++;
    my $value;
    if ($_ =~ /^(\d+)$/){
      $value = $1;
    }elsif($_ =~ /^(\d+\.\d+)/){
      $value = $1;
    }else{
      print RED, "\n\nFormat of data line not understood: $_\n\n", RESET;
      exit();
    }

    #If the same coordinate encountered multiple times (should not really happen), store the higher coverage value
    if (defined($data{$fixed_step_start})){
      if ($value > $data{$fixed_step_start}){
        $data{$fixed_step_start} = $value;
      }
    }else{
      $data{$fixed_step_start} = $value;
    }
  }

  my $data_stored = keys %data;
  print BLUE, "\n\tStored data for $data_stored coordinate positions", RESET;
  close(FILE);
  close(TMP);
  return(\%data);
}




########################################################################################################################################
#get input files to be merged
########################################################################################################################################
sub getFiles{
  my %args = @_;
  my $dir = $args{'-dir'};
  my %files;

  print BLUE, "\n\nSearching $dir for ucsc files to be merged", RESET;

  opendir(DIRHANDLE, "$dir") || die "\nCannot open directory: $dir\n\n";
  my @test_files = readdir(DIRHANDLE);

  my $file_count = 0;
  foreach my $test_file (sort @test_files){
    my $file_path = "$dir"."$test_file";

    #Skip directories within the specified directory
    if (-e $file_path && -d $file_path){
      print YELLOW, "\n\t$file_path  is a directory - skipping", RESET;
      next();
    }

    #If the results file is compressed uncompress it
    unless ($file_path =~ /(.*)\.gz$/){
      print RED, "\nFound an uncompressed file: $file_path\n\n\tMake sure all files are compressed before proceeding\n\n\t- A mix of compressed and uncompressed files may indicate a problem\n\n", RESET;
      exit();
    }

    if ($test_file =~ /(chr[\d|X|Y|MT]+\_\d+)/){
      $file_count++;
      $files{$1}{path} = $file_path;
    } 
  }
  print BLUE, "", RESET;
  foreach my $file (sort {$a cmp $b} keys %files){
    print YELLOW, "\n\tFound: $file\t $files{$file}{path}", RESET;
  }
  return(\%files);
}


