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

my $input_dir = '';
my $output_dir = '';
my $library_name = '';
my $color_set = '';
my $global_priority = '';
my $force = '';

GetOptions ('input_dir=s'=>\$input_dir, 'output_dir=s'=>\$output_dir, 'library_name=s'=>\$library_name, 'color_set=i'=>\$color_set, 'priority=i'=>\$global_priority, 'force=s'=>\$force);


#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tThis script merges UCSC files in a directory and places the merged versions in a new directory", RESET;
print GREEN, "\n\tSpecify the input directory using: --input_dir", RESET;
print GREEN, "\n\tSpecify the output directory using: --output_dir", RESET;
print GREEN, "\n\tSpecify the name of the library using:  --library_name", RESET;
print GREEN, "\n\tSpecify a color set using: --color_set=1 or --color_set=2", RESET;
print GREEN, "\n\tSpecify a priority value to control display order of expression wig tracks using:  --priority (use 80-100; smaller numbers displayed first)", RESET; 
print GREEN, "\n\tTo force cleaning of output directories use:  --force=yes", RESET;
print GREEN, "\n\nExample: mergeUcscTrackFiles.pl  --input_dir=/projects/malachig/alexa_seq/temp/website/5FU/ucsc/HS04391/  --output_dir=/projects/malachig/alexa_seq/temp/website/5FU/ucsc/HS04391/  --library_name='MIP101'  --color_set=1  --priority=80\n\n", RESET;

unless ($input_dir && $output_dir && $library_name && $color_set && $global_priority){
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  exit();
}
if ($input_dir eq $output_dir){
  print RED, "\ninput and output directories must be different!\n\n", RESET;
  exit();
}

my $color;
if ($color_set == 1){
  $color = "102,51,204";  #Dark purple
}else{
  $color = "153,102,204"; #Light purple
}

$input_dir = &checkDir('-dir'=>$input_dir, '-clear'=>"no");
if ($force){
  $output_dir = &checkDir('-dir'=>$output_dir, '-clear'=>"yes", '-force'=>$force);
}else{
  $output_dir = &checkDir('-dir'=>$output_dir, '-clear'=>"yes");
}


my $files_ref = &getFiles('-dir'=>$input_dir);

#Go through file of each region set
foreach my $region (sort {$a cmp $b} keys %{$files_ref}){
  my @file_list = @{$files_ref->{$region}->{file_list}};
  my %tracks;

  print MAGENTA, "\n\nProcessing: $region\n", RESET;

  foreach my $file (@file_list){

    open (UCSC, "zcat $file |") || die "\nCould not open input file: $file\n\n";

    my $new_track = 0;
    my $current_track_name = 0;
    while(<UCSC>){
      chomp($_);
      #Skip comment lines
      if ($_ =~ /^\#/){
        next();
      }
      #Strip out non-valid lines - assumes a combination of gff and fixedStep wiggle tracks only
      unless ($_ =~ /^track|^fixedStep|chr|^\d+/){
        next();
      }

      #Watch out for track def lines
      if ($_ =~ /^track\s+name\=(\S+)\s+.*\s+priority\=(\d+)/){
        my $track_name = $1;
        my $old_track_name = $1;
        my $priority = $2;
        my $library = "Lib";
        if ($track_name =~ /\_RAW\_/){
          #print YELLOW, "\nRAW: $_", RESET;
          if ($track_name =~ /^(.*)\_/){
            $library = $1;
          }

          $track_name = "$library"."_RAW";
          $tracks{$track_name}{type} = "wiggle";
          $tracks{$track_name}{priority} = 79;
          $tracks{$track_name}{track_def_line} = "track name=$track_name description=\"RAW Coverage. ($library_name)\" type=wiggle_0 color=$color yLineMark=0.0 yLineOnOff=on visibility=hide autoScale=on graphType=bar smoothingWindow=off maxHeightPixels=120:80:10 priority=79";

        }elsif($track_name =~ /\_N1\_/){
          #print YELLOW, "\nN1: $_", RESET;
          if ($track_name =~ /^(.*)\_N1\_/){
            $library = $1;
          }
          $track_name = "$library"."_N1";
          $tracks{$track_name}{type} = "wiggle";
          $tracks{$track_name}{priority} = $global_priority;
          $tracks{$track_name}{track_def_line} = "track name=$track_name description=\"Normalized Coverage. ($library_name)\" type=wiggle_0 color=$color yLineMark=0.0 yLineOnOff=on visibility=full autoScale=on graphType=bar smoothingWindow=off maxHeightPixels=120:80:10 priority=$global_priority";

        }else{
          $tracks{$track_name}{type} = "gff";
          $tracks{$track_name}{priority} = $priority;
          $tracks{$track_name}{track_def_line} = $_;
        }
        $new_track = 1;
        $current_track_name = $track_name;
        print YELLOW, "\n\t\tFound track: $old_track_name (will be stored as: $track_name in merged file at priority = $tracks{$track_name}{priority}", RESET;
        next();
      }

      #Unless the first track of this file has started, skip
      unless ($current_track_name){
        next();
      }

      if ($tracks{$current_track_name}{records}){
        push(@{$tracks{$current_track_name}{records}}, $_);
      }else{
        my @tmp;
        push(@tmp, $_);
        $tracks{$current_track_name}{records} = \@tmp;
      }
    }
    close(UCSC);
  }

  #Go through the tracks stored and print out a new merged ucsc file 
  my $new_ucsc_file = "$output_dir"."$region"."_merged.txt";
  print BLUE, "\n\n\tPrinting new merged file: $new_ucsc_file\n", RESET;
  open(OUT, ">$new_ucsc_file") || die "\nCould not open new UCSC file: $new_ucsc_file\n\n";

  #Print the browser line
  print OUT "browser hide all\nbrowser full knownGene\nbrowser pack multiz28way\n\n";
  foreach my $track_name (sort {$tracks{$a}->{priority} <=> $tracks{$b}->{priority}} keys %tracks){

    my $records_ref = $tracks{$track_name}{records};
    my $count = 0;

    #make sure there are some records to be printed
    if ($records_ref){
      $count = scalar(@{$records_ref});
      print BLUE, "\n\t\tPrinting $count records for $track_name to the merge file (priority = $tracks{$track_name}{priority})", RESET;
      print OUT "$tracks{$track_name}{track_def_line}\n";
      foreach my $record (@{$records_ref}){
        print OUT "$record\n";
      }
    }
   }
  close(OUT);

  my $gzip_cmd = "gzip $new_ucsc_file";
  system($gzip_cmd);
}

print "\n\nSCRIPT COMPLETE\n\n";

exit();


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

    $file_count++;
    #print BLUE, "\n\t$file_path was added to the list of files to be processed", RESET;
    #$mapped_read_files{$file_count}{path} = $file_path;

    #Get the region name
    if ($test_file =~ /(chr[\d|X|Y|MT]+\_\d+)/){
      if ($files{$1}){
        push(@{$files{$1}{file_list}}, $file_path);
      }else{
        my @tmp;
        push(@tmp, $file_path);
        $files{$1}{file_list} = \@tmp;
      }
    }
  }
  print BLUE, "", RESET;
  foreach my $region (sort {$a cmp $b} keys %files){
    my $fc = scalar (@{$files{$region}{file_list}});
    print YELLOW, "\n\tFound: $region with $fc files to be merged", RESET;
  }
  my $total = keys %files;
  print YELLOW, "\n\n\tFound a total of $total regions with files to be to be merged (UCSC recognized chromosomes only)\n", RESET;

  return(\%files);
}


