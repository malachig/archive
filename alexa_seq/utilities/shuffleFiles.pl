#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;
use Benchmark;
use Tie::File;
use List::Util 'shuffle';
use IO::File;

#Load the ALEXA libraries
BEGIN {
  use Cwd 'abs_path';
  if (abs_path($0) =~ /(.*)\/.*\/.*\.pl/){
    push (@INC, $1);
  }
}
use utilities::utility qw(:all);

my $input_dir = '';
my $output_dir = '';
my $output_file_count = '';

GetOptions ('input_dir=s'=>\$input_dir,'output_dir=s'=>\$output_dir, 'output_file_count=i'=>\$output_file_count);

#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tspecify a path to the files to be shuffled together using: --input_dir", RESET;
print GREEN, "\n\tspecify a path to the files to be shuffled together using: --output_dir", RESET;
print GREEN, "\n\tspecify the desired number of output files using: --output_file_count", RESET;
print GREEN, "\n\nExample: shuffleFiles.pl --input_dir=/projects/malachig/solexa/read_records/HS04401/ENST_v53/  --output_dir=/projects/malachig/solexa/read_records/HS04401/ENST_v53/shuffle/  --output_file_count=16\n\n", RESET;

unless ($input_dir && $output_dir && $output_file_count){
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  exit();
}

$input_dir = &checkDir('-dir'=>$input_dir, '-clear'=>"no");
$output_dir = &checkDir('-dir'=>$output_dir, '-clear'=>"no");

if ($input_dir eq $output_dir){
  print RED, "\n\nInput and output dirs can not be the same!!\n\n", RESET;
  exit();
}

#Get the input files to be processed
my %infiles;
my $master_header;
&getDataFiles('-input_dir'=>$input_dir);

#Now create N output file handles.  Each file will be read randomly and writen even across these output files
#Print the header to each file
my %outfiles;
for (my $i = 1; $i <= $output_file_count; $i++){
  my $file_name = "Temp_"."$i".".txt";
  $outfiles{$i}{name} = $file_name;
  my $path = "$output_dir"."$file_name";
  $outfiles{$i}{path} = "$path";
  my $fh = IO::File->new(">$path") || die "\nCould not open output file: $path\n\n";
  $outfiles{$i}{fh} = $fh;
  print $fh "$master_header\n";
}

#Go through each input file now and dump the output to 
my $temp_file = "$output_dir"."temp.txt";
my $current = 0;
foreach my $fc (sort {$infiles{$a}->{random} <=> $infiles{$b}->{random}} keys %infiles){
  my $file_name = $infiles{$fc}{file_name};
  my $file_path = $infiles{$fc}{file_path};
  print BLUE "\nProcessing: $file_name", RESET;  

  #Decompress file to a temp file
  my $cmd = "zcat $file_path > $temp_file";
  print YELLOW, "\n\t$cmd", RESET;
  system($cmd);

  print YELLOW, "\n\t\tRandomly drawing from file shuffling to $output_file_count files", RESET;
  tie my @file, 'Tie::File', $temp_file or die "\n\nCould not tie infile: $file_path\n\n";
  my $counter = 0;
  for (shuffle 1 .. $#file){
    $counter++;

    if($counter == 100000){
      $counter = 0;
      $| = 1; print YELLOW, ".", RESET; $| = 0;
    }
    my $file_target = &cycle('-current'=>$current, '-max'=>$output_file_count);
    $current = $file_target;

    my $fh = $outfiles{$file_target}{fh};
    print $fh "$file[$_]\n";
  }

  untie(@file);
}

#Delete the temp file
system("rm -f $temp_file");

#Close the file handles for output files
foreach my $i (keys %outfiles){
  my $fh = $outfiles{$i}{fh};
  close($fh);
}

#Compress all the temp files produced
foreach my $i (keys %outfiles){
  my $path = $outfiles{$i}{path};
  my $cmd = "gzip -f $path";
  system($cmd);
}



exit();

###########################################################################################################
#Cycle through list of numbers                                                                            #
###########################################################################################################
sub cycle{
  my %args=@_;
  my $c = $args{'-current'};
  my $max = $args{'-max'};

  $c++;
  if ($c > $max){
    $c = 1;
  }

  return($c);
}


###########################################################################################################
#Get data files and the headers                                                                           #
###########################################################################################################
sub getDataFiles{
  my %args = @_;
  my $dir = $args{'-input_dir'};

  my $dh = opendir(DIR, $dir) || die "\nCould not open directory: $dir\n\n";
  my @files = readdir(DIR);
  
  #Assign a random number to each file to allow random order processing of lane files
  srand();

  my $count = 0;
  my %headers;
  foreach my $file (@files){
    my $header = 1;
    chomp($file);
    unless ($file =~ /\.txt\.gz$/){
      next();
    }
    if (-d $file){
      next();
    }
    $count++;
    $infiles{$count}{file_name} = $file;
    $infiles{$count}{file_path} = $dir.$file;
    $infiles{$count}{random} = rand();
    #Get the header values for this file
    open (FILE, "zcat $dir$file |") || die "\nCould not open file: $dir$file";
    while(<FILE>){
      chomp($_);
      my @line = split("\t", $_);
      #Parse the column names and positions.  Check against a hard coded list of required columns before proceeding
      if ($header == 1){
        $infiles{$count}{header} = $_;
        $headers{$_}=1;
        $master_header=$_;
        my $col_count = 0;
        last();
      }
    }
    close(FILE);
  }
  closedir(DIR);

  my $files_count = keys %infiles;
  print BLUE, "\n\nFound $files_count files to be processed (all .txt.gz files in the specified directory)", RESET;

  my $headers_found = keys %headers;
  if ($headers_found > 1){
    print RED, "\n\nThese files do not all have the same header!!", RESET;
    exit();
  }

  return();
}




