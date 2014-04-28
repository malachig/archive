#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#The purpose of this script is to modify settings within existing custom UCSC tracks

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);

#Load the ALEXA modules
BEGIN {
  use Cwd 'abs_path';
  if (abs_path($0) =~ /(.*)\/.*\/.*\.pl/){
    push (@INC, $1);
  }
}
use utilities::utility qw(:all);

my $ucsc_dir = '';
my $new_dir = '';

GetOptions('ucsc_dir=s'=>\$ucsc_dir, 'new_dir=s'=>\$new_dir);

print GREEN, "\n\nUsage:  updateCustomTracks.pl  --ucsc_dir=/home/malachig/www/public/htdocs/solexa/DE/MM0472/  --new_dir=/home/malachig/www/public/htdocs/solexa/DE/MM0472/temp/ ", RESET;

unless ($ucsc_dir && $new_dir){
  print RED, "\n\nRequired parameter missing!\n\n", RESET;
  exit();
}

#Check input files
$ucsc_dir = &checkDir('-dir'=>$ucsc_dir, '-clear'=>"no");
$new_dir = &checkDir('-dir'=>$new_dir, '-clear'=>"no");

#Get all custom ucsc tracks
opendir(DIRHANDLE,"$ucsc_dir") || die "\nCannot open directory: $ucsc_dir\n\n";
my @files = readdir(DIRHANDLE);
closedir(DIRHANDLE);
my %files;
my $c = 0;
foreach my $file(@files){
  my $file_path = "$ucsc_dir"."$file";
  if ($file =~ /(.*\.txt)\.gz/){
    $c++;
    $files{$c}{old_file} = $file_path;
    $files{$c}{base} = $1;
    my $new_file_path = "$new_dir"."$1";
    $files{$c}{new_file} = $new_file_path;
  }
}
my $fc = keys %files;
print BLUE, "\n\nFound $fc files.  Writing new versions to: $new_dir", RESET;

#Parse through each file and write a new file in the new dir - watch for track and browser lines and update setting as you go
foreach my $c (sort {$a <=> $b} keys %files){
  print BLUE, "\n$c: Processing $files{$c}{base}", RESET;

  open (OLD, "zcat $files{$c}{old_file} |") || die "\nCould not open old file: $files{$c}{old_file}\n\n";
  open (NEW, ">$files{$c}{new_file}") || die "\nCould not open new file: $files{$c}{new_file} \n\n";
  while(<OLD>){
    if ($_ =~ /browser\shide\sall/){
      print NEW "#browser hide all\n";
    }elsif($_ =~ /(.*)(maxHeightPixels\=\S+)(.*)/){
      #maxHeightPixels=120:80:10
      print NEW "$1"."maxHeightPixels=120:40:10"."$3\n";
    }else{
      print NEW "$_";
    }
  }
  close(OLD);
  close(NEW);

  #Compress the new file:
  my $cmd = "gzip $files{$c}{new_file}";
  system($cmd);

}



exit();
