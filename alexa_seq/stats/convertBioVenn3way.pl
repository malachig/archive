#!/usr/bin/perl -w
#Written by Malachi Griffith and Obi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#The purpose of this script is to reduce the size of datasets to be imported to BioVenn: http://www.cmbi.ru.nl/biovenn/
#BioVenn takes as input up to three lists of IDs, it then determines the two-way or three-way overlap between these IDs and plots a Venn Diagram
#Unfortunately it does not accept datasets with more than 100,000 entries
#The purpose of this script is to take two/three datafiles containing IDs to be analyzed and reduce the size of the data while ensuring that the overlap relationships are maintained

#1.) Get IDs from file 'x' and file 'y'.  Assume these IDs are the first non-whitespace field encountered
#2.) Make a hash of all IDs.  eg. %h{$id}{x} = 1, %h{$id}{y} = 1, %h{$id}{z} = 1
#3.) Get the counts for each overlap class: 'x', 'y', 'z', 'x+y', 'y+z', 'x+z', and 'x+y+z'
#4.) Get the proportions for each overlap class: 'x', 'y', 'z', 'x+y', 'y+z', 'x+z', and 'x+y+z'
#5.) Using these proportions, create a dummy dataset with a total of N distinct IDs (where N is less than 100,000).  Maintain the proportions of each overlap class identified above
#6.) Create a subdirectory in the working directory for results files
#7.) Within this directory, write the dummy data to a new 'x', 'y' and 'z' file.
#8.) Also write a summary file to store the source files, names, and original and adjusted counts for: 'x', 'y', 'z', 'x+y', 'y+z', 'x+z', and 'x+y+z'

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);

#Load the ALEXA modules
my $script_dir;
BEGIN {
  use Cwd 'abs_path';
  if (abs_path($0) =~ /(.*)\/.*\/.*\.pl/){
    $script_dir = $1;
    push (@INC, $1);
  }
}
use utilities::utility qw(:all);

my $file_x = '';
my $file_y = '';
my $file_z = '';
my $name_x = '';
my $name_y = '';
my $name_z = '';
my $working_dir = '';
my $comparison_name = '';
my $n = '';

GetOptions('file_x=s'=>\$file_x, 'file_y=s'=>\$file_y, 'file_z=s'=>\$file_z,
           'name_x=s'=>\$name_x, 'name_y=s'=>\$name_y, 'name_z=s'=>\$name_z,
           'working_dir=s'=>\$working_dir, 'comparison_name=s'=>\$comparison_name, 'n=i'=>\$n);

print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tThis script formats data for input to BioVenn: http://www.cmbi.ru.nl/cdd/biovenn/", RESET;
print GREEN, "\n\tBasically it analyzes your IDs and creates a smaller set of IDs with the same overlap characteristics to convenient input to BioVenn", RESET;
print GREEN, "\n\tThis script assumes that IDs to be analyzed are the first non white-space elements in each row of each file", RESET;
print GREEN, "\n\tSpecify the input files containing the IDs to be formated using:  --file_x  --file_y  and  --file_z", RESET;
print GREEN, "\n\tSpecify a name corresponding to each of these using:  --name_x  --name_y  and  --name_z", RESET;
print GREEN, "\n\tSpecify a working directory where results will be stored using: --working_dir", RESET;
print GREEN, "\n\tSpecify a name for the comparison using:  --comparison_name", RESET;
print GREEN, "\n\tSpecify the total distinct IDs that will be allowed using:  --n", RESET;
print GREEN, "\n\nExample: convertBioVenn2way.pl  --file_x=alexa_seq_junctions_list.txt  --file_y=cufflink_junctions_list.txt  --file_z=scripture_junctions_list.txt  --name_x=ALEXA-Seq  --name_y=Cufflinks  --name_z=Scripture  --working_dir=/projects/alexa2/alexa_seq/read_records/HS04391/Junctions_v53/Summary/temp/  --comparison_name=ALEXA-Seq_vs_Cufflinks_vs_Scripture  --n=1000\n\n", RESET;

#Check user supplied options
unless ($file_x && $file_y && $file_z && $name_x && $name_y && $name_z && $working_dir && $comparison_name && $n){
  print RED, "\n\nRequired parameter missing!\n\n", RESET;
  exit();
}
$working_dir = &checkDir('-dir'=>$working_dir, '-clear'=>"no");

unless (-e $file_x && -e $file_y && -e $file_z){
  print RED, "\n\nInput files not found\n\n", RESET;
  exit();
}

#1.) Get IDs from file 'x' and file 'y'.  Assume these IDs are the first non-whitespace field encountered
#2.) Make a hash of all IDs.  eg. %h{$id}{x} = 1, %h{$id}{y} = 1, %h{$id}{z} = 1
open (X, "$file_x") || die "\n\nCould not open input file: $file_x\n\n";
open (Y, "$file_y") || die "\n\nCould not open input file: $file_y\n\n";
open (Z, "$file_z") || die "\n\nCould not open input file: $file_z\n\n";

my %h;
while(<X>){
  chomp($_);
  if ($_ =~ /^(\S+)/){
    $h{$_}{x}=1;
  }
}
while(<Y>){
  chomp($_);
  if ($_ =~ /^(\S+)/){
    $h{$_}{y}=1;
  }
}
while(<Z>){
  chomp($_);
  if ($_ =~ /^(\S+)/){
    $h{$_}{z}=1;
  }
}
close(X);
close(Y);
close(Z);


#3.) Get the counts for each overlap class: 'x', 'y', 'z', 'x+y', 'y+z', 'x+z', and 'x+y+z'
my $x_count = 0;
my $y_count = 0;
my $z_count = 0;
my $xy_count = 0;
my $xz_count = 0;
my $yz_count = 0;
my $xyz_count = 0;

my $grand_count = 0;
while (my ($id) = each %h){
  if($h{$id}{x} && $h{$id}{y} && $h{$id}{z}){
    $xyz_count++; $grand_count++;
  }elsif ($h{$id}{x} && $h{$id}{y}){
    $xy_count++; $grand_count++;
  }elsif ($h{$id}{x} && $h{$id}{z}){
    $xz_count++; $grand_count++;
  }elsif ($h{$id}{y} && $h{$id}{z}){
    $yz_count++; $grand_count++;
  }elsif($h{$id}{x}){
    $x_count++; $grand_count++;
  }elsif($h{$id}{y}){
    $y_count++; $grand_count++;
  }elsif($h{$id}{z}){
    $z_count++; $grand_count++;
  }
}

#4.) Get the proportions for each overlap class: 'x', 'y', 'z', 'x+y', 'y+z', 'x+z', and 'x+y+z'
my $x_p = $x_count/$grand_count;
my $y_p = $y_count/$grand_count;
my $z_p = $z_count/$grand_count;
my $xy_p = $xy_count/$grand_count;
my $xz_p = $xz_count/$grand_count;
my $yz_p = $yz_count/$grand_count;
my $xyz_p = $xyz_count/$grand_count;

#5.) Using these proportions, create a dummy dataset with a total of N distinct IDs (where N is less than 100,000).  Maintain the proportions of each overlap class identified above
my $x_adj_count = sprintf("%.0f", ($n*$x_p));
my $y_adj_count = sprintf("%.0f", ($n*$y_p));
my $z_adj_count = sprintf("%.0f", ($n*$z_p));
my $xy_adj_count = sprintf("%.0f", ($n*$xy_p));
my $xz_adj_count = sprintf("%.0f", ($n*$xz_p));
my $yz_adj_count = sprintf("%.0f", ($n*$yz_p));
my $xyz_adj_count = sprintf("%.0f", ($n*$xyz_p));
my $all_adj_count = $x_adj_count+$y_adj_count+$z_adj_count+$xy_adj_count+$xz_adj_count+$yz_adj_count+$xyz_adj_count;

#6.) Create a subdirectory in the working directory for results files
my $subdir = "$working_dir"."$comparison_name"."_BioVenn/";
unless (-e $subdir && -d $subdir){
  mkdir($subdir);
}

#7.) Within this directory, write the dummy data to a new 'x', 'y', and 'z' file.
my $new_x_file = "$subdir"."$name_x".".txt";
my $new_y_file = "$subdir"."$name_y".".txt";
my $new_z_file = "$subdir"."$name_z".".txt";
my $dummy_id = 0;
open (X2, ">$new_x_file") || die "\n\nCould not open output file: $new_x_file\n\n";
open (Y2, ">$new_y_file") || die "\n\nCould not open output file: $new_y_file\n\n";
open (Z2, ">$new_z_file") || die "\n\nCould not open output file: $new_z_file\n\n";

for (my $i = 0; $i < $xyz_adj_count; $i++){
  $dummy_id++;
  print X2 "$dummy_id\n";
  print Y2 "$dummy_id\n";
  print Z2 "$dummy_id\n";
}
for (my $i = 0; $i < $xy_adj_count; $i++){
  $dummy_id++;
  print X2 "$dummy_id\n";
  print Y2 "$dummy_id\n";
}
for (my $i = 0; $i < $xz_adj_count; $i++){
  $dummy_id++;
  print X2 "$dummy_id\n";
  print Z2 "$dummy_id\n";
}
for (my $i = 0; $i < $yz_adj_count; $i++){
  $dummy_id++;
  print Y2 "$dummy_id\n";
  print Z2 "$dummy_id\n";
}
for (my $i = 0; $i < $x_adj_count; $i++){
  $dummy_id++;
  print X2 "$dummy_id\n";
}
for (my $i = 0; $i < $y_adj_count; $i++){
  $dummy_id++;
  print Y2 "$dummy_id\n";
}
for (my $i = 0; $i < $z_adj_count; $i++){
  $dummy_id++;
  print Z2 "$dummy_id\n";
}
close(X2);
close(Y2);
close(Z2);


#8.) Also write a summary file to store the source files, names, and original and adjusted counts for 'x', 'y', 'z', 'x+y', 'y+z', 'x+z', and 'x+y+z'
my $summary_file = "$subdir"."Summary.txt";
my $x_pf = sprintf("%.3f", $x_p);
my $y_pf = sprintf("%.3f", $y_p);
my $z_pf = sprintf("%.3f", $z_p);
my $xy_pf = sprintf("%.3f", $xy_p);
my $xz_pf = sprintf("%.3f", $xz_p);
my $yz_pf = sprintf("%.3f", $yz_p);
my $xyz_pf = sprintf("%.3f", $xyz_p);
my $all_pf = $x_pf+$y_pf+$z_pf+$xy_pf+$xz_pf+$yz_pf+$xyz_pf;

open (SUMMARY, ">$summary_file") || die "\n\nCould not open output file: $summary_file\n\n";
print SUMMARY "Input data\n";
print SUMMARY "X = $name_x ($file_x)\n";
print SUMMARY "Y = $name_y ($file_y)\n";
print SUMMARY "Z = $name_z ($file_z)\n";
print SUMMARY "\nOriginal counts\n";
print SUMMARY "All (distinct ids) = $grand_count\n";
print SUMMARY "X = $x_count ($name_x only)\n";
print SUMMARY "Y = $y_count ($name_y only)\n";
print SUMMARY "Z = $z_count ($name_z only)\n";
print SUMMARY "XY = $xy_count ($name_x + $name_y)\n";
print SUMMARY "XZ = $xz_count ($name_x + $name_z)\n";
print SUMMARY "YZ = $yz_count ($name_y + $name_z)\n";
print SUMMARY "XYZ = $xyz_count ($name_x + $name_y + name_z)\n";
print SUMMARY "\nProportions (rounded)\n";
print SUMMARY "All = $all_pf\n";
print SUMMARY "X = $x_pf ($name_x only)\n";
print SUMMARY "Y = $y_pf ($name_y only)\n";
print SUMMARY "Z = $z_pf ($name_z only)\n";
print SUMMARY "XY = $xy_pf ($name_x + $name_y)\n";
print SUMMARY "XZ = $xz_pf ($name_x + $name_z)\n";
print SUMMARY "YZ = $yz_pf ($name_y + $name_z)\n";
print SUMMARY "XYZ = $xyz_pf ($name_x + $name_y + $name_z)\n";

print SUMMARY "\nAdjusted counts (for N = $n)\n";
print SUMMARY "All (distinct ids) = $all_adj_count\n";
print SUMMARY "X = $x_adj_count ($name_x only)\n";
print SUMMARY "Y = $y_adj_count ($name_y only)\n";
print SUMMARY "Z = $z_adj_count ($name_z only)\n";
print SUMMARY "XY = $xy_adj_count ($name_x + $name_y)\n";
print SUMMARY "XZ = $xz_adj_count ($name_x + $name_z)\n";
print SUMMARY "YZ = $yz_adj_count ($name_y + $name_z)\n";
print SUMMARY "XYZ = $xyz_adj_count ($name_x + $name_y + $name_z)\n";
close(SUMMARY);

print BLUE, "\n\nInput data\n", RESET;
print YELLOW, "X = $name_x ($file_x)\n", RESET;
print YELLOW, "Y = $name_y ($file_y)\n", RESET;
print YELLOW, "Z = $name_z ($file_z)\n", RESET;
print BLUE, "\nOriginal counts\n", RESET;
print YELLOW, "All (distinct ids) = $grand_count\n", RESET;
print YELLOW, "X = $x_count ($name_x only)\n", RESET;
print YELLOW, "Y = $y_count ($name_y only)\n", RESET;
print YELLOW, "Z = $z_count ($name_z only)\n", RESET;
print YELLOW, "XY = $xy_count ($name_x + $name_y)\n", RESET;
print YELLOW, "XZ = $xz_count ($name_x + $name_z)\n", RESET;
print YELLOW, "YZ = $yz_count ($name_y + $name_z)\n", RESET;
print YELLOW, "XYZ = $xyz_count ($name_x + $name_y + $name_z)\n", RESET;
print BLUE, "\nProportions (rounded)\n", RESET;
print YELLOW, "All = $all_pf\n", RESET;
print YELLOW, "X = $x_pf ($name_x only)\n", RESET;
print YELLOW, "Y = $y_pf ($name_y only)\n", RESET;
print YELLOW, "Z = $z_pf ($name_z only)\n", RESET;
print YELLOW, "XY = $xy_pf ($name_x + $name_y)\n", RESET;
print YELLOW, "XZ = $xz_pf ($name_x + $name_z)\n", RESET;
print YELLOW, "YZ = $yz_pf ($name_y + $name_z)\n", RESET;
print YELLOW, "XYZ = $xyz_pf ($name_x + $name_y + $name_z)\n", RESET;
print BLUE, "\nAdjusted counts (for N = $n)\n", RESET;
print YELLOW, "All (distinct ids) = $all_adj_count\n", RESET;
print YELLOW, "X = $x_adj_count ($name_x only)\n", RESET;
print YELLOW, "Y = $y_adj_count ($name_y only)\n", RESET;
print YELLOW, "Z = $z_adj_count ($name_z only)\n", RESET;
print YELLOW, "XY = $xy_adj_count ($name_x + $name_y)\n", RESET;
print YELLOW, "XZ = $xz_adj_count ($name_x + $name_z)\n", RESET;
print YELLOW, "YZ = $yz_adj_count ($name_y + $name_z)\n", RESET;
print YELLOW, "XYZ = $xyz_adj_count ($name_x + $name_y + $name_z)\n", RESET;
print BLUE, "\n\nResults stored here:\n$subdir\n\n", RESET;

exit();

