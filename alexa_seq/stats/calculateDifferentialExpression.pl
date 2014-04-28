#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#This script takes a pair of input files containing expression value columns (one for each library to be compared)
#It creates a temp file containing three descriptive columns and the two data columns
#It this feeds this temp file into an R script which calculates differential expression values and assigns p-values and q-values to them

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
    push (@INC, $script_dir);
  }
}
use utilities::utility qw(:all);

my $dataname = '';            #Description of the libraries and data being used to create these results
my $fileA = '';               #File containing expression values for library A
my $fileB = '';               #File containing expression values for library B
my $columnA = '';             #Column from file1 containing data of interest 
my $columnB = '';             #Column from file2 containing data of interest
my $results_dir = '';         #Directory for output text files and figures
my $short_name = '';
my $image_type = '';
my $libraryA_name = '';
my $libraryB_name = '';
my $comp_id = '';

GetOptions('dataname=s'=>\$dataname, 'fileA=s'=>\$fileA, 'fileB=s'=>\$fileB, 'columnA=s'=>\$columnA, 'columnB=s'=>\$columnB, 'results_dir=s'=>\$results_dir, 
           'short_name=s'=>\$short_name, 'image_type=s'=>\$image_type, 'libraryA_name=s'=>\$libraryA_name, 'libraryB_name=s'=>\$libraryB_name, 'comp_id=s'=>\$comp_id);

print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tThis script generates DE results for Genes, Exons, Junctions, etc.", RESET;
print GREEN, "\n\tSpecify a descriptive name for the results using: --dataname (will be used for file names, etc.)", RESET;
print GREEN, "\n\tSpecify a comparison ID to be used in temp file naming using:  --comp_id", RESET;
print GREEN, "\n\tSpecify the complete path to a file containing expression data for 'library A' using: --fileA", RESET;
print GREEN, "\n\tSpecify the complete path to a file containing expression data for 'library B' using: --fileB", RESET;
print GREEN, "\n\tSpecify a name for library A using:  --libraryA_name", RESET;
print GREEN, "\n\tSpecify a name for library B using:  --libraryB_name", RESET;
print GREEN, "\n\tSpecify the data column to be used for 'library A' using: --columnA", RESET;
print GREEN, "\n\tSpecify the data column to be used for 'library B' using: --columnB", RESET;
print GREEN, "\n\tSpecify the complete path to an output directory for results using: --results_dir", RESET;
print GREEN, "\n\tSpecify a short name describing the data type to be used in figure legends using:  --short_name", RESET;
print GREEN, "\n\tSpecify the image type to generate using: --image_type='JPEG' or --image_type='TIFF' or --image_type='SVG'", RESET;

print GREEN, "\n\nUsage: calculateDifferentialExpression.pl  --dataname='MIP101_vs_MIP5FUR_AvgCoverage'  --comp_id=HS04391_vs_HS04401  --fileA=/projects/malachig/solexa/read_records/HS04401/ENST_v53/Summary/HS04401_GeneExpression_v53.txt  --fileB=/projects/malachig/solexa/read_records/HS04391/ENST_v53/Summary/HS04391_GeneExpression_v53.txt  --columnA='Average_Coverage_RAW'   --columnB='Average_Coverage_RAW' --results_dir=/projects/malachig/solexa/figures_and_stats/DE/5FU/ENST_v53/  --short_name='Gene'  --image_type='JPEG'  --libraryA_name='MIP/5FU'  --libraryB_name='MIP101'\n\n", RESET;


#Check user supplied options
unless ($dataname && $comp_id && $fileA && $fileB && $libraryA_name && $libraryB_name && $columnA && $columnB && $results_dir && $short_name && ($image_type =~ /tiff|jpeg|svg/i)){
  print RED, "\n\nRequired parameter missing!\n\n", RESET;
  exit();
}
chomp($dataname, $fileA, $fileB, $columnA, $columnB, $results_dir);
unless($columnA eq $columnB){
  print RED, "\n\nNOTE: columnA and columnB are not the same! ($columnA & $columnB) - continuing anyway ...\n\n", RESET;
}
$results_dir = &checkDir('-dir'=>$results_dir, '-clear'=>"no");

#Check length of each input file and abort if they are not the same
#Make sure the IDs of A match that of B (assumes that $descriptive_columns[0] contains unique IDs)
#Create the required temp file
my $temp_data_file = &processFiles();

#Now submit the required R job using calculateDifferentialExpression.R
print GREEN, "\nRunning R code:\n", RESET;

my $r_script;
if ($image_type =~ /tiff/i){
  $r_script = "$script_dir/R_bin/calculateDifferentialExpression_TIFF.R";
}elsif($image_type =~ /jpeg/i){
  $r_script = "$script_dir/R_bin/calculateDifferentialExpression_JPEG.R";
}elsif($image_type =~ /svg/i){
  $r_script = "$script_dir/R_bin/calculateDifferentialExpression_SVG.R";
}else{
  print RED, "\nImage type option not recognized\n\n", RESET;
  exit();
}
my $cmd_r = "$r_script $dataname $temp_data_file $results_dir $short_name $libraryA_name $libraryB_name";
print BLUE, "\n$cmd_r\n\n", RESET;
system($cmd_r);

#Remove the temp data file
my $cmd_rm = "rm -f $temp_data_file";
print BLUE, "\n$cmd_rm\n\n", RESET;
system($cmd_rm);

my $message = &memoryUsage;
print YELLOW, "\n\n$message\n\n", RESET;

exit();

##############################################################################################
#Process input files                                                                         #
##############################################################################################
sub processFiles{
  print GREEN, "\nProcessing input files:\n", RESET;

  my %columns_A;
  my %columns_B;

  my %ids_A;
  my %ids_B;

  #FILE_A
  open(FILE_A, "$fileA") || die "\nCould not open fileA: $fileA\n\n", RESET;
  my $header = 1;
  while(<FILE_A>){
    chomp($_);
    my @line = split("\t", $_);

    if ($header == 1){
      $header = 0;
      my $pos = 1;
      foreach my $head (@line){
        $columns_A{$head}{column_position} = $pos;
        $pos++;
      }
      next();
    }
    my $id = $line[0];
    $ids_A{$id}=1;
  }
  close(FILE_A);

  #FILE B
  open(FILE_B, "$fileB") || die "\nCould not open fileB: $fileB\n\n", RESET;
  $header = 1;
  while(<FILE_B>){
    chomp($_);
    my @line = split("\t", $_);

    if ($header == 1){
      $header = 0;
      my $pos = 1;
      foreach my $head (@line){
        $columns_B{$head}{column_position} = $pos;
        $pos++;
      }
      next();
    }
    my $id = $line[0];
    $ids_B{$id}=1;
    unless($ids_A{$id}){
      print RED, "\nFound an ID in file B that was not in file A!!\n\n", RESET;
      exit();
    }
  }
  close(FILE_B);

  #Make sure both files had the same number of entries
  my $a_count = 0;
  my $b_count = 0;
  while (my ($aid) = each %ids_A){
    $a_count++;
  }
  while (my ($bid) = each %ids_B){
    $b_count++;
  }

  unless ($a_count == $b_count){
    print RED, "\nFile A and B do not have the same number of entries!! (A = $a_count   B = $b_count)\n\n", RESET;
    exit();
  }

  #Make sure the requested columns for each file were found
  unless(($columns_A{$columnA}) && ($columns_B{$columnB})){
    print RED, "\nDid not found both ColumnA ($columnA) and ColumnB ($columnB) in FileA and FileB respectively!!\n\n", RESET;
    exit();
  }

  #Cut the required columns out of these files
  my $temp_file = "$results_dir"."$short_name"."_"."$comp_id"."_"."temp_data.txt";
  my $tmp_fileA = "$results_dir"."$short_name"."_"."$comp_id"."_"."tmpA.txt";
  my $tmp_fileB = "$results_dir"."$short_name"."_"."$comp_id"."_"."tmpB.txt";
  my $data_columnA = $columns_A{$columnA}{column_position};
  my $data_columnB = $columns_B{$columnB}{column_position};
  my $exp_columnA = $columns_A{'Expressed'}{column_position};
  my $exp_columnB = $columns_B{'Expressed'}{column_position};
  my $desc1_column = $columns_A{'FID'}{column_position};
  my $desc2_column = $columns_A{'Seq_Name'}{column_position};

  my $cmd_cutA = "cut -f 1,$desc1_column,$desc2_column,$data_columnA,$exp_columnA $fileA > $tmp_fileA";
  my $cmd_cutB = "cut -f $data_columnB,$exp_columnB $fileB > $tmp_fileB";
  my $cmd_paste = "paste $tmp_fileA $tmp_fileB > $temp_file";
  my $cmd_rm1 = "rm -f $tmp_fileA $tmp_fileB";

  print BLUE, "\n$cmd_cutA", RESET;
  system($cmd_cutA);
  print BLUE, "\n$cmd_cutB", RESET;
  system($cmd_cutB);

  #Make sure the two files have the same number of lines before attempting to paste
  my $resultA = `wc -l $tmp_fileA`;
  my $resultB = `wc -l $tmp_fileB`;
  chomp($resultA);
  chomp($resultB);
  my $linesA;
  my $linesB;
  if ($resultA =~ /^(\d+)/){
    $linesA = $1;
  }
  if ($resultB =~ /^(\d+)/){
    $linesB = $1;
  }
  unless ($linesA == $linesB){
    print RED, "\n\nDid not get the same number of lines of data from the two input files - deleting temp files and aborting...\n\n", RESET;
    system($cmd_rm1);
    exit();
  }

  print BLUE, "\n$cmd_paste", RESET;
  system($cmd_paste);
  print BLUE, "\n$cmd_rm1\n\n", RESET;
  system($cmd_rm1);

  return($temp_file);
}

