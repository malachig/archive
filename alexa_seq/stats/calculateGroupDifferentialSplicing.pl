#!/usr/bin/perl -w
#Written by Malachi Griffith and Obi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#This script takes a pair of input files containing expression value columns (one for each library to be compared)
#It creates a temp file containing three descriptive columns and the two data columns
#It this feeds this temp file into an R script which calculates differential splicing values and assigns p-values and q-values to them

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
my $feat_data_file = '';      #Expression matrix file for feature being analysed
my $feat_expressed_file = ''; #Matrix file with status of whether feature is expressed or not (1 | 0)
my $gene_data_file = '';      #Expression matrix file for gene-level
my $gene_expressed_file = ''; #Matrix file with status of whether gene is expressed or not (1 | 0)
my $group_member_lists = '';  #List of groups and their members (e.g., GroupA:Lib1,Lib2,Lib3;GroupB:Lib4,Lib5,Lib6)
my $results_dir = '';         #Directory for output text files and figures
my $short_name = '';
my $image_type = '';
my $comp_id = '';
my $debug = '';

GetOptions('dataname=s'=>\$dataname, 'feat_data_file=s'=>\$feat_data_file, 'feat_expressed_file=s'=>\$feat_expressed_file, 'gene_data_file=s'=>\$gene_data_file, 'gene_expressed_file=s'=>\$gene_expressed_file, 'group_member_lists=s'=>\$group_member_lists, 'results_dir=s'=>\$results_dir, 'short_name=s'=>\$short_name, 'image_type=s'=>\$image_type, 'comp_id=s'=>\$comp_id, 'debug=s'=>\$debug);

print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tThis script generates SI results for Genes, Exons, Junctions, etc.", RESET;
print GREEN, "\n\tSpecify a descriptive name for the results using: --dataname (will be used for file names, etc.)", RESET;
print GREEN, "\n\tSpecify the complete path to an input expression matrix file for the correct feature type using: --feat_data_file", RESET;
print GREEN, "\n\tSpecify the complete path to an input expressed (1/0) matrix file for the correct feature type using: --feat_expressed_file", RESET;
print GREEN, "\n\tSpecify the complete path to an input expression matrix file for the gene-level using: --gene_data_file", RESET;
print GREEN, "\n\tSpecify the complete path to an input expressed (1/0) matrix file for the gene-level using: --gene_expressed_file", RESET;
print GREEN, "\n\tSpecify a comparison ID to be used in temp file naming using:  --comp_id", RESET;
print GREEN, "\n\tSpecify the groups and their members that will be compared  --group_member_lists", RESET;
print GREEN, "\n\t\tFormat:  \"GroupA:Lib1,Lib2,Lib3;GroupB:Lib4,Lib5,Lib6\"", RESET;
print GREEN, "\n\tSpecify the complete path to an output directory for results using: --results_dir", RESET;
print GREEN, "\n\tSpecify a short name describing the data type to be used in figure legends using:  --short_name", RESET;
print GREEN, "\n\tSpecify the image type to generate using: --image_type='JPEG' or --image_type='TIFF' or --image_type='SVG'", RESET;
print GREEN, "\n\tUse --debug=1 to prevent clean up of temp files", RESET;

print GREEN, "\n\nUsage: calculateGroupDifferentialSplicing.pl  --dataname=\"Basal_vs_BasalNorm_Gene\"  --feat_data_file=Matrix_TranscriptExpression_v53.txt.gz  --feat_expressed_file=Expressed_TranscriptExpression_v53.txt.gz  --gene_data_file=Matrix_GeneExpression_v53.txt.gz  --gene_expressed_file=Expressed_GeneExpression_v53.txt.gz  --group_member_lists=Basal:HCC1954,HCC3153,HCC1569,HCC70,HCC1500;BasalNorm:MCF10A,MCF10F --results_dir=/scratch/obig/ALEXA/analysis/figures_and_stats/SI/BCCL2/ENST_v53/  --short_name=Gene  --image_type=svg  --comp_id=Basal_vs_BasalNorm\n\n", RESET;

#Check user supplied options
unless ($dataname && $feat_data_file && $feat_expressed_file && $gene_data_file && $gene_expressed_file && $comp_id && $group_member_lists && $results_dir && $short_name && ($image_type =~ /tiff|jpeg|svg/i)){
  print RED, "\n\nRequired parameter missing!\n\n", RESET;
  exit();
}
chomp($dataname, $results_dir, $group_member_lists);

$results_dir = &checkDir('-dir'=>$results_dir, '-clear'=>"no");

#Unzip input matrix file if necessary.
my $temp_data_file = "$results_dir"."$short_name"."_"."$comp_id"."_"."temp_data_GSI.txt";
my $temp_expressed_file = "$results_dir"."$short_name"."_"."$comp_id"."_"."temp_expressed_GSI.txt";
my $temp_gene_data_file = "$results_dir"."$short_name"."_"."$comp_id"."_"."temp_gene_data_GSI.txt";
my $temp_gene_expressed_file = "$results_dir"."$short_name"."_"."$comp_id"."_"."temp_gene_expressed_GSI.txt";

if ($feat_data_file =~ /\.txt\.gz$/){
  my $zcat_cmd = "zcat $feat_data_file > $temp_data_file";
  system ($zcat_cmd);
}else{
  my $cp_cmd = "cp $feat_data_file $temp_data_file";
  system ($cp_cmd);
}

if ($feat_expressed_file =~ /\.txt\.gz$/){
  my $zcat_cmd2 = "zcat $feat_expressed_file > $temp_expressed_file";
  system ($zcat_cmd2);
}else{
  my $cp_cmd2 = "cp $feat_expressed_file $temp_expressed_file";
  system ($cp_cmd2);
}

if ($gene_data_file =~ /\.txt\.gz$/){
  my $zcat_cmd3 = "zcat $gene_data_file > $temp_gene_data_file";
  system ($zcat_cmd3);
}else{
  my $cp_cmd3 = "cp $gene_data_file $temp_gene_data_file";
  system ($cp_cmd3);
}

if ($gene_expressed_file =~ /\.txt\.gz$/){
  my $zcat_cmd4 = "zcat $gene_expressed_file > $temp_gene_expressed_file";
  system ($zcat_cmd4);
}else{
  my $cp_cmd4 = "cp $gene_expressed_file $temp_gene_expressed_file";
  system ($cp_cmd4);
}

#Now submit the required R job using calculateDifferentialExpression.R
print GREEN, "\nRunning R code:\n", RESET;

my $r_script;
if($image_type =~ /svg/i){
  $r_script = "$script_dir/R_bin/calculateGroupDifferentialSplicing_SVG.R";
}else{
  print RED, "\nImage type not currently supported\n\n", RESET;
  exit();
}

my $cmd_r = "$r_script $dataname $temp_data_file $temp_expressed_file $temp_gene_data_file $temp_gene_expressed_file $results_dir \"$short_name\" \"$group_member_lists\"";
print BLUE, "\n$cmd_r\n\n", RESET;
system($cmd_r);

#Remove the temp data files
my $cmd_rm = "rm -f $temp_data_file";
my $cmd_rm2 = "rm -f $temp_expressed_file";
my $cmd_rm3 = "rm -f $temp_gene_data_file";
my $cmd_rm4 = "rm -f $temp_gene_expressed_file";

if ($debug){
  print YELLOW, "\nLeaving $temp_data_file\n\n", RESET;
  print YELLOW, "\nLeaving $temp_expressed_file\n\n", RESET;
  print YELLOW, "\nLeaving $temp_gene_data_file\n\n", RESET;
  print YELLOW, "\nLeaving $temp_gene_expressed_file\n\n", RESET;
}else{
  print BLUE, "\n$cmd_rm\n\n", RESET;
  print BLUE, "\n$cmd_rm2\n\n", RESET;
  print BLUE, "\n$cmd_rm3\n\n", RESET;
  print BLUE, "\n$cmd_rm4\n\n", RESET;
  system($cmd_rm);
  system($cmd_rm2);
  system($cmd_rm3);
  system($cmd_rm4);
}

exit();

