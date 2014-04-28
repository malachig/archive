#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009-2011 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#Given a root directory for the expression data of a library, summarize the proportion of reads consumed by the top 0.1%, 1%, 2%, 5%, etc. of genes

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);

#Load the ALEXA libraries
my $script_dir;
BEGIN {
  use Cwd 'abs_path';
  if (abs_path($0) =~ /(.*)\/.*\/.*\.pl/){
    push (@INC, $1);
    $script_dir = $1;
  }
}
use utilities::utility qw(:all);

my $library_name = '';
my $ensembl_version = '';
my $analysis_dir = '';

GetOptions('library_name=s'=>\$library_name, 'ensembl_version=s'=>\$ensembl_version, 'analysis_dir=s'=>\$analysis_dir);

print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the library name using: --library_name", RESET;
print GREEN, "\n\tSpecify the EnsEMBL verion number using: --ensembl_version", RESET;
print GREEN, "\n\tSpecify the complete path to a root library expression directory using: --analysis_dir", RESET;
print GREEN, "\n\nUsage: summarizeLibrarySkew.pl  --library_name=HS118652  --ensembl_version=60  --analysis_dir=/gscmnt/gc2142/techd/analysis/alexa_seq/\n\n", RESET;

#Check user supplied options
unless ($library_name && $ensembl_version && $analysis_dir){
  print RED, "\n\nRequired parameter missing!\n\n", RESET;
  exit();
}

$analysis_dir = &checkDir('-dir'=>$analysis_dir, '-clear'=>"no");
my $library_dir = "$analysis_dir"."read_records/$library_name/";
my $working_dir = "$analysis_dir"."figures_and_stats/$library_name/temp/";
my $results_dir = "$analysis_dir"."figures_and_stats/$library_name/Expression_v"."$ensembl_version/";


#Check input directories
$library_dir = &checkDir('-dir'=>$library_dir, '-clear'=>"no");
$working_dir = &checkDir('-dir'=>$working_dir, '-clear'=>"no");
$results_dir = &checkDir('-dir'=>$results_dir, '-clear'=>"no");

#Get expression files
my %files;
my $file_count = 0;
my @data_names = qw (GeneExpression);
my @short_names = qw (Gene);

unless (scalar(@data_names) == scalar(@short_names)){
  print RED, "\n\nFull and abreviated name lists are not the same length!\n\n", RESET;
  exit();
}

foreach my $name (@data_names){
  my $file = &findFile('-root_dir'=>$library_dir, '-version'=>$ensembl_version, '-name'=>$name);
  if ($file){
    $file_count++;
    $files{$file_count}{name} = $name;
    $files{$file_count}{path} = $file;
    $files{$file_count}{temp_file} = "$working_dir"."tmp"."$file_count".".txt";
    my $short_name = shift (@short_names);
    $files{$file_count}{short_name} = $short_name;
  }
}

#Summarize files found for the user
print BLUE, "\n\nFound the following files from which data will be summarized", RESET;
foreach my $file (sort {$a <=> $b} keys %files){
  print BLUE, "\n\t$file\t$files{$file}{name}\t\t$files{$file}{path}", RESET;
}

#Grab the array of read counts
my @rc;
foreach my $file (sort {$a <=> $b} keys %files){
  my $col_name = "Read_Count";
  print BLUE, "\n\nGetting $col_name data from these files", RESET;
  open (GENE, $files{$file}{path}) || die "\n\nCould not open file: $files{$file}{path}\n\n";
  my %columns;
  my $header = 1;
  while(<GENE>){
    chomp($_);
    my @line = split("\t", $_);
    if ($header){
      my $c = 0;
      foreach my $col (@line){
        $columns{$col}{position} = $c;
        $c++;
      }
      $header = 0;
      next();
    }
    my $read_count = @line[$columns{$col_name}{position}];
    push(@rc, $read_count);
  }
  close(GENE);
}

#Sort the array of read count
my @rc_sort = sort {$b <=> $a} @rc;


#Create a file that summarizes the library expression skew as reflected by the proportion of all reads consumed by the top 0.01%, 0.1%, 1%, ... 10%, ..., 50% of genes
#PercentGenes NumberGenes  ReadCount  PercentReads

#Define the percent gene values to summarize
my $total_reads = 0;
my @percent_genes = qw (0.001 0.01 0.05 0.1 0.5 1 2 3 4 5 10 20 30 40 50);
my %cutoffs;
my %gene_counts;
my $gene_count = scalar(@rc);
my $non_zero_gene_count = 0;
foreach my $rc (@rc){
  $total_reads+=$rc;
  if ($rc > 0){
    $non_zero_gene_count++;
  }
}
for my $pg (@percent_genes){
  my $gene_number = sprintf("%0.f", ($pg/100)*$non_zero_gene_count);
  $gene_number++;
  $cutoffs{$pg}{gene_number} = $gene_number;
  $cutoffs{$pg}{gene_percent} = $pg;
  $cutoffs{$pg}{category} = "Percent Genes";
  $gene_counts{$gene_number} = 1;
}

#Also define some absolute gene counts and add to the list
my @abs_genes = qw (1 2 5 10 50 100 500 1000);
foreach my $abs_gene (@abs_genes){
  my $gene_number = $abs_gene;
  my $gene_percent = sprintf("%0.3f", (($gene_number/$non_zero_gene_count)*100));
  unless ($gene_counts{$gene_number} || defined($cutoffs{$gene_percent})){
    $cutoffs{$gene_percent}{gene_number} = $gene_number;
    $cutoffs{$gene_percent}{gene_percent} = $gene_percent;
    $cutoffs{$gene_percent}{category} = "Absolute Genes";
  }
}

#Also define some gene counts that correspond to particular round number percents of all reads consumed...
my @read_cutoffs = qw (0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.75 0.8 0.9);
my %rc_percents;
my $g = 0;
my $cum_rc = 0;
foreach my $rc (@rc_sort){
  $g++;
  $cum_rc += $rc;
  my $current_percent = $cum_rc/$total_reads;
  foreach my $r_cutoff (@read_cutoffs){
    if ($current_percent >= $r_cutoff){
      unless (defined($rc_percents{$r_cutoff})){
        $rc_percents{$r_cutoff} = 1;
        my $gene_percent = sprintf("%0.3f", (($g/$non_zero_gene_count)*100));
        unless ($gene_counts{$g} || defined($cutoffs{$gene_percent})){
          $cutoffs{$gene_percent}{gene_number} = $g;
          $cutoffs{$gene_percent}{gene_percent} = $gene_percent;
          $cutoffs{$gene_percent}{category} = "Percent Reads";
        }
      }
    }
  }
}


my $skew_file = "$results_dir"."LibraryExpressionSkew.txt";
open(SKEW, ">$skew_file") || die "\nCould not open skew file: $skew_file\n\n";

print SKEW "PercentGenes\tNumberGenes\tReadCount\tPercentReads\tNonZeroGenes\n";
foreach my $i (sort {$a <=> $b} keys %cutoffs){
  my $gene_percent = $cutoffs{$i}{gene_percent};
  my $gene_number = $cutoffs{$i}{gene_number};
  my $category = $cutoffs{$i}{category};
  my @array_slice = @rc_sort[0 .. $gene_number-1];
  my $count = scalar(@array_slice);
  my $read_sum = 0;
  foreach my $rc (@array_slice){
    $read_sum += $rc;
  }
  my $percent_reads = sprintf("%.2f", (($read_sum / $total_reads)*100));
  print BLUE, "\nPercentGenes: $gene_percent\tNumberGenes: $gene_number\tReadCount: $read_sum\tPercentReads: $percent_reads\tNonZeroGene: $non_zero_gene_count\tCategory: $category", RESET;
  print SKEW "$gene_percent\t$gene_number\t$read_sum\t$percent_reads\t$non_zero_gene_count\n";
}

close (SKEW);

print "\n\nPrinted output to: $skew_file\n\n";

exit();



######################################################################################################################################################################
#find expression file                                                                                                                                                #
######################################################################################################################################################################
sub findFile{
  my %args = @_;
  my $root_dir = $args{'-root_dir'};
  my $version = $args{'-version'};
  my $name = $args{'-name'};

  my $query = "$root_dir"."*"."/Summary/"."$library_name"."*"."_$name"."*"."_v"."$version".".txt";
  my $file = `ls $query`;
  chomp($file);
  return($file);
}


