#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#The purpose of this script is to figure out the overall average difference between exon-region expression values and exon-exon junctions or exon boundaries
#Note that exon region expression values are extracted from transcript mappings and exon junction/boundary expression values are from mappings directly to these sequences 
#For this reason, the exon region estimates have an unfair advantage in terms of coverage because coverage for a particular exon can come from reads overlaping from adjacent exons
#A junction/boundary on the other hand must have a read mapping within it and actually overlapping the junction/boundary point by at least 5-10 bases
#For these reasons, junction/boundary expression values are biased toward a lower level than exon regions

#Attempt to correct for this bias
#For each GENE, get the average normalized expression value for all exon regions
#For this same gene, get the average normalized expression value for all exon-exon junctions/boundaries (Canonical only!)
#Compare these two averages.  Determine the ratio of exon_region/exon_junction and exon_region/exon_boundary
#Now get the average (or median) of this ratio across all genes in the data set

#This final average ratio can be used to correct the exon junction/boundary expression values
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
use utilities::Descriptive;
use utilities::utility qw(:all);

my $exon_region_file = '';
my $exon_junction_file = '';
my $exon_boundary_file = '';
my $er_expression_cutoff = '';
my $library = '';
my $library_name_file = '';

GetOptions('exon_region_file=s'=>\$exon_region_file, 'exon_junction_file=s'=>\$exon_junction_file, 'exon_boundary_file=s'=>\$exon_boundary_file, 'er_expression_cutoff=f'=>\$er_expression_cutoff, 'library=s'=>\$library, 'library_name_file=s'=>\$library_name_file);

print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify a file containg exon region expression data using: --exon_region_file", RESET;
print GREEN, "\n\tSpecify a file containg exon junction expression data using: --exon_junction_file", RESET;
print GREEN, "\n\tSpecify a file containg exon boundary expression data using: --exon_boundary_file", RESET;
print GREEN, "\n\tSpecify a minimum average exon region expression cutoff for a gene to be used using: --er_expression_cutoff", RESET;
print GREEN, "\n\tSpecify the library ID using: --library", RESET;
print GREEN, "\n\tSpecify the library name config file using: --library_name_file", RESET;
print GREEN, "\n\nExample: compareExonAndJunctionExpression.pl  --exon_region_file=/projects/malachig/solexa/read_records/HS04391/ENST_v49/Summary/HS04391_Lanes1-23_ExonRegionExpression_v49.txt  --exon_junction_file=/projects/malachig/solexa/read_records/HS04391/Junctions_v49/Summary/HS04391_Lanes1-23_JunctionExpression_v49.txt  --exon_boundary_file=/projects/malachig/solexa/read_records/HS04391/Boundaries_v49/Summary/HS04391_Lanes1-23_BoundaryExpression_v49.txt  --er_expression_cutoff=4.2  --library=HS04391  --library_name_file=/projects/malachig/solexa/batch_jobs/5FU/5FU_Lib_Names.txt\n\n", RESET;

#Check user supplied options
unless ($exon_region_file && $exon_junction_file && $exon_boundary_file && ($er_expression_cutoff =~ /\d+/) && $library && $library_name_file){
  print RED, "\n\nRequired parameter missing!\n\n", RESET;
  exit();
}
unless((-e $exon_region_file) && (-e $exon_junction_file) && (-e $exon_boundary_file)){
  print RED, "\nOne of the specified input files is not found!\n\n", RESET;
  exit();
}
unless(-e $library_name_file){
  print RED, "\nLibrary name file not found!\n\n", RESET;
  exit();
}
my $er_ref = &getExpressionData('-file'=>$exon_region_file);
my $j_ref = &getExpressionData('-file'=>$exon_junction_file);
my $b_ref = &getExpressionData('-file'=>$exon_boundary_file);

#Go through each gene and get averages for each type of seq
my @er_versus_j_raw;
my @er_versus_b_raw;
my @er_versus_j_norm;
my @er_versus_b_norm;

foreach my $gene_id (sort {$a <=> $b} keys %{$er_ref}){
  my @er_data_raw = @{$er_ref->{$gene_id}->{raw_data}};
  my @er_data_norm = @{$er_ref->{$gene_id}->{norm_data}};

  my @j_data_raw;
  my @j_data_norm;
  if ($j_ref->{$gene_id}->{norm_data}){
    @j_data_raw = @{$j_ref->{$gene_id}->{raw_data}};
    @j_data_norm = @{$j_ref->{$gene_id}->{norm_data}};
  }
  my @b_data_raw;
  my @b_data_norm;
  if ($b_ref->{$gene_id}->{norm_data}){
    @b_data_raw = @{$b_ref->{$gene_id}->{raw_data}};
    @b_data_norm = @{$b_ref->{$gene_id}->{norm_data}};
  }

  my $er_mean_raw = 0;
  my $j_mean_raw = 0;
  my $b_mean_raw = 0;
  my $er_mean_norm = 0;
  my $j_mean_norm = 0;
  my $b_mean_norm = 0;

  if (scalar(@er_data_norm) > 0){
    my $er_stat_raw = Statistics::Descriptive::Full->new();
    $er_stat_raw->add_data(@er_data_raw);
    $er_mean_raw = $er_stat_raw->mean();

    my $er_stat_norm = Statistics::Descriptive::Full->new();
    $er_stat_norm->add_data(@er_data_norm);
    $er_mean_norm = $er_stat_norm->mean();
  }

  if (scalar(@j_data_norm) > 0){
    my $j_stat_raw = Statistics::Descriptive::Full->new();
    $j_stat_raw->add_data(@j_data_raw);
    $j_mean_raw = $j_stat_raw->mean();

    my $j_stat_norm = Statistics::Descriptive::Full->new();
    $j_stat_norm->add_data(@j_data_norm);
    $j_mean_norm = $j_stat_norm->mean();
  }

  if (scalar(@b_data_norm) > 0){
    my $b_stat_raw = Statistics::Descriptive::Full->new();
    $b_stat_raw->add_data(@b_data_raw);
    $b_mean_raw = $b_stat_raw->mean();

    my $b_stat_norm = Statistics::Descriptive::Full->new();
    $b_stat_norm->add_data(@b_data_norm);
    $b_mean_norm = $b_stat_norm->mean();
  }

  #Make sure some data was found for both exon regions and junctions/boundaries
  #Also, use the average exon region expression level to apply a basic filter, so that the estimates derived from this script are only based on genes that are really likely to be expressed above background
  if($er_mean_norm > $er_expression_cutoff && $j_mean_norm > 0){
    my $er_versus_j_raw = $er_mean_raw/$j_mean_raw;
    push(@er_versus_j_raw, $er_versus_j_raw);

    my $er_versus_j_norm = $er_mean_norm/$j_mean_norm;
    push(@er_versus_j_norm, $er_versus_j_norm);
  }

  if($er_mean_norm > $er_expression_cutoff && $b_mean_norm > 0){
    my $er_versus_b_raw = $er_mean_raw/$b_mean_raw;
    push(@er_versus_b_raw, $er_versus_b_raw);

    my $er_versus_b_norm = $er_mean_norm/$b_mean_norm;
    push(@er_versus_b_norm, $er_versus_b_norm);
  }
}

my $er_v_j_stat_raw = Statistics::Descriptive::Full->new();
$er_v_j_stat_raw->add_data(@er_versus_j_raw);
my $er_v_j_count_raw = $er_v_j_stat_raw->count();
my $er_v_j_mean_raw = $er_v_j_stat_raw->mean();
my $er_v_j_median_raw = $er_v_j_stat_raw->median();
my $er_v_j_min_raw = $er_v_j_stat_raw->min();
my $er_v_j_max_raw = $er_v_j_stat_raw->max();
my $er_v_j_std_raw = $er_v_j_stat_raw->standard_deviation();
print MAGENTA, "\n\nExonRegion/ExonJunction (RAW):\n\tCount = $er_v_j_count_raw\n\tMean = $er_v_j_mean_raw\n\tMedian = $er_v_j_median_raw\n\tRange = $er_v_j_min_raw - $er_v_j_max_raw\n\tStd = $er_v_j_std_raw\n", RESET; 

my $er_v_j_stat_norm = Statistics::Descriptive::Full->new();
$er_v_j_stat_norm->add_data(@er_versus_j_norm);
my $er_v_j_count_norm = $er_v_j_stat_norm->count();
my $er_v_j_mean_norm = $er_v_j_stat_norm->mean();
my $er_v_j_median_norm = $er_v_j_stat_norm->median();
my $er_v_j_min_norm = $er_v_j_stat_norm->min();
my $er_v_j_max_norm = $er_v_j_stat_norm->max();
my $er_v_j_std_norm = $er_v_j_stat_norm->standard_deviation();
print MAGENTA, "\n\nExonRegion/ExonJunction (NORM):\n\tCount = $er_v_j_count_norm\n\tMean = $er_v_j_mean_norm\n\tMedian = $er_v_j_median_norm\n\tRange = $er_v_j_min_norm - $er_v_j_max_norm\n\tStd = $er_v_j_std_norm\n", RESET; 


my $er_v_b_stat_raw = Statistics::Descriptive::Full->new();
$er_v_b_stat_raw->add_data(@er_versus_b_raw);
my $er_v_b_count_raw = $er_v_b_stat_raw->count();
my $er_v_b_mean_raw = $er_v_b_stat_raw->mean();
my $er_v_b_median_raw = $er_v_b_stat_raw->median();
my $er_v_b_min_raw = $er_v_b_stat_raw->min();
my $er_v_b_max_raw = $er_v_b_stat_raw->max();
my $er_v_b_std_raw = $er_v_b_stat_raw->standard_deviation();
print MAGENTA, "\n\nExonRegion/ExonBoundary (RAW):\n\tCount = $er_v_b_count_raw\n\tMean = $er_v_b_mean_raw\n\tMedian = $er_v_b_median_raw\n\tRange = $er_v_b_min_raw - $er_v_b_max_raw\n\tStd = $er_v_b_std_raw\n", RESET; 

my $er_v_b_stat_norm = Statistics::Descriptive::Full->new();
$er_v_b_stat_norm->add_data(@er_versus_b_norm);
my $er_v_b_count_norm = $er_v_b_stat_norm->count();
my $er_v_b_mean_norm = $er_v_b_stat_norm->mean();
my $er_v_b_median_norm = $er_v_b_stat_norm->median();
my $er_v_b_min_norm = $er_v_b_stat_norm->min();
my $er_v_b_max_norm = $er_v_b_stat_norm->max();
my $er_v_b_std_norm = $er_v_b_stat_norm->standard_deviation();
print MAGENTA, "\n\nExonRegion/ExonBoundary (NORM):\n\tCount = $er_v_b_count_norm\n\tMean = $er_v_b_mean_norm\n\tMedian = $er_v_b_median_norm\n\tRange = $er_v_b_min_norm - $er_v_b_max_norm\n\tStd = $er_v_b_std_norm\n", RESET; 

#Open the library name file and update the correction factor value in the last column
print BLUE, "\n\nWriting correction factor (median exon_region/exon_junction RAW) to specified library names file: $library_name_file\n", RESET;
open(LIBS, "$library_name_file") || die "\n\nCould not open library name file: $library_name_file\n\n";
my %libs;
my $lc = 0;
while(<LIBS>){
  chomp($_);
    if ($_ =~ /^\w+/){
    $lc++;
    my @line = split("\t", $_);
    my $lib = $line[0];
    $libs{$lib}{lc} = $lc;
    $libs{$lib}{line} = \@line;
  }
}
close(LIBS);

my $new_library_name_file = "$library_name_file".".tmp";
unless($libs{$library}){
  print RED, "\n\nCould not find library: $library in the library names file!\n\n", RESET;
  exit();
}
open(NEW, ">$new_library_name_file") || die "\n\nCould not open new library names file: $new_library_name_file\n\n"; 
foreach my $lib (sort {$libs{$a}->{lc} <=> $libs{$b}->{lc}} keys %libs){
  my @line = @{$libs{$lib}{line}};
  if ($lib eq $library){
    $line[3] = $er_v_j_median_raw;
  }
  print NEW "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\n";
}
close(NEW);

my $cmd = "mv $new_library_name_file $library_name_file";
system($cmd);

exit();


######################################################################################################################################################
#parse input expression file and store expression values for EnsEMBL supported sequences                                                             #
######################################################################################################################################################
sub getExpressionData{
  my %args = @_;
  my $file = $args{'-file'};
  
  my %data;

  my %columns;
  my $header = 1;

  my $data_count = 0;

  print BLUE, "\n\nImporting data from: $file", RESET;

  open (FILE, "$file") || die "\nCould not open input file: $file";
  while(<FILE>){
    my $line = $_;
    chomp($line);
    my @line = split("\t", $line);

    if ($header == 1){
      $header = 0;
      my $pos = 0;
      foreach my $head (@line){
        $columns{$head}{column_position} = $pos;
        $pos++;
      }
      next();
    }

    my $gene_id = $line[$columns{'Gene_ID'}{column_position}];
    my $ensembl_support = $line[$columns{'Supporting_EnsEMBL_Count'}{column_position}];
    
    my $raw_coverage = $line[$columns{'Average_Coverage_RAW'}{column_position}];
    my $norm_coverage = $line[$columns{'Average_Coverage_NORM1'}{column_position}];

    #Unless this exonRegion/junction/boundary corresponds to a known ensembl transcript, skip it
    unless ($ensembl_support >= 1){
      next();
    }
    $data_count++;

    if ($data{$gene_id}){
      push(@{$data{$gene_id}{raw_data}}, $raw_coverage);
      push(@{$data{$gene_id}{norm_data}}, $norm_coverage);
    }else{
      my @tmp1;
      my @tmp2;
      push (@tmp1, $raw_coverage);
      push (@tmp2, $norm_coverage);
      $data{$gene_id}{raw_data} = \@tmp1;
      $data{$gene_id}{norm_data} = \@tmp2;
    }
  }
  close(FILE);

  my $gene_count = keys %data;
  print BLUE, "\nStored $data_count data values corresponding to $gene_count genes\n\n", RESET;

  return(\%data);
}


