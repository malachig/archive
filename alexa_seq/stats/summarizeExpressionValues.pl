#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#Given a root directory for the expression data of a library, generate graphs and critical statistics for the expression of:
#-  Genes, Exon Region, Introns, Junctions, Boundaries, Introns (introns, active intron regions, silent intron regions) and Intergenics (intergenics, active intergenic regions, and silent intergenic regions)
#Generate graphs and summary statistics from the following data types:
#-  Normalized average coverage values, percent coverage 1x

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
my $library_dir = '';
my $working_dir = '';
my $results_dir = '';
my $image_type = '';

GetOptions('library_name=s'=>\$library_name, 'ensembl_version=s'=>\$ensembl_version, 'library_dir=s'=>\$library_dir, 'working_dir=s'=>\$working_dir, 'results_dir=s'=>\$results_dir, 'image_type=s'=>\$image_type);

print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the library name using: --library_name", RESET;
print GREEN, "\n\tSpecify the EnsEMBL verion number using: --ensembl_version", RESET;
print GREEN, "\n\tSpecify the complete path to a root library expression directory using: --library_dir", RESET;
print GREEN, "\n\tSpecify the complete path to a temporary working directory using: --working_dir", RESET;
print GREEN, "\n\tSpecify the complete path to a results directory using: --results_dir", RESET;
print GREEN, "\n\tSpecify the image type to generate using: --image_type=tiff or --image_type=jpeg or --image_type=svg", RESET;
print GREEN, "\n\nUsage: summarizeExpressionValues.pl  --library_name=HS04391  --ensembl_version=49  --library_dir=/projects/malachig/solexa/read_records/HS04391/  --working_dir=/projects/malachig/solexa/figures_and_stats/HS04391/temp/  --results_dir=/projects/malachig/solexa/figures_and_stats/HS04391/Expression/  --image_type=tiff\n\n", RESET;

#Check user supplied options
unless ($library_name && $ensembl_version && $library_dir && $working_dir && $results_dir && ($image_type =~ /jpeg|svg/i)){
  print RED, "\n\nRequired parameter missing!\n\n", RESET;
  exit();
}

#Check input directories
$library_dir = &checkDir('-dir'=>$library_dir, '-clear'=>"no");
$working_dir = &checkDir('-dir'=>$working_dir, '-clear'=>"no");
$results_dir = &checkDir('-dir'=>$results_dir, '-clear'=>"no");

#Get expression files
my %files;
my $file_count = 0;
my @data_names = qw (GeneExpression TranscriptExpression ExonRegionExpression JunctionExpression KnownJunctionExpression NovelJunctionExpression BoundaryExpression KnownBoundaryExpression NovelBoundaryExpression IntronExpression ActiveIntronRegionExpression SilentIntronRegionExpression IntergenicExpression ActiveIntergenicRegionExpression SilentIntergenicRegionExpression);
my @short_names = qw (Gene Transcript ExonRegion Junction KnownJunction NovelJunction Boundary KnownBoundary NovelBoundary Intron ActiveIntronRegion SilentIntronRegion Intergenic ActiveIntergenicRegion SilentIntergenicRegion);

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

#1.) Create a file to allow comparison of expression values across the various sequence types (gene, transcript, exon region, junction, etc.)
#Extract the Normalized expression data from each of these files and join into a single file
#Skip 'NA' values and events with an expression value of 0

#At the same time, calculate the grand average coverage for each sequence type
#Grand Average Coverage = Grand Cumulative Coverage / Grand Base Count (unmasked_bases)
my $grand_coverage_file = "$results_dir"."GrandAverageCoverageValues.txt";
open(GRAND_COVERAGE, ">$grand_coverage_file") || die "\nCould not open grand coverage file: $grand_coverage_file\n\n";

my $col_name = "Average_Coverage_NORM1";
print BLUE, "\n\nGetting $col_name data from these files", RESET;
my @temp_files;
my $max_lines = 0;
foreach my $file (sort {$a <=> $b} keys %files){
  my $path = $files{$file}{path};
  my $name = $files{$file}{short_name};
  my $temp_file = $files{$file}{temp_file};
  my $data_count = &extractColumnData('-path'=>$path, '-column'=>$col_name, '-temp_file'=>$temp_file, '-name'=>$name, '-file_num'=>$file);
  push(@temp_files, $temp_file);
  if ($data_count > $max_lines){
    $max_lines = $data_count;
  }
}
close (GRAND_COVERAGE);

#Join all tmp files into a single file with one column for each data set
my $joined_file = "$working_dir"."Average_Coverage_NORM1_joined.txt";
my $cmd_paste = "paste -d \',\' @temp_files > $joined_file";
system($cmd_paste);

#Clean up the temp files
my $cmd_rm = "rm -f $working_dir"."tmp*";
system($cmd_rm);


#Now do the same thing for RAW coverage data
open(GRAND_COVERAGE, ">$grand_coverage_file") || die "\nCould not open grand coverage file: $grand_coverage_file\n\n";

$col_name = "Average_Coverage_RAW";
print BLUE, "\n\nGetting $col_name data from these files", RESET;
@temp_files = ();
$max_lines = 0;
foreach my $file (sort {$a <=> $b} keys %files){
  my $path = $files{$file}{path};
  my $name = $files{$file}{short_name};
  my $temp_file = $files{$file}{temp_file};
  my $data_count = &extractColumnData('-path'=>$path, '-column'=>$col_name, '-temp_file'=>$temp_file, '-name'=>$name, '-file_num'=>$file);
  push(@temp_files, $temp_file);
  if ($data_count > $max_lines){
    $max_lines = $data_count;
  }
}
close (GRAND_COVERAGE);

#Join all tmp files into a single file with one column for each data set
my $joined_file2 = "$working_dir"."Average_Coverage_RAW_joined.txt";
$cmd_paste = "paste -d \',\' @temp_files > $joined_file2";
system($cmd_paste);

#Clean up the temp files
$cmd_rm = "rm -f $working_dir"."tmp*";
system($cmd_rm);


#2.) Create a file that will allow a direct comparison between the expression of 'Silent' Intronic regions and the Genes they reside in
#    Answer the question.  Is the expression of intronic sequence correlated to the expression level of the gene (suggesting the presence of pre-mRNA contamination)?
#    - e.g. Gene ID, Gene expression, Intron expression

#Similarly, compare this correlation to that between exon regions and genes
#This should be well correlated...
my %gene_data;
my %exon_region_data;
my %silent_intron_data;

#2-A.) Get gene expression data
foreach my $file (keys %files){
  unless ($files{$file}{name} eq "GeneExpression"){
    next(); 
  }

  my %columns = %{$files{$file}{columns}};
  open (GENE, "$files{$file}{path}") || die "\nCould not open file: $files{$file}{path}\n\n";

  my $header = 1;
  while(<GENE>){
    chomp($_);
    my @line = split ("\t", $_);
    if ($header == 1){
      $header = 0;
      next();
    }

    my $gene_id = $line[$columns{'Gene_ID'}{column_position}];
    my $expression = $line[$columns{$col_name}{column_position}];
    $gene_data{$gene_id}{expression} = $expression;
  }
  close(GENE);
}
my $gene_count = keys %gene_data;
print BLUE, "\n\nFound expression data for $gene_count genes\n", RESET;

#2-B.) Get exon region data
foreach my $file (keys %files){
  unless ($files{$file}{name} eq "ExonRegionExpression"){
    next(); 
  }

  my %columns = %{$files{$file}{columns}};
  open (ER, "$files{$file}{path}") || die "\nCould not open file: $files{$file}{path}\n\n";

  my $header = 1;
  while(<ER>){
    chomp($_);
    my @line = split ("\t", $_);
    if ($header == 1){
      $header = 0;
      next();
    }

    my $er_id = $line[$columns{'ExonRegion_ID'}{column_position}];
    my $gene_id = $line[$columns{'Gene_ID'}{column_position}];
    my $expression = $line[$columns{$col_name}{column_position}];
    $exon_region_data{$er_id}{er_expression} = $expression;
    $exon_region_data{$er_id}{gene_expression} = $gene_data{$gene_id}{expression};
  }
  close(ER);
}
my $er_count = keys %exon_region_data;
print BLUE, "\n\nFound expression data for $er_count exon regions\n", RESET;


#2-C.) Get silent intron region data
foreach my $file (keys %files){
  unless ($files{$file}{name} eq "SilentIntronRegionExpression"){
    next();
  }

  my %columns = %{$files{$file}{columns}};
  open (SR, "$files{$file}{path}") || die "\nCould not open file: $files{$file}{path}\n\n";

  my $header = 1;
  while(<SR>){
    chomp($_);
    my @line = split ("\t", $_);
    if ($header == 1){
      $header = 0;
      next();
    }

    my $sr_id = $line[$columns{'Silent_Region_ID'}{column_position}];
    my $gene_id = $line[$columns{'Gene_ID'}{column_position}];
    my $expression = $line[$columns{$col_name}{column_position}];

    #Skip silent intron regions that correspond to multiple genes
    unless ($gene_data{$gene_id}){
      next();
    }

    $silent_intron_data{$sr_id}{sr_expression} = $expression;
    $silent_intron_data{$sr_id}{gene_expression} = $gene_data{$gene_id}{expression};
  }
  close(SR);
}
my $sr_count = keys %silent_intron_data;
print BLUE, "\n\nFound expression data for $sr_count silent intron regions\n", RESET;

#2-D.) Go through the exon regions and silent intron regions and print out a temp file for importing to R
#      - Only print out introns/exon regions that are expressed above 0 and where the corresponding gene was also expressed above 0
my $er_file = "$working_dir"."Average_Coverage_NORM1_ER.txt";
my $sr_file = "$working_dir"."Average_Coverage_NORM1_SIR.txt";

open (ER_OUT, ">$er_file") || die "\nCould not open output file: $er_file\n\n";
print ER_OUT "ExonRegion_ID\tGene_Expression\tER_Expression\n";
foreach my $er_id (sort {$a cmp $b} keys %exon_region_data){
  my $gene_expression = $exon_region_data{$er_id}{gene_expression};
  my $er_expression =  $exon_region_data{$er_id}{er_expression};
  if ($gene_expression > 0 && $er_expression > 0){
    print ER_OUT "$er_id\t$gene_expression\t$er_expression\n";
  }
}
close(ER_OUT);

open (SR_OUT, ">$sr_file") || die "\nCould not open output file: $sr_file\n\n";
print SR_OUT "Silent_Region_ID\tGene_Expression\tSR_Expression\n";
foreach my $sr_id (sort {$a cmp $b} keys %silent_intron_data){
  my $gene_expression = $silent_intron_data{$sr_id}{gene_expression};
  my $sr_expression =  $silent_intron_data{$sr_id}{sr_expression};
  if ($gene_expression > 0 && $sr_expression > 0){
    print SR_OUT "$sr_id\t$gene_expression\t$sr_expression\n";
  }
}
close(SR_OUT);


#3.) Execute the R code - which will summarize the distribution of expression values for all sequence types, calculate percentiles, generate graphs, etc.
#    - This code will also print out a stats text file with basic statistics
#    - Included in this text fill will be two results which will be used to decided whether a sequence element is expressed or not:
#    - A.) The 90th percentile of all silent INTERGENIC region expression estimates
#    - B.) The coefficiencts of a linear fit to sliding window estimates of the 90th percentile for silent INTRONIC region expression estimates
my $r_script;

if($image_type =~ /jpeg/i){
  $r_script = "$script_dir/R_bin/summarizeExpressionValues_JPEG.R";
}else{
  $r_script = "$script_dir/R_bin/summarizeExpressionValues_SVG.R";
}
my $r_cmd = "$r_script $joined_file $joined_file2 $er_file $sr_file $results_dir $library_name";
print BLUE, "\n\nExecuting: $r_cmd", RESET;
system($r_cmd);


#4.) Get values from the stats file produced by the R script and use this to come up with an expression cutoff level for every gene and for non-genic elements
#    - For all genic elements (exons, junctions, boundaries, etc.) set a cutoff based on a linear fit to the 90th percetile of silent intron region expression estimates
#    - This gives a cutoff that is dependent on gene expression level.  
#    - If the gene expression value is 0 use the 95th percentile of silent intergenic region expression estimates
#    - If the gene expression low and results in a cutoff estimate lower than that from the intergenic estimate, use the intergenic estimate instead

#4-A.) First get the neccessary values generated by R in the stats file
my $stats_file = "$results_dir"."$library_name"."_stats.txt";
open (STATS, "$stats_file") || die "\nCould not open stats file: $stats_file\n\n";

my $intergenic_cutoff;
my $slope;
my $intercept;

print BLUE, "\n\nRetrieving intergenic cutoff, slope and intercept values from: $stats_file", RESET;
while(<STATS>){
  chomp($_);
  if ($_ =~ /percentile.*intergenic.*is\:\s+(\S+)\s+/){
    $intergenic_cutoff = $1;
    print BLUE, "\n\tFound intergenic cutoff: $intergenic_cutoff", RESET; 
  }
  if ($_ =~ /^Slope.*\=\s+(\S+)/){
    $slope = $1;
    print BLUE, "\n\tFound slope: $slope", RESET; 
  }
  if ($_ =~ /^Intercept.*\=\s+(\S+)/){
    $intercept = $1;
    print BLUE, "\n\tFound intercept: $intercept\n\n", RESET; 
  }
}
close(STATS);

unless (defined($intergenic_cutoff) && defined($slope) && defined($intercept)){
  print RED, "\nCould not find neccessary coefficients from stats file!\n\n", RESET;
  exit();
}

#4-B.) Now go through all genes and use the values retrieved above to set an expression cutoff for every gene.  
#      - Remember that this cutoff was determined by fitting data to log2 values but the values stored in the cutoffs file are NOT on a log2 scale
my $cutoffs_file = "$results_dir"."$library_name"."_NORM1_average_coverage_cutoffs.txt";
print BLUE, "\n\nWriting gene-by-gene cutoff values for $library_name to: $cutoffs_file\n\n", RESET;
open(CUT, ">$cutoffs_file") || die "\nCould not open cutoffs file: $cutoffs_file\n\n";
print CUT "Gene_ID\t$col_name\tIntronicCutoffEstimate\tCutoff\n";
print CUT "0\tN/A\tN/A\t$intergenic_cutoff\n";
foreach my $gene_id (sort {$a <=> $b} keys %gene_data){

  my $expression = $gene_data{$gene_id}{expression};
  my $expression_log2 = 0;
  my $fitted_y = 0;
  my $fitted_y_exp = 0;

  my $expression_cutoff = $intergenic_cutoff;
  if ($expression > 0){
    $expression_log2 =  log($expression)/log(2);

    #Use line equation from R analysis to estimate cutoff from linear fit of 95th percentiles of silent intron region expression values
    #Remember that this line was fit to data on a log2 scale
    #fitted.Y = mX + b  (m = slope, b = intercept, X = gene expression, Y = silent intron 95th percentile, fitted.Y = estimate of 95th percentile from linear regression fit)
    $fitted_y = ($slope*$expression_log2) + $intercept;

    #Convert from log2 scale back to normal scale
    $fitted_y_exp = sprintf("%.5f", 2**$fitted_y);
    if ($fitted_y_exp > $expression_cutoff){
      $expression_cutoff = $fitted_y_exp;
    }
  }
  #print YELLOW, "\n$gene_id\texp: $gene_data{$gene_id}{expression}\texp_log2: $expression_log2\tcutoff: $expression_cutoff\tslope: $slope\tintercept: $intercept\tfitted.y = $fitted_y\tfitted.y.exp = $fitted_y_exp", RESET;
  #print YELLOW, "\n$gene_id\texp: $gene_data{$gene_id}{expression}\texp_log2: $expression_log2\tfitted.y.exp = $fitted_y_exp\tcutoff: $expression_cutoff", RESET;
  print CUT "$gene_id\t$expression\t$fitted_y_exp\t$expression_cutoff\n";
}

close(CUT);

#Clean up the temp files
my $cmd_rm2 = "rm -f $working_dir"."Average*";
system($cmd_rm2);


#Generate a histogram to show the distribution of cutoff values
$r_cmd = "$script_dir/R_bin/summarizeCutoffs.R $cutoffs_file $results_dir";
print BLUE, "\n\nExecuting: $r_cmd", RESET;
system($r_cmd);

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


######################################################################################################################################################################
#find expression file                                                                                                                                                #
######################################################################################################################################################################
sub extractColumnData{
  my %args = @_;
  my $path = $args{'-path'};
  my $column = $args{'-column'};
  my $temp_file = $args{'-temp_file'};
  my $name = $args{'-name'};
  my $file_num = $args{'-file_num'};

  print MAGENTA, "\n\n$name", RESET;
  print BLUE, "\nProcessing: $path", RESET;

  my $lines_found = 0;
  my $data_found = 0;
    
  my %columns;
  my $header = 1;
  
  open (FILE, "$path") || die "\nCould not open file: $path";
  open (OUT, ">$temp_file") || die "\nCould not open temp_file: $temp_file";

  print OUT "$name\n";

  my $grand_cumulative_coverage = 0;
  my $grand_base_count = 0;

  my $grand_cumulative_coverage_known = 0;
  my $grand_base_count_known = 0;

  while(<FILE>){
    $lines_found++;
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

      #Make sure the desired column was found
      unless ($columns{$column}){
        print RED, "\nCould not file column: $column in file: $path\n\n", RESET;
        exit();
      }
      $files{$file_num}{columns} = \%columns;
      next();
    }

    my $data = $line[$columns{$column}{column_position}];
    #Skip NA values
    if ($data eq"NA"){
      next();
    }

    #Get the cumulative coverage and sequence length for this element
    my $cumulative_coverage = $line[$columns{'Cumulative_Coverage'}{column_position}];
    my $base_count = $line[$columns{'Base_Count'}{column_position}];
     
    $grand_cumulative_coverage += $cumulative_coverage;
    $grand_base_count += $base_count;

    if ($data > 0){
      print OUT "$data\n";
      $data_found++;
    }
  }

  close(FILE);
  close(OUT);

  my $grand_average_coverage = sprintf("%.5f", ($grand_cumulative_coverage/$grand_base_count));
  print BLUE, "\n\tScaned $lines_found lines for $column data and wrote $data_found data values to a temp file (non-NA and non-zero data only)", RESET;
  print BLUE, "\n\t\tGRAND AVERAGE COVERAGE:\tGrand Cumulative Coverage = $grand_cumulative_coverage\tGrand Base Count = $grand_base_count\tGrand Average Coverage = $grand_average_coverage", RESET; 

  print GRAND_COVERAGE "GRAND AVERAGE COVERAGE FOR: $name\n";
  print GRAND_COVERAGE "Grand Cumulative Coverage = $grand_cumulative_coverage\tGrand Base Count = $grand_base_count\tGrand Average Coverage = $grand_average_coverage\n\n";

  return($data_found);
}



