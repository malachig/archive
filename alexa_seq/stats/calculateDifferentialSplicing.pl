#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#This script takes expression data for two libraries ('A' and 'B') as input, calculates splicing index values, applies filters, and generates lists of significant differentially spliced 
#These calculations can be applied to the following seq types: transcripts, exon regions, junctions, boundaries, introns, silent intron regions, and active introns regions

#Four input files are required:
#i.)   Library A Gene expression file
#ii.)  Library B Gene expression file
#iii.) Library A Seq expression file
#iv.)  Library B Seq expression file

#Steps:
#1.) Import expression data for both libraries.  Check that the number of gene and seq records are equal between the two libraries
#2.) Print a temp files to a working directory
#    - Do not print records unless both the gene AND seq were found to be expressed in at least one of the two libraries

#3.) Submit these files to an R script that will do the splicing index calculations

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use IO::File;
use Tie::File;

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
my $gene_fileA = '';          #File containing gene expression values for library A
my $gene_fileB = '';          #File containing gene expression values for library B
my $exon_fileA = '';          #File containing exon expression values for library A
my $exon_fileB = '';          #File containing exon expression values for library B
my $seq_fileA = '';           #File containing seq expression values for library A
my $seq_fileB = '';           #File containing seq expression values for library B
my $data_column = '';         #Column name for data of interest 
my $results_dir = '';         #Directory for output text files and figures
my $short_name = '';
my $image_type = '';
my $comp_id = '';

GetOptions('dataname=s'=>\$dataname, 'comp_id=s'=>\$comp_id,
           'gene_fileA=s'=>\$gene_fileA, 'gene_fileB=s'=>\$gene_fileB, 
           'exon_fileA=s'=>\$exon_fileA, 'exon_fileB=s'=>\$exon_fileB, 
           'seq_fileA=s'=>\$seq_fileA, 'seq_fileB=s'=>\$seq_fileB, 'data_column=s'=>\$data_column, 
           'results_dir=s'=>\$results_dir, 'short_name=s'=>\$short_name, 'image_type=s'=>\$image_type);

print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tThis script generates SI results for Transcripts, Exons, Junctions, etc.", RESET;
print GREEN, "\n\tSpecify a descriptive name for the results using: --dataname (will be used for file names, etc.)", RESET;
print GREEN, "\n\tSpecify a comparison ID to be used in temp file naming using:  --comp_id", RESET;
print GREEN, "\n\tSpecify the complete path to a file containing GENE expression data for 'library A' using: --gene_fileA", RESET;
print GREEN, "\n\tSpecify the complete path to a file containing GENE expression data for 'library B' using: --gene_fileB", RESET;
print GREEN, "\n\tSpecify the complete path to a file containing EXON expression data for 'library A' using: --exon_fileA", RESET;
print GREEN, "\n\tSpecify the complete path to a file containing EXON expression data for 'library B' using: --exon_fileB", RESET;
print GREEN, "\n\tSpecify the complete path to a file containing SEQ expression data for 'library A' using: --seq_fileA", RESET;
print GREEN, "\n\tSpecify the complete path to a file containing SEQ expression data for 'library B' using: --seq_fileB", RESET;
print GREEN, "\n\tSpecify the data column to be used using: --data_column", RESET;
print GREEN, "\n\tSpecify the complete path to an output directory for results using: --results_dir", RESET;
print GREEN, "\n\tSpecify a short name describing the data type to be used in figure legends using:  --short_name", RESET; 
print GREEN, "\n\tSpecify the image type to generate using: --image_type='JPEG' or --image_type='TIFF' or --image_type='SVG'", RESET;


print GREEN, "\n\nUsage: calculateDifferentialSplicing.pl  --dataname='MIP101_vs_MIP5FUR_AvgCoverage'  --comp_id=HS04391_vs_HS04401  --gene_fileA=/projects/malachig/solexa/read_records/HS04391/ENST_v53/Summary/HS04391_GeneExpression_v53.txt  --gene_fileB=/projects/malachig/solexa/read_records/HS04401/ENST_v53/Summary/HS04401_GeneExpression_v53.txt  --exon_fileA=/projects/malachig/solexa/read_records/HS04391/ENST_v53/Summary/HS04391_ExonRegionExpression_v53.txt  --exon_fileB=/projects/malachig/solexa/read_records/HS04401/ENST_v53/Summary/HS04401_ExonRegionExpression_v53.txt  --seq_fileA=/projects/malachig/solexa/read_records/HS04391/Junctions_v53/Summary/HS04391_JunctionExpression_v53.txt   --seq_fileB=/projects/malachig/solexa/read_records/HS04401/Junctions_v53/Summary/HS04401_JunctionExpression_v53.txt  --data_column='Average_Coverage_RAW'   --results_dir=/projects/malachig/solexa/read_records/SI/Junctions_v53/  --short_name='Junction'  --image_type='SVG'\n\n", RESET;

#Check user supplied options
unless ($dataname && $comp_id && $gene_fileA && $gene_fileB && $exon_fileA && $exon_fileB && $seq_fileA && $seq_fileB && $data_column && $results_dir && $short_name && ($image_type =~ /tiff|jpeg|svg/i)){
  print RED, "\n\nRequired parameter missing!\n\n", RESET;
  exit();
}

$results_dir = &checkDir('-dir'=>$results_dir, '-clear'=>"no");

#Check length of each input file and abort if they are not the same
#Make sure the IDs of A match that of B
#Create the required temp file
my $temp_data_file = &processFiles();
my $message = &memoryUsage;
print YELLOW, "\n\n$message\n\n", RESET;

#Submit input file to R script which will do the calculation of SI values, p-values, multiple testing correction, etc.
my $r_script;
if ($image_type =~ /tiff/i){
  $r_script = "$script_dir/R_bin/calculateDifferentialSplicing_TIFF.R";
}elsif($image_type =~ /jpeg/i){
  $r_script = "$script_dir/R_bin/calculateDifferentialSplicing_JPEG.R";
}elsif($image_type =~ /svg/i){
  $r_script = "$script_dir/R_bin/calculateDifferentialSplicing_SVG.R";
}else{
  print RED, "\nImage type option not recognized\n\n", RESET;
  exit();
}

my $r_cmd = "$r_script $dataname $temp_data_file $results_dir $short_name";

print BLUE, "\nExecuting: $r_cmd\n\n", RESET;
system($r_cmd);

system("rm -f $temp_data_file");

$message = &memoryUsage;
print YELLOW, "\n\n$message\n\n", RESET;

exit();




##############################################################################################
#Process input files                                                                         #
##############################################################################################
sub processFiles{
  print BLUE, "\n\nProcessing input files:", RESET;

  my %gene_columns_A;
  my %gene_columns_B;
  my %seq_columns_A;
  my %seq_columns_B;
  my %exon_columns_A;
  my %exon_columns_B;

  my %genes_A;
  my %genes_B;
  my %seqs_A;
  my %seqs_B;
  my %gene_exons_A;
  my %gene_exons_B;

  my @header_names;

  #GENE FILE_A - Store gene data for library A
  print BLUE, "\n\t$gene_fileA", RESET;
  open(GENE_FILE_A, "$gene_fileA") || die "\nCould not open gene file A: $gene_fileA\n\n", RESET;
  my $header = 1;
  while(<GENE_FILE_A>){
    chomp($_);
    my @line = split("\t", $_);

    if ($header == 1){
      $header = 0;
      my $pos = 0;
      foreach my $head (@line){
        $gene_columns_A{$head}{column_position} = $pos;
        $pos++;
      }
      next();
    }

    my $id = $line[0];
    $genes_A{$id}{data} = $line[$gene_columns_A{$data_column}{column_position}];
    $genes_A{$id}{gene_id} = $line[$gene_columns_A{'Gene_ID'}{column_position}];
    $genes_A{$id}{expressed} = $line[$gene_columns_A{'Expressed'}{column_position}];
  }
  close(GENE_FILE_A);

  #GENE FILE_B - Store gene data for library A
  print BLUE, "\n\t$gene_fileB", RESET;
  open(GENE_FILE_B, "$gene_fileB") || die "\nCould not open gene file B: $gene_fileB\n\n", RESET;
  $header = 1;
  while(<GENE_FILE_B>){
    chomp($_);
    my @line = split("\t", $_);

    if ($header == 1){
      $header = 0;
      my $pos = 0;
      foreach my $head (@line){
        $gene_columns_B{$head}{column_position} = $pos;
        $pos++;
      }
      next();
    }

    my $id = $line[0];
    $genes_B{$id}{data} = $line[$gene_columns_B{$data_column}{column_position}];
    $genes_B{$id}{gene_id} = $line[$gene_columns_B{'Gene_ID'}{column_position}];
    $genes_B{$id}{expressed} = $line[$gene_columns_B{'Expressed'}{column_position}];
  }
  close(GENE_FILE_B);


  #EXON FILE_A - Store the number of exons expresed for each gene for library A
  print BLUE, "\n\t$exon_fileA", RESET;
  open(EXON_FILE_A, "$exon_fileA") || die "\nCould not open exon file A: $exon_fileA\n\n", RESET;
  $header = 1;
  while(<EXON_FILE_A>){
    chomp($_);
    my @line = split("\t", $_);

    if ($header == 1){
      $header = 0;
      my $pos = 0;
      foreach my $head (@line){
        $exon_columns_A{$head}{column_position} = $pos;
        $pos++;
      }
      next();
    }

    my $id = $line[0];
    my $gene_id = $line[$exon_columns_A{'Gene_ID'}{column_position}];
    my $expressed = $line[$exon_columns_A{'Expressed'}{column_position}];

    unless (defined($gene_exons_A{$gene_id})){
      $gene_exons_A{$gene_id}{expressed_exons} = 0;
    }
    if ($expressed == 1){
      $gene_exons_A{$gene_id}{expressed_exons}++;
    }
  }
  close(EXON_FILE_A);

  #EXON FILE_B - Store the number of exons expresed for each gene for library B
  print BLUE, "\n\t$exon_fileB", RESET;
  open(EXON_FILE_B, "$exon_fileB") || die "\nCould not open exon file B: $exon_fileB\n\n", RESET;
  $header = 1;
  while(<EXON_FILE_B>){
    chomp($_);
    my @line = split("\t", $_);

    if ($header == 1){
      $header = 0;
      my $pos = 0;
      foreach my $head (@line){
        $exon_columns_B{$head}{column_position} = $pos;
        $pos++;
      }
      next();
    }

    my $id = $line[0];
    my $gene_id = $line[$exon_columns_B{'Gene_ID'}{column_position}];
    my $expressed = $line[$exon_columns_B{'Expressed'}{column_position}];

    unless (defined($gene_exons_B{$gene_id})){
      $gene_exons_B{$gene_id}{expressed_exons} = 0;
    }
    if ($expressed == 1){
      $gene_exons_B{$gene_id}{expressed_exons}++;
    }
  }
  close(EXON_FILE_B);


  #SEQ_FILE_A and SEQ_FILE_B
  #Tie both files as an array
  print BLUE, "\n\t$seq_fileA", RESET;
  print BLUE, "\n\t$seq_fileB", RESET;
  my @seq_fileA;
  tie @seq_fileA, 'Tie::File', "$seq_fileA", mode => O_RDONLY;
  my @seq_fileB;
  tie @seq_fileB, 'Tie::File', "$seq_fileB", mode => O_RDONLY;

  #Get the header values
  my $header_A = $seq_fileA[0];
  my @line_A = split("\t", $header_A);
  my @heads_A;
  my $pos = 0;
  foreach my $head (@line_A){
    $seq_columns_A{$head}{column_position} = $pos;
    $pos++;
    push(@heads_A, $head);
  }
  push(@header_names, $heads_A[0]);
  push(@header_names, "Gene_ID");
  push(@header_names, "Seq_Name");
  push(@header_names, "FID");
  push(@header_names, "A_GENE_RAW");
  push(@header_names, "B_GENE_RAW");
  push(@header_names, "A_SEQ_RAW");
  push(@header_names, "B_SEQ_RAW");
  my $header_B = $seq_fileB[0];
  my @line_B = split("\t", $header_B);
  $pos = 0;
  foreach my $head (@line_B){
    $seq_columns_B{$head}{column_position} = $pos;
    $pos++;
  }

  #Get the line count
  my $line_count = $#seq_fileA;

  for (my $i = 1; $i <= $line_count; $i++){
    my $stringA = $seq_fileA[$i];
    my $stringB = $seq_fileB[$i];
    chomp($stringA);
    chomp($stringB);
    my @lineA = split("\t", $stringA);
    my @lineB = split("\t", $stringB);

    my $idA = $lineA[0];
    my $dataA = $lineA[$seq_columns_A{$data_column}{column_position}];
    my $gene_id = $lineA[$seq_columns_A{'Gene_ID'}{column_position}];
    my $desc1 = $lineA[$seq_columns_A{'Seq_Name'}{column_position}];
    my $desc2 = $lineA[$seq_columns_A{'FID'}{column_position}];
    my $expressedA = $lineA[$seq_columns_A{'Expressed'}{column_position}];
    
    my $idB = $lineB[0];
    my $dataB = $lineB[$seq_columns_B{$data_column}{column_position}];
    my $expressedB = $lineB[$seq_columns_B{'Expressed'}{column_position}];

    unless ($idA eq $idB){
      print RED, "\n\nSEQ IDs from each file did not match ($idA vs $idB)\n\n", RESET;
      exit();
    }

    #Some intronic elements will be defined for multiple overlapping genes - skip these
    unless ($genes_A{$gene_id}){
      next();
    }

    #If there are NA values, skip these
    if ($genes_A{$gene_id}{expressed} eq "NA" || $expressedA eq "NA" || $genes_B{$gene_id}{expressed} eq "NA" || $expressedB eq "NA"){
      next();
    }

    #Note: For very highly DE genes, we can get 'false' AE events being picked up even though the gene is effectively off in one case
    #We could require that the gene is expressed in BOTH conditions to be considered for AE, but this is probably too stringent.
    #For example, we could have a minor isoform that uses just a small subset of possible exons, this could result in the gene being not expressed overall but we dont want to lose these events
    #Therefore, require that both genes have to have at least 1 exon expressed above background.
    unless($gene_exons_A{$gene_id}{expressed_exons} > 0 && $gene_exons_B{$gene_id}{expressed_exons} > 0){
      next();
    }

    #Apply basic expression filter before printing records
    #The GENE must be considered 'expressed' in BOTH libraries
    #The SEQ must be considered 'expressed' in at least one of the two libraries
    unless (($genes_A{$gene_id}{expressed} == 1 && $genes_B{$gene_id}{expressed} == 1) && ($expressedA == 1 || $expressedB == 1)){
      next();
    }

    #Store values if they made it past the filters
    $seqs_A{$idA}{data} = $dataA;
    $seqs_A{$idA}{gene_id} = $gene_id;
    $seqs_A{$idA}{desc1} = $desc1;
    $seqs_A{$idA}{desc2} = $desc2;
    $seqs_A{$idA}{expressed} = $expressedA;

    $seqs_B{$idB}{data} = $dataB;
    $seqs_B{$idB}{expressed} = $expressedB;
  }

  #Untie both files
  untie @seq_fileA;
  untie @seq_fileB;


  #Make sure the gene, exon and seq data objects have the same number of records
  my $g_count_A = 0;
  my $g_count_B = 0;
  my $s_count_A = 0;
  my $s_count_B = 0;
  my $ge_count_A = 0;
  my $ge_count_B = 0;
  while (my ($id) = each %genes_A){$g_count_A++;}
  while (my ($id) = each %genes_B){$g_count_B++;}
  while (my ($id) = each %seqs_A){$s_count_A++;}
  while (my ($id) = each %seqs_B){$s_count_B++;}
  while (my ($id) = each %gene_exons_A){$ge_count_A++;}
  while (my ($id) = each %gene_exons_B){$ge_count_B++;}

  unless ($g_count_A == $g_count_B){
    print RED, "\n\nNumber of GENE records found for library A and B are not equal (A = $g_count_A   B = $g_count_B)\n\n", RESET;
    exit();
  }
  unless ($s_count_A == $s_count_B){
    print RED, "\n\nNumber of SEQ records found for library A and B are not equal (A = $s_count_A   B = $s_count_B)\n\n", RESET;
    exit();
  }
  unless ($ge_count_A == $ge_count_B){
    print RED, "\n\nNumber of GeneExon records found for library A and B are not equal (A = $ge_count_A   B = $ge_count_B)\n\n", RESET;
    exit();
  }

  #Write out data to the temp file
  my $temp_file = "$results_dir"."$short_name"."_"."$comp_id"."_"."temp_data.txt";

  open(TEMP, ">$temp_file") || die "\nCould not open temp file for writing: $temp_file\n\n";

  $"="\t"; print TEMP "@header_names\n"; $"=" ";

  my $total_records = 0;
  my $expressed_records_printed = 0;
  foreach my $seq_A_id (sort keys %seqs_A){
    $total_records++;
    my $gene_id = $seqs_A{$seq_A_id}{gene_id};

    $expressed_records_printed++;
    print TEMP "$seq_A_id\t$gene_id\t$seqs_A{$seq_A_id}{desc1}\t$seqs_A{$seq_A_id}{desc2}\t$genes_A{$gene_id}{data}\t$genes_B{$gene_id}{data}\t$seqs_A{$seq_A_id}{data}\t$seqs_B{$seq_A_id}{data}\n";
  }

  close(TEMP);
 
  print BLUE, "\nPrinted $expressed_records_printed (of $total_records possible) records where the Gene AND Seq were expressed in at least one of the two libraries\n\n";
  
  return($temp_file);
}










