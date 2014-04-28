#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#This script takes a list of differentially or alternatively expressed exons as input (a DE or SI file)
#It then gets the sequence for these exons and prints as a fasta file to allow primer design for validation experiments
#Place PRIMER3 brackets so that the primers designed will flank the boundary between the target exon and an adjacent one 

#A number of options should be considered
#1.) Minimum exon size - should be large enough to allow design of a qPCR primer pair (this could be very limiting) - may need to consider pasting on an adjacent exon on either side (not always the valid thing to do...)
#    - Might be a good idea to skip records for small exons anyway, because these could have higher false positives (more difficult to measure)
#2.) Reciprocity - focus on exons that are DE in the opposite direction to the gene overall?
#3.) Fold-change - prioritize according to the DE of the exon?
#4.) Percent Seq DE - prioritize according to how much of the DE going on is specific to the exon versus the gene overall?
#5.) Absolute expression level of the exon in either cell line?
#6.) Percent gene expression of the exon in either cell line?

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
use utilities::ALEXA_DB qw(:all);

my $infile = '';
my $database = '';
my $server = '';
my $user = '';
my $password = '';
my $ensembl_version = '';
my $annotation_dir = '';
my $analysis_dir = '';
my $libraryA = '';
my $libraryB = '';
my $outfile = '';
my $type = '';
my $positive_only = '';
my $negative_only = '';

GetOptions ('infile=s'=>\$infile, 'database=s'=>\$database,'server=s'=>\$server, 'user=s'=>\$user, 'password=s'=>\$password, 'ensembl_version=i'=>\$ensembl_version, 
           'annotation_dir=s'=>\$annotation_dir, 'analysis_dir=s'=>\$analysis_dir,'libraryA=s'=>\$libraryA, 'libraryB=s'=>\$libraryB, 'outfile=s'=>\$outfile, 'type=s'=>\$type,
           'positive_only=i'=>\$positive_only, 'negative_only=i'=>\$negative_only);

#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tThis script produces a fasta file of sequences formated for primer design to validate DE or AE exons - supplied as a DE or SI file", RESET;
print GREEN, "\n\tSpecify the input file contain DE or AE exon regions using: --infile", RESET;
print GREEN, "\n\tSpecify whether this file is a DE or AE file using: --type=DE or --type=AE", RESET;
print GREEN, "\n\tSpecify the ALEXA database and server to query using: --database and --server", RESET;
print GREEN, "\n\tSpecify the ALEXA user and password for access using: --user and --password", RESET;
print GREEN, "\n\tSpecify the ensembl version used for this analyssi using: --ensembl_version", RESET;
print GREEN, "\n\tSpecify the directory containing ALEXA-seq annotation files using: --annotation_dir", RESET;
print GREEN, "\n\tSpecify the directory containing analysis files using: --analysis_dir", RESET;
print GREEN, "\n\tSpecify the name of libraryA using: --libraryA", RESET;
print GREEN, "\n\tSpecify the name of libraryB using: --libraryB", RESET;
print GREEN, "\n\tSpecify the name for the output file using: --outfile", RESET;
print GREEN, "\n\tTo limit targets to those with only +ve or -ve fold-changes use:  --positive_only=1  OR --negative_only=1", RESET;
print GREEN, "\n\nExample: getSequenceForPcrDesign_ExonRegion.pl  --infile=/projects/malachig/solexa/figures_and_stats/SI/5FU/ENST_v53/Mip101_vs_Mip5FuR_ExonRegion_SI_Values_Sorted_Cutoff.txt  --type=AE  --database=ALEXA_hs_53_36o  --server=jango.bcgsc.ca  --user=viewer  --password=viewer  --ensembl_version=53  --annotation_dir=/projects/malachig/sequence_databases/hs_53_36o/  --analysis_dir=/projects/malachig/solexa/  --libraryA=HS04391  --libraryB=HS04401  --outfile=/home/malachig/alternativeExons_validationSequences.txt\n\n", RESET;


unless ($infile && $database && $server && $user && $password && $ensembl_version && $annotation_dir && $analysis_dir && $libraryA && $libraryB && $outfile && $type){
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  exit();
}
unless($type =~ /DE|AE/i){
  print RED, "\nValue for parameter: --type not understood (must be DE or AE)\n\n", RESET;
  exit();
}

#Check file and directories
$annotation_dir = &checkDir('-dir'=>$annotation_dir, '-clear'=>"no");
$analysis_dir = &checkDir('-dir'=>$analysis_dir, '-clear'=>"no");


#Now go through each potential target
#First skip exon regions that do not meet various criteria. (size, reciprocal, fold-change cutoff, etc.)
my $reciprocal_test = 0;       #Event must be reciprocal or not
my $min_size = 50;             #Minimum size of the exon itself
my $min_target_size = 120;     #Minimum size of the exon or exon plus flanking exon
my $min_seq_fold_change = 3;   #Minimum fold change of the exon
my $max_targets_per_gene = 10; #Max number of exon targets allowed for a single gene
my $print_center = 1;          #Print sequence targeting the exon region only
my $print_left = 0;            #Print sequence where the target exon in joined to the exon to its left 
my $print_right = 0;           #Print sequence where the target exon in joined to the exon to its right 


#Load the input file
my %target_exon_regions;
&loadInfile('-file'=>$infile);

#Load annotations for the Exon Regions
my $er_annotation_file = "$annotation_dir"."exonRegions/exonRegions_annotated.txt.gz";
&loadAnnotations('-file'=>$er_annotation_file);

#Get the exon region (exon content blocks) sequences from the ALEXA DB - limit to those genes represented in the input file
&getGeneData();
my $gene_exon_content_ref;
my $genes_ref;

#If it passes the above filters: get the target sequence for the exon region (i.e. the exon content block that it corresponds to)
#Then grab the exon content block to the left and right
#If both are available, use the larger
open (OUT, ">$outfile") || die "\n\nCould not open output file: $outfile\n\n";

print BLUE, "\n\nConsidering each potential target exon region:\n",RESET;
my $pass_count = 0;
my $seqs_printed = 0;
foreach my $id (sort {$target_exon_regions{$b}->{si_abs} <=> $target_exon_regions{$a}->{si_abs}} keys %target_exon_regions){
  my $fid = $target_exon_regions{$id}{fid};
  my $gene_name = $target_exon_regions{$id}{gene_name};
  my $gene_id = $target_exon_regions{$id}{gene_id};
  my $seq_name = $target_exon_regions{$id}{seq_name};
  my $size = $target_exon_regions{$id}{unmasked_base_count};
  my $reciprocal = $target_exon_regions{$id}{reciprocal};
  my $seq_fold_change = $target_exon_regions{$id}{seq_fold_change};
  my $si = $target_exon_regions{$id}{si};

  #Get the sequence for the exon region
  my $start = $target_exon_regions{$id}{unit1_start};
  my $end = $target_exon_regions{$id}{unit1_end};
  my $temp_size = ($end-$start)+1;
  my $gene_seq = $genes_ref->{$gene_id}->{sequence};
  my $seq = substr ($gene_seq, $start-1, $temp_size);
  $target_exon_regions{$id}{sequence} = $seq;
  my $er_seq = $target_exon_regions{$id}{sequence};
  my $er_target_size = length($er_seq);

  print BLUE, "\n\n$id", RESET;

  if ($size >= $min_size && $reciprocal >= $reciprocal_test && abs($seq_fold_change) >= $min_seq_fold_change){
    $pass_count++;
    print BLUE, "\n\t$pass_count ($gene_name - $seq_name): Passed filter (Size: $size\tReciprocal: $reciprocal\tSeqFoldChange: $seq_fold_change)", RESET;

    #Get the exon content order from the seq name
    my $target_ec_number;
    if ($seq_name =~ /ER(\d+)\w+/){
      $target_ec_number = $1;      
    }else{
      print RED, "\n\nCould not get exon content order from Seq name: $seq_name\n\n", RESET;
      exit();
    }

    #Get the target left and right sequences
    my $left_ec_number = $target_ec_number-1;
    my $right_ec_number = $target_ec_number+1;
    my $target_ec_seq = '';
    my $left_ec_seq = '';
    my $right_ec_seq = '';
    my $exon_content_ref = $gene_exon_content_ref->{$gene_id}->{exon_content};
    foreach my $ec_id (sort {$exon_content_ref->{$a}->{start} <=> $exon_content_ref->{$b}->{start}} keys %{$exon_content_ref}){
      if ($exon_content_ref->{$ec_id}->{exon_order} == $target_ec_number){
        $target_ec_seq = $exon_content_ref->{$ec_id}->{sequence};
      }
      if ($exon_content_ref->{$ec_id}->{exon_order} == $left_ec_number){
        $left_ec_seq = $exon_content_ref->{$ec_id}->{sequence};
      }
      if ($exon_content_ref->{$ec_id}->{exon_order} == $right_ec_number){
        $right_ec_seq = $exon_content_ref->{$ec_id}->{sequence};
      }
    }

    $genes_ref->{$gene_id}->{target_seqs}++;
    if ($genes_ref->{$gene_id}->{target_seqs} > $max_targets_per_gene){
      print YELLOW, "\n\t\tAlready exceeded max targets allowed for this gene ($max_targets_per_gene)", RESET;
      next();
    }


    #Join the target sequence with the left and right sequences.  Then pick the pair that results in a larger sequence
    unless ($target_ec_seq && ($left_ec_seq || $right_ec_seq)){
      print YELLOW, "\n\nCould not retrieve exon content sequences!!\n\n", RESET;
      next();
    }
    my $target_size = length($target_ec_seq);
    my $left_join = "$left_ec_seq"."[N]"."$target_ec_seq";
    my $left_size = length($left_join);
    my $right_join = "$target_ec_seq"."[N]"."$right_ec_seq";
    my $right_size = length($right_join);
    print BLUE, "\n\t\tTarget_seq_size: $target_size\tLeft_seq_size: $left_size\tRight_seq_size: $right_size", RESET;

    #Print out larger side only
    #my $target_seq = $left_join;
    #if ($right_size > $left_size){
    #$target_seq = $right_join;
    #$seqs_printed++;
    #}

    #Print out BOTH left and right designs if available as well as the target exon region itself - depending on which options are turned on above
    if (($left_size > $target_size+3) && ($left_size > $min_target_size) && $print_left){
      $seqs_printed++;
      print OUT ">$fid\t$gene_name\t$seq_name\t$seq_fold_change\t$si\tLEFT_JOIN\n$left_join\n";
    }
    if (($right_size > $target_size+3) && ($right_size > $min_target_size) && $print_right){
      $seqs_printed++;
      print OUT ">$fid\t$gene_name\t$seq_name\t$seq_fold_change\t$si\tRIGHT_JOIN\n$right_join\n";
    }
    if (($er_target_size > $min_target_size) && $print_center){
      $seqs_printed++;
      print OUT ">$fid\t$gene_name\t$seq_name\t$seq_fold_change\t$si\tTARGET_ONLY\n$er_seq\n";
    }

    
  }else{
    print YELLOW, "\n\tDid not pass filter (Size: $size\tReciprocal: $reciprocal\tSeqFoldChange: $seq_fold_change)", RESET;
  }
}
close(OUT);

print BLUE, "\n\nPrinted a total of $seqs_printed sequences to the output file for PCR design\n\n", RESET;

exit();

####################################################################################################################
#Load input file                                                                                                   #
####################################################################################################################
sub loadInfile{
  my %args = @_;
  my $file = $args{'-file'};

  open(IN, "$file") || die "\nCould not open input file\n\n";
  my $header_line = 1;
  my %columns;

  while(<IN>){
    chomp($_);
    my @line = split("\t", $_);
    if ($header_line == 1){
      my $col_count = 0;
      foreach my $column (@line){
        $columns{$column}{position} = $col_count;
        $col_count++;
      }
      $header_line = 0;
      next();
    }
    if ($type =~ /de/i){
      if ($positive_only){
        unless ($line[$columns{Fold_Change}{position}] > 0){
          next();
        }
      }
      if ($negative_only){
        unless ($line[$columns{Fold_Change}{position}] < 0){
          next();
        }
      }

      my $id = $line[$columns{ExonRegion_ID}{position}];
      $target_exon_regions{$id}{seq_name} = $line[$columns{Seq_Name}{position}];
      $target_exon_regions{$id}{seq_fold_change} = $line[$columns{Fold_Change}{position}];
      $target_exon_regions{$id}{si} = $line[$columns{Fold_Change}{position}];
      $target_exon_regions{$id}{si_abs} = abs($line[$columns{Fold_Change}{position}]);
      $target_exon_regions{$id}{reciprocal} = 1;
    }
    if ($type =~ /ae/i){

      if ($positive_only){
        unless ($line[$columns{SEQ_Fold_Change}{position}] > 0){
          next();
        }
      }
      if ($negative_only){
        unless ($line[$columns{SEQ_Fold_Change}{position}] < 0){
          next();
        }
      }

      my $id = $line[$columns{ExonRegion_ID}{position}];
      $target_exon_regions{$id}{seq_name} = $line[$columns{Seq_Name}{position}];
      $target_exon_regions{$id}{a_gene_norm} = $line[$columns{A_GENE_Norm}{position}];
      $target_exon_regions{$id}{b_gene_norm} = $line[$columns{B_GENE_Norm}{position}];
      $target_exon_regions{$id}{a_seq_norm} = $line[$columns{A_SEQ_Norm}{position}];
      $target_exon_regions{$id}{b_seq_norm} = $line[$columns{B_SEQ_Norm}{position}];
      $target_exon_regions{$id}{gene_fold_change} = $line[$columns{GENE_Fold_Change}{position}];
      $target_exon_regions{$id}{seq_fold_change} = $line[$columns{SEQ_Fold_Change}{position}];
      $target_exon_regions{$id}{si} = $line[$columns{SI}{position}];
      $target_exon_regions{$id}{si_abs} = abs($line[$columns{SI}{position}]);
      $target_exon_regions{$id}{reciprocal} = $line[$columns{Reciprocal}{position}];
      $target_exon_regions{$id}{reciprocity} = $line[$columns{Reciprocity}{position}];
      $target_exon_regions{$id}{percent_seq_log2_de} = $line[$columns{percent_SEQ_Log2_DE}{position}];
    }

  }
  close(IN);
  return();
}


####################################################################################################################
#Load annoations file                                                                                              #
####################################################################################################################
sub loadAnnotations{
  my %args = @_;
  my $file = $args{'-file'};

  open(IN, "zcat $file |") || die "\nCould not open input file\n\n";
  my $header_line = 1;
  my %columns;

  while(<IN>){
    chomp($_);
    my @line = split("\t", $_);
    if ($header_line == 1){
      my $col_count = 0;
      foreach my $column (@line){
        $columns{$column}{position} = $col_count;
        $col_count++;
      }
      $header_line = 0;
      next();
    }
    my $id = $line[$columns{ExonRegion_ID}{position}];
    #Only store data for exon regions that might be targeted
    if ($target_exon_regions{$id}){
      $target_exon_regions{$id}{gene_id} = $line[$columns{Gene_ID}{position}];
      $target_exon_regions{$id}{gene_name} = $line[$columns{Gene_Name}{position}];
      $target_exon_regions{$id}{base_count} = $line[$columns{Base_Count}{position}];
      $target_exon_regions{$id}{unmasked_base_count} = $line[$columns{UnMasked_Base_Count}{position}];
      $target_exon_regions{$id}{fid} = $line[$columns{FID}{position}];
      my $gene_id = $line[$columns{Gene_ID}{position}];
      my $start = $line[$columns{Unit1_start}{position}];
      my $end = $line[$columns{Unit1_end}{position}];
      my $size = ($end-$start)+1;
      $target_exon_regions{$id}{unit1_start} = $start;
      $target_exon_regions{$id}{unit1_end} = $end;


    }
  }
  close(IN);
  return();
}


####################################################################################################################
#Get gene info                                                                                                     #
####################################################################################################################
sub getGeneData{

  my $alexa_dbh = &connectDB('-database'=>$database, '-server'=>$server, '-user'=>$user, '-password'=>$password);
  my %tmp;
  foreach my $id (keys %target_exon_regions){
    $tmp{$target_exon_regions{$id}{gene_id}}=1;
  }
  my @gene_ids = keys %tmp;
  $gene_exon_content_ref = &getExonContent ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids);
  $genes_ref = &getGeneInfo ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids,'-sequence'=>"yes");
  $alexa_dbh->disconnect();

  foreach my $gene_id (sort {$a <=> $b} keys %{$gene_exon_content_ref}){

    $genes_ref->{$gene_id}->{target_seqs} = 0;
    my $exon_content_ref = $gene_exon_content_ref->{$gene_id}->{exon_content};
    my $exon_content_blocks = keys %{$exon_content_ref};
    my $gene_seq = $genes_ref->{$gene_id}->{sequence};

    #1.) Assign 'exon order' to each exon content block
    my $exon_order = 0;
    foreach my $ec_id (sort {$exon_content_ref->{$a}->{start} <=> $exon_content_ref->{$b}->{start}} keys %{$exon_content_ref}){
      #Note that since we are dealing with GENE coordinates which are always stored 5' to 3', the strand does no matter at this stage
      $exon_order++;
      $exon_content_ref->{$ec_id}->{exon_order} = $exon_order;

      my $start = $exon_content_ref->{$ec_id}->{start};
      my $end = $exon_content_ref->{$ec_id}->{end};
      my $size = ($end-$start)+1;
      my $seq = substr ($gene_seq, $start-1, $size);
      $exon_content_ref->{$ec_id}->{sequence} = $seq;
     
    }
  }
  return();
}





