#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#Generate expression values for transcripts by combining the expression estimates for exon regions and exon-exon junctions specific to each transcript

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use BerkeleyDB;

#Load the ALEXA modules
BEGIN {
  use Cwd 'abs_path';
  if (abs_path($0) =~ /(.*)\/.*\/.*\.pl/){
    push (@INC, $1);
  }
}
use utilities::Descriptive;
use utilities::utility qw(:all);
use utilities::ALEXA_DB qw(:all);
use utilities::mapping qw(:all);

my $ensembl_version = '';
my $database = '';
my $server = '';
my $user = '';
my $password = '';
my $annotation_dir = '';
my $exon_region_file = '';
my $exon_junction_file = '';
my $results_dir = '';
my $library_name = '';
my $cutoffs_file = '';

GetOptions('ensembl_version=s'=>\$ensembl_version, 'database=s'=>\$database,'server=s'=>\$server, 'user=s'=>\$user, 'password=s'=>\$password, 'annotation_dir=s'=>\$annotation_dir,
           'exon_region_file=s'=>\$exon_region_file, 'exon_junction_file=s'=>\$exon_junction_file,
           'results_dir=s'=>\$results_dir, 'library_name=s'=>\$library_name, 'cutoffs_file=s'=>\$cutoffs_file);

print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the EnsEMBL verion number using: --ensembl_version", RESET;
print GREEN, "\n\tSpecify the ALEXA database and server to query using: --database and --server", RESET;
print GREEN, "\n\tSpecify the ALEXA user and password for access using: --user and --password\n", RESET;
print GREEN, "\n\tSpecify a directory containing gene and exon region annotations using: --annotation_dir", RESET;
print GREEN, "\n\tSpecify a file containg exon region expression data using: --exon_region_file", RESET;
print GREEN, "\n\tSpecify a file containg exon junction expression data using: --exon_junction_file", RESET;
print GREEN, "\n\tSpecify the path to a directory for results files using: --results_dir", RESET;
print GREEN, "\n\tSpecify a library name prefix to be used for output file using:  --library_name", RESET;
print GREEN, "\n\tSpecify the path to a file containing expression cutoffs values using:  --cutoffs_file", RESET;
print GREEN, "\n\t\tIf these have not been calculated yet, use: --cutoffs=0", RESET;


print GREEN, "\n\nExample: generateExpressionValues_Transcript.pl  --ensembl_version=53  --database=ALEXA_hs_53_36o  --server=jango.bcgsc.ca  --user=viewer  --password=viewer  --annotation_dir=/projects/malachig/sequence_databases/hs_53_36o/  --exon_region_file=/projects/malachig/solexa/read_records/HS04391/ENST_v53/Summary/HS04391_Lanes1-23_ExonRegionExpression_v53.txt  --exon_junction_file=/projects/malachig/solexa/read_records/HS04391/Junctions_v53/Summary/HS04391_Lanes1-23_JunctionExpression_v53.txt  --results_dir=/projects/malachig/solexa/read_records/HS04391/Transcripts_v53/Summary/  --library_name='HS04391_Lanes1-23'  --cutoffs_file=/projects/malachig/solexa/figures_and_stats/HS04391/Expression_v53/HS04391_NORM1_average_coverage_cutoffs.txt\n\n", RESET;

#Check user supplied options
unless ($ensembl_version && $database && $server && $user && $password && $annotation_dir && $exon_region_file && $exon_junction_file && $results_dir && $library_name && ($cutoffs_file || $cutoffs_file eq '0')){
  print RED, "\n\nRequired parameter missing!\n\n", RESET;
  exit();
}

unless((-e $exon_region_file) && (-e $exon_junction_file)){
  print RED, "\nOne of the specified input files is not found!\n\n", RESET;
  exit();
}
$results_dir = &checkDir('-dir'=>$results_dir, '-clear'=>"no");
$annotation_dir = &checkDir('-dir'=>$annotation_dir, '-clear'=>"no");


#0.) If specified, get the gene-by-gene expression cutoffs.
#    - These will be used to decide whether a particular sequence is expressed or not
my $gene_cutoffs_ref;
if ($cutoffs_file && -e $cutoffs_file){
  $gene_cutoffs_ref = &importExpressionCutoffs ('-cutoffs_file'=>$cutoffs_file);
}else{
  $cutoffs_file = 0;
  print YELLOW, "\nCutoffs file not specified - or not found, expression will be evaluated by percent coverage only\n\n", RESET;
}


#1.) Get all the exon region and exon-junction expression values that correspond to specific trancsripts
my $er_ref = &getExpressionData('-file'=>$exon_region_file);
my $j_ref = &getExpressionData('-file'=>$exon_junction_file);

#2.) Import gene and transcript info from ALEXA 
my $alexa_dbh = &connectDB('-database'=>$database, '-server'=>$server, '-user'=>$user, '-password'=>$password);

my @gene_ids = @{&getAllGenes ('-dbh'=>$alexa_dbh, '-gene_type'=>'All')};
my $g_storable_name = "$database"."_AllGenes_GeneInfo_NoSeq.storable";
my $genes_ref = &getGeneInfo ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids, '-sequence'=>"no", '-storable'=>$g_storable_name);

my $trans_storable = "$database"."_AllGenes_TranscriptInfo_NoSeq.storable";
my $gene_transcripts_ref = &getTranscripts('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids, '-sequence'=>"no", '-storable'=>$trans_storable);

$alexa_dbh->disconnect();

#3.) Calculate an expression value for each transcript by averaging the expression values for the transcript specific exon regions and/or junctions

#3-A.) First store the expression data for each transcript specific sequence within the appropriate transcript object
foreach my $er (keys %{$er_ref}){
  my $raw_expression = $er_ref->{$er}->{raw_expression};
  my $norm_expression = $er_ref->{$er}->{norm_expression};
  my $gene_id = $er_ref->{$er}->{gene_id};
  my $trans_id = $er_ref->{$er}->{trans_id};
  my $bases_covered_1x = $er_ref->{$er}->{bases_covered_1x};

  my $cumulative_coverage = $er_ref->{$er}->{cumulative_coverage};

  my $trans_ref = $gene_transcripts_ref->{$gene_id}->{transcripts};
  if ($trans_ref->{$trans_id}->{norm_data}){
    push(@{$trans_ref->{$trans_id}->{raw_data}}, $raw_expression);
    push(@{$trans_ref->{$trans_id}->{norm_data}}, $norm_expression);
    $trans_ref->{$trans_id}->{specific_seq_count}++;
    $trans_ref->{$trans_id}->{specific_exon_regions}++;
    $trans_ref->{$trans_id}->{bases_covered_1x} += $bases_covered_1x;
    $trans_ref->{$trans_id}->{cumulative_coverage} += $cumulative_coverage;

  }else{
    my @tmp1;
    my @tmp2;
    push(@tmp1, $raw_expression);
    push(@tmp2, $norm_expression);
    $trans_ref->{$trans_id}->{raw_data} = \@tmp1;
    $trans_ref->{$trans_id}->{norm_data} = \@tmp2;
    $trans_ref->{$trans_id}->{specific_seq_count} = 1;
    $trans_ref->{$trans_id}->{specific_exon_regions} = 1;
    $trans_ref->{$trans_id}->{specific_exon_junctions} = 0;
    $trans_ref->{$trans_id}->{bases_covered_1x} = $bases_covered_1x;
    $trans_ref->{$trans_id}->{cumulative_coverage} = $cumulative_coverage;
  }
}

foreach my $j (keys %{$j_ref}){
  my $raw_expression = $j_ref->{$j}->{raw_expression};
  my $norm_expression = $j_ref->{$j}->{norm_expression};
  my $gene_id = $j_ref->{$j}->{gene_id};
  my $trans_id = $j_ref->{$j}->{trans_id};
  my $bases_covered_1x = $j_ref->{$j}->{bases_covered_1x};
  my $cumulative_coverage = $j_ref->{$j}->{cumulative_coverage};

  my $trans_ref = $gene_transcripts_ref->{$gene_id}->{transcripts};
  if ($trans_ref->{$trans_id}->{norm_data}){
    push(@{$trans_ref->{$trans_id}->{raw_data}}, $raw_expression);
    push(@{$trans_ref->{$trans_id}->{norm_data}}, $norm_expression);
    $trans_ref->{$trans_id}->{specific_seq_count}++;
    $trans_ref->{$trans_id}->{specific_exon_junctions}++;
    $trans_ref->{$trans_id}->{bases_covered_1x} += $bases_covered_1x;
    $trans_ref->{$trans_id}->{cumulative_coverage} += $cumulative_coverage;
  }else{
    my @tmp1;
    my @tmp2;
    push(@tmp1, $raw_expression);
    push(@tmp2, $norm_expression);
    $trans_ref->{$trans_id}->{raw_data} = \@tmp1;
    $trans_ref->{$trans_id}->{norm_data} = \@tmp2;
    $trans_ref->{$trans_id}->{specific_seq_count} = 1;
    $trans_ref->{$trans_id}->{specific_exon_junctions} = 1;
    $trans_ref->{$trans_id}->{specific_exon_regions} = 0;
    $trans_ref->{$trans_id}->{bases_covered_1x} = $bases_covered_1x;
    $trans_ref->{$trans_id}->{cumulative_coverage} = $cumulative_coverage;
  }
}


#Import the transcripts database file 
my $transcript_file = "$annotation_dir"."transcripts/transcripts_annotated.txt.gz";
unless (-e $transcript_file){
  print RED, "\nCould not find transcript annotation files:\n\t$transcript_file\n\n", RESET;
  exit();
}

#Import transcript data
my $transcript_header;
my %trans_info;
my $trans_info_ref = \%trans_info;
&importTranscripts('-infile'=>$transcript_file);




#3-B.) Now go through each transcript of each gene, calculate the average expression value
my $grand_expressed_count = 0;
foreach my $g (sort {$a <=> $b} keys %{$genes_ref}){
  my $trans_ref = $gene_transcripts_ref->{$g}->{transcripts};

  my $trans_count = keys %{$trans_ref};
  my $measurable_trans_count = 0;
  my $expressed_trans_count = 0;

  foreach my $t (sort {$a <=> $b} keys %{$trans_ref}){
    if ($trans_ref->{$t}->{norm_data}){
      $measurable_trans_count++;

      #Calculate the average expression by combining exon region and exon junction values
      my @raw_data = @{$trans_ref->{$t}->{raw_data}};
      my $raw_stat = Statistics::Descriptive::Full->new();
      $raw_stat->add_data(@raw_data);
      my $raw_mean = $raw_stat->mean();
      $trans_ref->{$t}->{raw_expression} = sprintf("%.10f", $raw_mean);

      my @norm_data = @{$trans_ref->{$t}->{norm_data}};
      my $norm_stat = Statistics::Descriptive::Full->new();
      $norm_stat->add_data(@norm_data);
      my $norm_mean = $norm_stat->mean();
      $trans_ref->{$t}->{norm_expression} = sprintf("%.10f", $norm_mean);

      #Get the average coverage for all of these elements
      my $percent_coverage_1x = sprintf("%.2f", (($trans_ref->{$t}->{bases_covered_1x}/$trans_info_ref->{$t}->{base_count})*100));
      $trans_ref->{$t}->{percent_coverage_1x} = $percent_coverage_1x;
      $trans_ref->{$t}->{expressed} = 0;

      #Determine whether this element should be considered as expressed above background
      my $cutoff_test = 1;
      my $percent_gene_expression = 0;
      if ($cutoffs_file){
        my @result = @{&testExpression('-cutoffs_ref'=>$gene_cutoffs_ref, '-gene_id'=>$g, '-norm_expression_value'=>$norm_mean, '-raw_expression_value'=>$raw_mean, '-percent_gene_expression_cutoff'=>0)};
        $cutoff_test = $result[0];
        $percent_gene_expression = $result[1];
      }
      $trans_ref->{$t}->{percent_gene_expression} = $percent_gene_expression;
      if ($percent_coverage_1x >= 75.0 && $cutoff_test == 1){
        $trans_ref->{$t}->{expressed} = 1;
        $expressed_trans_count++;
        $grand_expressed_count++;
      }
    }else{
      #Initialize variables for transcripts that did not have transcript specific elements to use...
      $trans_ref->{$t}->{specific_exon_regions} = "NA";
      $trans_ref->{$t}->{specific_exon_junctions} = "NA";
      $trans_ref->{$t}->{cumulative_coverage} = "NA";
      $trans_ref->{$t}->{raw_expression} = "NA";
      $trans_ref->{$t}->{norm_expression} = "NA";
      $trans_ref->{$t}->{bases_covered_1x} = "NA";
      $trans_ref->{$t}->{percent_coverage_1x} = "NA";
      $trans_ref->{$t}->{expressed} = "NA";
      $trans_ref->{$t}->{percent_gene_expression} = "NA";
    }
  }

  $genes_ref->{$g}->{trans_count} = $trans_count;
  $genes_ref->{$g}->{measurable_trans_count} = $measurable_trans_count;
  $genes_ref->{$g}->{expressed_trans_count} = $expressed_trans_count;

}


#4.) Print out a transcript and gene level summary

#4-A.) Gene-Level
#Gene_ID, EnsEMBL_Gene_ID, Gene_Name, Chromosome, Strand, Unit1_start_chr, Unit1_end_chr, Transcript_Count, Measurable_Transcript_Count, Expressed_Transcript_Count 
my $gene_output_file = "$results_dir"."$library_name"."_TranscriptGeneSummary_v"."$ensembl_version".".txt";
open (GENE, ">$gene_output_file") || die "\nCould not open gene output file: $gene_output_file\n\n";
print GENE "ALEXA_ID\tEnsEMBL_Gene_ID\tGene_Name\tChromosome\tStrand\tUnit1_start_chr\tUnit1_end_chr\tTranscript_Count\tMeasurable_Transcript_Count\tExpressed_Transcript_Count\n";

foreach my $g (sort {$a <=> $b} keys %{$genes_ref}){

  print GENE "$g\t$genes_ref->{$g}->{ensembl_g_id}\t$genes_ref->{$g}->{gene_name}\t$genes_ref->{$g}->{chromosome}\t$genes_ref->{$g}->{chr_strand}\t$genes_ref->{$g}->{chr_start}\t$genes_ref->{$g}->{chr_end}\t$genes_ref->{$g}->{trans_count}\t$genes_ref->{$g}->{measurable_trans_count}\t$genes_ref->{$g}->{expressed_trans_count}\n";

}
close (GENE);

#4-B.) Transcript level
my $transcript_output_file = "$results_dir"."$library_name"."_TranscriptExpression_v"."$ensembl_version".".txt";
open (TRANS, ">$transcript_output_file") || die "\nCould not open transcript output file: $transcript_output_file\n\n";
print TRANS "$transcript_header\tCumulative_Coverage\tAverage_Coverage_RAW\tAverage_Coverage_NORM1\tBases_Covered_1x\tPercent_Coverage_1x\tExpressed\tPercent_Gene_Expression\n";

foreach my $g (sort {$a <=> $b} keys %{$genes_ref}){
  my $trans_ref = $gene_transcripts_ref->{$g}->{transcripts};

  foreach my $t (sort {$a <=> $b} keys %{$trans_ref}){

    #Get chromosome coordinates for this transcript 
    my $start = $trans_ref->{$t}->{transcript_start};
    my $end = $trans_ref->{$t}->{transcript_end};
    my $coords_ref = &convertGeneCoordinates ('-gene_object'=>$genes_ref, '-gene_id'=>$g, '-start_pos'=>$start, '-end_pos'=>$end, '-ordered'=>"yes");
    my $chr_start = $coords_ref->{$g}->{chr_start};
    my $chr_end = $coords_ref->{$g}->{chr_end};

    print TRANS "$trans_info_ref->{$t}->{line}\t$trans_ref->{$t}->{cumulative_coverage}\t$trans_ref->{$t}->{raw_expression}\t$trans_ref->{$t}->{norm_expression}\t$trans_ref->{$t}->{bases_covered_1x}\t$trans_ref->{$t}->{percent_coverage_1x}\t$trans_ref->{$t}->{expressed}\t$trans_ref->{$t}->{percent_gene_expression}\n";

  }
}
close (TRANS);

#Summarize the total memory usage at close (since Perl doesnt usually release memory ... this should be the max used by the script):
my $pid = $$;
my $ps_query = `ps -p $pid -o pmem,rss`;
my @process_info = split ("\n", $ps_query);
my $memory_usage = '';
my $memory_usage_p = '';
if ($process_info[1] =~ /(\S+)\s+(\S+)/){
  $memory_usage_p = $1;
  $memory_usage = $2;
}
my $memory_usage_m = sprintf("%.1f", ($memory_usage/1024));
print YELLOW, "\n\nMemory usage at end of script: $memory_usage_m Mb ($memory_usage_p%)", RESET; 

print BLUE, "\n\nFound $grand_expressed_count expressed transcripts ... Wrote expression results to:\n\t$gene_output_file\n\t$transcript_output_file\n\n", RESET;

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

    my $seq_id;
    if ($columns{'ExonRegion_ID'}){
      $seq_id = $line[$columns{'ExonRegion_ID'}{column_position}];
    }elsif($columns{'Junction_ID'}){
      $seq_id = $line[$columns{'Junction_ID'}{column_position}];
    }else{
      print RED, "\nSeq ID not found\n\n", RESET;
      exit();
    }

    my $gene_id = $line[$columns{'Gene_ID'}{column_position}];
    my $raw_coverage = $line[$columns{'Average_Coverage_RAW'}{column_position}];
    my $norm_coverage = $line[$columns{'Average_Coverage_NORM1'}{column_position}];
    my $specific_trans_id = $line[$columns{'Specific_Trans_ID'}{column_position}];
    my $cumulative_coverage = $line[$columns{'Cumulative_Coverage'}{column_position}];
    my $bases_covered_1x = $line[$columns{'Bases_Covered_1x'}{column_position}];

    $data_count++;

    if ($specific_trans_id =~ /(\w+)\_(\d+)/){
      my $type = $1;
      my $trans_id = $2;
      $data{$seq_id}{gene_id} = $gene_id;
      $data{$seq_id}{raw_expression} = $raw_coverage;
      $data{$seq_id}{norm_expression} = $norm_coverage;
      $data{$seq_id}{type} = $type;
      $data{$seq_id}{trans_id} = $trans_id;
      $data{$seq_id}{cumulative_coverage} = $cumulative_coverage;
      $data{$seq_id}{bases_covered_1x} = $bases_covered_1x;
    }else{
      next();
    }
  }
  close(FILE);

  my $seq_count = keys %data;
  print BLUE, "\nScanned $data_count sequences and found $seq_count transcript specific sequences (exon regions or junctions)\n\n", RESET;

  return(\%data);
}


############################################################################################################################################
#Import the transcript annotation file 
############################################################################################################################################
sub importTranscripts{
  my %args = @_;
  my $transcript_file = $args{'-infile'};

  open (TRANS, "zcat $transcript_file |") || die "\nCould not open transcript annotation file: $transcript_file\n\n";  

  my $header = 1;
  my %columns;

  while(<TRANS>){
    chomp($_);
    my $line = $_;
    my @line = split("\t", $line);

    if ($header == 1){
      $transcript_header = $_;
      my $column_count = 0;
      foreach my $column (@line){
        $columns{$column}{column_pos} = $column_count;
        $column_count++;
      }
      $header = 0;
      next();
    }

    my $trans_id = $line[$columns{'Trans_ID'}{column_pos}];
    my $base_count = $line[$columns{'Base_Count'}{column_pos}];
    my $unmasked_base_count = $line[$columns{'UnMasked_Base_Count'}{column_pos}];
    my $coding_base_count = $line[$columns{'Coding_Base_Count'}{column_pos}];

    $trans_info_ref->{$trans_id}->{base_count} = $base_count;
    $trans_info_ref->{$trans_id}->{unmasked_base_count} = $unmasked_base_count;
    $trans_info_ref->{$trans_id}->{coding_base_count} = $coding_base_count;
    $trans_info_ref->{$trans_id}->{line} = $line;
  }
  close (TRANS);

  return();
}





