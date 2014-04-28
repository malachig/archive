#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#The purpose of this script is to summarize the expression of repeat elements obtained from 'RepBase'

#The input for this script are reads mapped to individual repeat element sequences (with a minumum length of 70 bases for simple repeats) e.g. (TGAA)n
#The outputs are a summary file describing the observed expression of each unique repeat element contained in this Repeat database

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;
use BerkeleyDB;

#Load the ALEXA modules
BEGIN {
  use Cwd 'abs_path';
  if (abs_path($0) =~ /(.*)\/.*\/.*\.pl/){
    push (@INC, $1);
  }
}
use utilities::utility qw(:all);
use utilities::mapping qw(:all);

my $repeat_library_file = '';
my $library = '';
my $read_record_dir = '';
my $mapped_reads_dir = '';
my $working_dir = '';
my $min_bit_score = '';
my $results_dir = '';
my $cutoffs_file = '';
my $log_file = '';

GetOptions ('repeat_library_file=s'=>\$repeat_library_file, 'library=s'=>\$library,
	    'mapped_reads_dir=s'=>\$mapped_reads_dir, 'read_record_dir=s'=>\$read_record_dir, 'working_dir=s'=>\$working_dir,
	    'min_bit_score=f'=>\$min_bit_score, 'results_dir=s'=>\$results_dir, 'cutoffs_file=s'=>\$cutoffs_file, 'log_file=s'=>\$log_file);

#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tThis script summarizes reads mapped to repeat elements", RESET;
print GREEN, "\n\tSpecify the complete path to a fasta file of all repeat elements, each with a header line describing them using: --repeat_library_file", RESET;
print GREEN, "\n\tSpecify the library name for read data using: --library", RESET;
print GREEN, "\n\tSpecify a directory containing read record files for this library using:  --read_record_dir", RESET;
print GREEN, "\n\tSpecify a directory contain COMPRESSED,  junction read mapping files using: --mapped_reads_dir", RESET;
print GREEN, "\n\tSpecify a working directory for temp files using:  --working_dir", RESET;
print GREEN, "\n\tSpecify a minimum bit score for each hit to be accepted using: --min_bit_score", RESET;
print GREEN, "\n\tSpecify a directory to write summary results files to using: --results_dir", RESET;
print GREEN, "\n\tSpecify the path to a file containing expression cutoffs values using:  --cutoffs_file", RESET;
print GREEN, "\n\t\tIf these have not been calculated yet, use: --cutoffs=0", RESET;
print GREEN, "\n\tSpecify a log file using: --log_file", RESET;

print GREEN, "\n\nExample: generateExpressionValues_Repeats.pl  --repeat_library_file=/projects/malachig/sequence_databases/repeat_masker_RepBase/RepBase13.05_HumanCustom/Homo_sapiens_all.txt  --library=HS04391  --read_record_dir=/projects/malachig/solexa/read_records/HS04391/  --mapped_reads_dir=/projects/malachig/solexa/read_records/HS04391/Repeats/  --working_dir=/projects/malachig/solexa/read_records/HS04391/Repeats/temp/  --min_bit_score=48.1  --results_dir=/projects/malachig/solexa/read_records/HS04391/Repeats/Summary/  --cutoffs_file=/projects/malachig/solexa/figures_and_stats/HS04391/Expression_v49/HS04391_NORM1_average_coverage_cutoffs.txt  --log_file=/projects/malachig/solexa/logs/HS04391/generateExpressionValues_Repeats_LOG.txt\n\n", RESET;

unless ($repeat_library_file && $library && $read_record_dir && $mapped_reads_dir && $working_dir && $min_bit_score && $results_dir && ($cutoffs_file || $cutoffs_file eq '0') && $log_file){
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  exit();
}

#Scale all expression values to a constant value (values will be scaled to what they would be expected to be if the library contained 10 billion mapped bases)
#If a library has less than this, the normalized expression values will be increased and vice versa
#The actual number of mapped bases is determined by parsing the 'read record' files and counting ALL bases that have been unambiguously assigned to one of the following classes:
#'ENST_U', 'INTRON_U', 'INTERGENIC_U', 'NOVEL_JUNCTION_U', 'NOVEL_BOUNDARY_U' 
my $library_normalization_value = 10000000000;
my $library_size;

open (LOG, ">$log_file") || die "\nCould not open log file: $log_file\n\n";

#Check input directories
$mapped_reads_dir = &checkDir('-dir'=>$mapped_reads_dir, '-clear'=>"no");
$working_dir = &checkDir('-dir'=>$working_dir, '-clear'=>"no");
$results_dir = &checkDir('-dir'=>$results_dir, '-clear'=>"no");

#0.) If specified, get the gene-by-gene expression cutoffs.
#    - These will be used to decide whether a particular sequence is expressed or not
my $gene_cutoffs_ref;
if ($cutoffs_file && -e $cutoffs_file){
  $gene_cutoffs_ref = &importExpressionCutoffs ('-cutoffs_file'=>$cutoffs_file);
}else{
  $cutoffs_file = 0;
  print YELLOW, "\nCutoffs file not specified - or not found, expression will be evaluated by percent coverage only\n\n", RESET;
}


#1.) Get a list of read record files corresponding to Junction mapped reads
my $mapped_reads_files_ref = &getReadFiles('-input_dir'=>$mapped_reads_dir);

#2.) Get a list of possible repeat elements, their classes and their sources
my $repeats_ref = &getRepeatElements('-repeats_file'=>$repeat_library_file);

#3.) Count read hits to each repeat element
my $repeat_hits = &parseRepeatHitRecords('-mapped_reads_files'=>$mapped_reads_files_ref, '-min_bit_score'=>$min_bit_score, '-working_dir'=>$working_dir);

#Get the size of the library by examining the read records files and counting reads that are successfully mapped
$library_size = &importLibrarySize ('-read_records_dir'=>$read_record_dir);

#4.) Summarize junction results
my $repeats_summary_file = "$results_dir"."$library"."_RepeatsExpression.txt";
&summarizeRepeatsResults('-outfile'=>$repeats_summary_file, '-library_size'=>$library_size, '-library_normalization_value'=>$library_normalization_value);

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
print LOG "\n\nMemory usage at end of script: $memory_usage_m Mb ($memory_usage_p%)"; 

print "\n\nSCRIPT COMPLETE\n\n";
print LOG "\n\nSCRIPT COMPLETE\n\n";

close(LOG);
exit();


####################################################################################################################################
#getReadFiles
####################################################################################################################################
sub getReadFiles{
  my %args = @_;
  my $input_dir = $args{'-input_dir'};

  #Get files from this directory
  print BLUE, "\n\nSearching: $input_dir for mapped reads files", RESET;
  print LOG "\n\nSearching: $input_dir for mapped reads files";

  my %read_files;
  opendir(DIRHANDLE, "$input_dir") || die "\nCannot open directory: $input_dir\n\n";
  my @test_files = readdir(DIRHANDLE);
  my $file_count = 0;
  closedir(DIRHANDLE);

  foreach my $test_file (sort @test_files){
    my $file_path = "$input_dir"."$test_file";

    #Skip directories within the specified directory
    if (-e $file_path && -d $file_path){
      print YELLOW, "\n\t$file_path  is a directory - skipping", RESET;
      print LOG "\n\t$file_path  is a directory - skipping";
      next();
    }
    $file_count++;

    #Make sure the results is compressed
    unless ($file_path =~ /(.*)\.gz$/){
      print RED, "\nFound an uncompressed file: $file_path\n\n\tMake sure all files are compressed before proceeding\n\n\t- A mix of compressed and uncompressed files may indicate a problem (i.e. you need to figure out which is complete and which might be partial!!)\n\n", RESET;
      exit();
    }

    #Check file name
    my $flowcell_lane;
    if ($test_file =~ /(\w+\_Lane\d+)/){
      $flowcell_lane = $1;
    }else{
      print RED, "\nFlowcell name in file name not understood: $test_file\n\n", RESET;
      exit();
    }
    $read_files{$file_count}{flowcell_lane} = "$flowcell_lane";
    $read_files{$file_count}{file_name} = "$test_file";
    $read_files{$file_count}{file_path} = "$file_path";

  }

  print BLUE, "\n\nFound the following $file_count files to be processed:", RESET;
  print LOG "\n\nFound the following $file_count files to be processed:";
  foreach my $file_num (sort {$read_files{$a}->{flowcell_lane} cmp $read_files{$b}->{flowcell_lane}} keys %read_files){
    print YELLOW, "\n\t$read_files{$file_num}{file_name}", RESET;
    print LOG "\n\t$read_files{$file_num}{file_name}";
  }

  #print BLUE, "\n\nDoes this seem correct?? (y/n): ", RESET;
  #my $answer = <>;
  #chomp ($answer);

  #unless ($answer =~ /^y$|^yes$/i){
  #  print RED, "\nAborting then ...\n\n", RESET;
  #  exit();
  #}

  return(\%read_files);
}


####################################################################################################################################
#Import repeats database
####################################################################################################################################
sub getRepeatElements{
  my %args = @_;
  my $infile = $args{'-repeats_file'};

  print BLUE, "\n\nParsing: $infile for basic repeat element data", RESET;
  print LOG "\n\nParsing: $infile for basic repeat element data";


  open(REPEATS, "$infile") || die "\nCould not open repeats file: $infile";
  my %repeats;
  my $repeats_count = 0;
  my $current_rep = '';

  while(<REPEATS>){
    chomp($_);
    if ($_ =~ /^\>(.*)/){
      $repeats_count++;
      my @repeat_data = split ("\t", $1);
      $repeats{$repeat_data[0]}{count} = $repeats_count;
      $repeats{$repeat_data[0]}{class} = $repeat_data[1];
      $repeats{$repeat_data[0]}{source} = $repeat_data[2];
      $repeats{$repeat_data[0]}{length} = 0;
      $repeats{$repeat_data[0]}{read_count} = 0;
      $repeats{$repeat_data[0]}{cumulative_coverage} = 0;
      $current_rep = $repeat_data[0];
    }else{
      $repeats{$current_rep}{length} += length($_);
    }
  }
  close(REPEATS);

  my $rep_count = keys %repeats;
  print BLUE, "\n\tFound $rep_count unique repeat element in this file\n\n", RESET;
  print LOG "\n\tFound $rep_count unique repeat element in this file\n\n";

  return(\%repeats);
}


#############################################################################################################################################
#parseRepeatHitRecords()                                                                                                                    #
#############################################################################################################################################
sub parseRepeatHitRecords{
  my %args = @_;
  my $files_ref = $args{'-mapped_reads_files'};
  my $min_bit_score = $args{'-min_bit_score'};
  my $working_dir = $args{'-working_dir'};

  my $total_read_hits = 0;
  my $passing_read_hits = 0;

  print BLUE, "\n\nProcessing repeat read hits from the list of input files\n", RESET;
  print LOG "\n\nProcessing repeat read hits from the list of input files\n";

  #Go through each mapped reads file and parse the count data for reads that pass quality criteria
  foreach my $file (sort {$files_ref->{$a}->{flowcell_lane} cmp $files_ref->{$b}->{flowcell_lane}} keys %{$files_ref}){
    my $file_path = $files_ref->{$file}->{file_path};

    print BLUE, "\n\tProcessing: $file_path\n\t", RESET;
    print LOG "\n\tProcessing: $file_path\n\t";

    #open a decompression file handle to this file:
    open (READS, "zcat $file_path |") || die "\nCould not open mapped reads file: $file_path\n\n";

    my $header = 1;
    my $counter = 0;

    my %columns;
    while(<READS>){
      chomp($_);
      my @line = split ("\t", $_);

      if ($header == 1){
	my $col_count = 0;
	foreach my $name (@line){
	  $columns{$name}{pos} = $col_count;
	  $col_count++;
	}
	$header = 0;
	next();
      }

      $counter++;
      if ($counter == 10000){
	$| = 1; print BLUE, ".", RESET;  $| = 0;
	print LOG ".";
	$counter = 0;
      }

      my $r1_repeat_name = $line[$columns{'R1_RepeatName'}{pos}];
      my $r2_repeat_name = $line[$columns{'R2_RepeatName'}{pos}];
      my $r1_hit_type = $line[$columns{'R1_HitType'}{pos}];
      my $r2_hit_type = $line[$columns{'R2_HitType'}{pos}];
      my $r1_strand = $line[$columns{'R1_Strand'}{pos}];
      my $r2_strand = $line[$columns{'R2_Strand'}{pos}];
      my $r1_bit_score = $line[$columns{'R1_BitScore'}{pos}];
      my $r2_bit_score = $line[$columns{'R2_BitScore'}{pos}];
      my $r1_alignment_length = $line[$columns{'R1_AlignmentLength'}{pos}];
      my $r2_alignment_length = $line[$columns{'R2_AlignmentLength'}{pos}];

      unless($r1_hit_type eq "None"){
        $total_read_hits++;
      }
      unless($r2_hit_type eq "None"){
        $total_read_hits++;
      }

      #Check format of repeat names
      if ($r1_repeat_name =~ /(.*)\|(.*)\|(.*)/){
        $r1_repeat_name = $1;
      }elsif($r1_repeat_name =~ /(.*)\|(.*)/){
        $r1_repeat_name = $1;
      }

      if ($r2_repeat_name =~ /(.*)\|(.*)\|(.*)/){
        $r2_repeat_name = $1;
      }elsif($r2_repeat_name =~ /(.*)\|(.*)/){
        $r2_repeat_name = $1;
      }

      #Deal with R1
      if ($r1_hit_type eq "Top_Hit" && $r1_bit_score >= $min_bit_score){

        #Make sure the repeat name parsed from the blast results can be found in the repeats database
        unless ($repeats_ref->{$r1_repeat_name}){
          print RED, "\nRepeat name: $r1_repeat_name from the parsed blast results does not exist (or name format has changed) in repeats database\n\n", RESET;
          exit();
        }
       
        $passing_read_hits++;
        $repeats_ref->{$r1_repeat_name}->{read_count}++;
        $repeats_ref->{$r1_repeat_name}->{cumulative_coverage} += $r1_alignment_length;
      }

      #Deal with R2
      if ($r2_hit_type eq "Top_Hit" && $r2_bit_score >= $min_bit_score){

        #Make sure the repeat name parsed from the blast results can be found in the repeats database
        unless ($repeats_ref->{$r2_repeat_name}){
          print RED, "\nRepeat name: $r2_repeat_name from the parsed blast results does not exist (or name format has changed) in repeats database\n\n", RESET;
          exit();
        }
    
        $passing_read_hits++;
        $repeats_ref->{$r2_repeat_name}->{read_count}++;
        $repeats_ref->{$r2_repeat_name}->{cumulative_coverage} += $r2_alignment_length;
      }

    }
    close(READS);
  }

  print BLUE, "\n\nA total of $total_read_hits read hits were imported and $passing_read_hits of these reads passed ('TOP_HIT' and bit_score >= $min_bit_score)\n", RESET;
  print LOG "\n\nA total of $total_read_hits read hits were imported and $passing_read_hits of these reads passed ('TOP_HIT' and bit_score >= $min_bit_score)\n";

  return(\$passing_read_hits);
}


#############################################################################################################################################
#summarizeRepeatsResults()                                                                                                                  #
#############################################################################################################################################
sub summarizeRepeatsResults{
  my %args = @_;
  my $outfile = $args{'-outfile'};
  my $library_size = $args{'-library_size'};
  my $norm_value = $args{'-library_normalization_value'};

  print BLUE, "\n\nPrinting summary file: $outfile\n", RESET;
  print LOG "\n\nPrinting summary file: $outfile\n";

  my $library_correction_ratio = 1;
  unless ($library_size == 0){
    $library_correction_ratio = $norm_value/$library_size;
  }

  #Calculate corrected read counts for all repeat elements (corrected to account for differences in library size)
  #Also calculate average coverage values (corrected for  read length) and a normalized version of this
  foreach my $repeat (keys %{$repeats_ref}){
    my $read_count_norm1 = ($repeats_ref->{$repeat}->{read_count})*$library_correction_ratio;
    $repeats_ref->{$repeat}->{read_count_norm1} = sprintf("%.2f", $read_count_norm1);

    my $average_coverage_raw = ($repeats_ref->{$repeat}->{cumulative_coverage}/$repeats_ref->{$repeat}->{length});
    my $average_coverage_norm1 = ($repeats_ref->{$repeat}->{cumulative_coverage}/$repeats_ref->{$repeat}->{length})*$library_correction_ratio;

    $repeats_ref->{$repeat}->{average_coverage_raw} = sprintf("%.10f", $average_coverage_raw);
    $repeats_ref->{$repeat}->{average_coverage_norm1} = sprintf("%.10f", $average_coverage_norm1);

    my $cutoff_test = 1;
    if ($cutoffs_file){
      my @result = @{&testExpression('-cutoffs_ref'=>$gene_cutoffs_ref, '-gene_id'=>'0', '-norm_expression_value'=>$repeats_ref->{$repeat}->{average_coverage_norm1}, '-raw_expression_value'=>$repeats_ref->{$repeat}->{average_coverage_raw}, '-percent_gene_expression_cutoff'=>0)};
      $cutoff_test = $result[0];
    }
    if ($cutoff_test == 1){
      $repeats_ref->{$repeat}->{expressed} = 1;
    }else{
      $repeats_ref->{$repeat}->{expressed} = 0;
    }
  }

  #Print out a summary of the expression observations for each repeat element
  open (OUT, ">$outfile") || die "\nCould not open outfile: $outfile\n\n";  

  #print header line
  print OUT "RepeatName\tCount\tClass\tSource\tLength\tReadCount\tReadCount_NORM1\tAverageCoverage\tAverageCoverage_NORM1\tExpressed\n";

  foreach my $repeat (sort {$repeats_ref->{$a}->{count} <=> $repeats_ref->{$b}->{count}} keys %{$repeats_ref}){

    print OUT "$repeat\t$repeats_ref->{$repeat}->{count}\t$repeats_ref->{$repeat}->{class}\t$repeats_ref->{$repeat}->{source}\t$repeats_ref->{$repeat}->{length}\t$repeats_ref->{$repeat}->{read_count}\t$repeats_ref->{$repeat}->{read_count_norm1}\t$repeats_ref->{$repeat}->{average_coverage_raw}\t$repeats_ref->{$repeat}->{average_coverage_norm1}\t$repeats_ref->{$repeat}->{expressed}\n";

  }


  close(OUT);

  return();
}


