=head1 NAME

mapping.pm - library modules that contains generic read mapping and parsing utilities

=head1 SYNOPSIS

use mapping qw(:all);

=head2 NOTE

currently located in '~/utilities'

=head2 RECENT CHANGES

Created March 2009

=head1 DESCRIPTION

Generic utility for mapping reads and parsing results

=head1 EXAMPLES

use lib '/home/malachig/svn/solexa';

use utilities::mapping qw(:all);

=head1 SEE ALSO

None

=head1 BUGS

Contact author via email

=head1 AUTHOR

Written by Malachi Griffith (malachig@bcgsc.ca)

=head1 ACKNOWLEDGEMENTS

University of British Columbia Graduate Studies

Michael Smith Foundation for Health Research

Natural Sciences and Engineering Research Council

Genome British Columbia

=head1 AFFLIATIONS

Malachi Griffith is supervised by Marco A. Marra

Genome Sciences Centre, BC Cancer Research Centre, BC Cancer Agency, UBC Faculty of Medicine - Medical Genetics

=head1 SUBROUTINES

=cut

#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

package utilities::mapping;
require Exporter;

@ISA = qw( Exporter );
@EXPORT = qw();

@EXPORT_OK = qw(&parseBlastFiles &processReadPairBlastHits &joinRecords &convertGeneCoordinatesGeneric &removeFailingReads &importExpressionCutoffs &testExpression &importLibrarySize);

%EXPORT_TAGS = (
     all => [qw(&parseBlastFiles &processReadPairBlastHits &joinRecords &convertGeneCoordinatesGeneric &removeFailingReads &importExpressionCutoffs &testExpression &importLibrarySize)]
);

use strict;
use Data::Dumper;
use Term::ANSIColor qw(:constants);
use BerkeleyDB;

#Load the ALEXA libraries
BEGIN {
  use Cwd 'abs_path';
  if (abs_path($0) =~ /(.*)\/.*\/.*\.pl/){
    push (@INC, $1);
  }
  if (abs_path($0) =~ /(.*)\/.*\.pl/){
    push (@INC, $1);
  }
}
use utilities::utility qw(:all);
use utilities::ALEXA_DB qw(:all);


=head2 parseBlastFiles

=over 3

=item Function:

Parse a directory of compressed blast files

=item Return:

Count of top hits stored

=item Args:

'-seq_type' - Sequence type

'-blast_results_dir' - Directory containing blast results

'-top_hits_db' - Reference to Berkeley DB object

'-reads_analyzed_per_block' - User specified number of reads processed for each blast file

'-total_reads_analyzed' - User specified total number of reads analyzed for the specified lane

'-log_file_handle' - Handle for LOG file

=item Example(s):

my $top_hits_stored = &parseBlastFiles('-blast_results_dir'=>$blast_results_dir, '-top_hits_db'=>\%top_hits_db, '-min_bit_score'=>$min_bit_score, '-reads_analyzed_per_block'=>$reads_analyzed_per_block, '-total_reads_analyzed'=>$total_reads_analyzed, '-log_file_handle'=>$log_fh);

=back

=cut

sub parseBlastFiles{
  my %args = @_;
  my $seq_type = $args{'-seq_type'};
  my $blast_results_dir = $args{'-blast_results_dir'};
  my $top_hits_db_ref = $args{'-top_hits_db'};
  my $min_bit_score = $args{'-min_bit_score'};
  my $reads_analyzed_per_block = $args{'-reads_analyzed_per_block'};
  my $total_reads_analyzed = $args{'-total_reads_analyzed'};
  my $log = $args{'-log_file_handle'};
  my $seq_db_ref=$args{'-seq_db'};

  #1.) Parse blast results files one at a time
  #2.)  Print the best hit to the 'top_hits' file
  my $unambig_readpair_reads;
  my $ambig_readpair_reads;
  my $unambig_individual_reads;
  my $ambig_individual_reads;

  #Grand total counters
  my $gt_blast_results_processed = 0;
  my $gt_redundant_hits = 0;
  my $gt_blast_results_stored = 0;
  my $gt_read_ids_processed = 0;
  my $gt_unambig_readpair_reads = 0;
  my $gt_ambig_readpair_reads = 0;
  my $gt_unambig_individual_reads = 0;
  my $gt_ambig_individual_reads = 0;

  #Get files from this directory
  print BLUE, "\n\nSearching $blast_results_dir for blast results files", RESET;
  print $log "\n\nSearching $blast_results_dir for blast results files";

  my %blast_files;
  opendir(DIRHANDLE, "$blast_results_dir") || die "\nCannot open directory: $blast_results_dir\n\n";
  my @test_files = readdir(DIRHANDLE);
  my $file_count = 0;

  FILE_TEST:foreach my $test_file (sort @test_files){
      my $file_path = "$blast_results_dir"."$test_file";
      my $first_read_pair_id;

      #Skip directories within the specified directory
      if (-e $file_path && -d $file_path){
	print YELLOW, "\n\t$file_path  is a directory - skipping", RESET;
	print $log "\n\t$file_path  is a directory - skipping";
	next();
      }

      #If the results file is compressed uncompress it
      unless ($file_path =~ /(.*)\.gz$/){
        print RED, "\nFound an uncompressed file: $file_path\n\n\tMake sure all files are compressed before proceeding\n\n\t- A mix of compressed and uncompressed files may indicate a problem (i.e. you need to figure out which is complete and which might be partial!!)\n\n", RESET;
        exit();
      }

      #Check for basic format of a blast results file (12 columns, tab-delimited)
      open(TEST, "zcat $file_path | ") || die "\nCould not read file: $file_path\n\n";
      my $line_count = 0;

      while(<TEST>){
        if ($_ =~ /^\s*$/){next();}        
	if ($line_count > 0){
	  last();
	}
	chomp($_);
	my @line = split("\t", $_);

	my $read_id = $line[0];

	if ($read_id =~ /([a-zA-Z0-9]+_\d+_\d+_\d+_\d+)_(R\d)/){
	  $first_read_pair_id = $1;
	}else{
	  print RED, "\nCould not extract first read pair ID from the individual read ID: $read_id\n\n", RESET;
	  exit();
	}

	my $line_elements = scalar(@line);

	unless ($line_elements == 12){
	  print YELLOW, "\n\t$file_path does not appear to be a blast tabular output file!!\n\n", RESET;
	  print $log "\n\t$file_path does not appear to be a blast tabular output file!!\n\n";
	  close(TEST);
	  next(FILE_TEST);
	}
	$line_count++;
      }
      close(TEST);
      $file_count++;
      print BLUE, "\n\t$file_path was added to the list of files to be processed", RESET;
      print $log "\n\t$file_path was added to the list of files to be processed";

      $blast_files{$file_count}{path} = $file_path;
      $blast_files{$file_count}{first_read_pair_id} = $first_read_pair_id;
    }

  #Parse blast files one at a time and write top hits to berkley DB file as each read is processed
  my %temp_blast_results;
  my %stored_blast_results;

  my $num_files = keys %blast_files;

  unless ($num_files > 0){
    print RED, "\n\nFound no valid blast files in results directory! - Aborting\n\n", RESET;
    exit();
  }

  print BLUE, "\n\nBegin parsing $num_files blast results files", RESET;
  print $log "\n\nBegin parsing $num_files blast results files";

  foreach my $file_count (sort {$blast_files{$a}{path} cmp $blast_files{$b}{path}} keys %blast_files){

    my $blast_results_processed = 0;
    my $redundant_hits = 0;
    my $blast_results_stored = 0;
    my $read_ids_processed = 0;

    $unambig_readpair_reads = 0;
    $ambig_readpair_reads = 0;
    $unambig_individual_reads = 0;
    $ambig_individual_reads = 0;

    my $results_file = $blast_files{$file_count}{path};
    my $working_read_pair_id = $blast_files{$file_count}{first_read_pair_id};
    my $read_pair_id;

    print YELLOW, "\n\n\tParsing $results_file for blast results", RESET;
    print $log "\n\n\tParsing $results_file for blast results";

    #2a.) Go through the hits for each single read
    #    - Store the best hit
    #    - If the best and second best hit have the same alignment length and score, consider this a tie and make note of it as an ambiguous hit
    open (BLAST_RESULTS, "zcat $results_file |") || die "\nCould not open blast results file: $results_file\n\n";

    while (<BLAST_RESULTS>){
      if ($_ =~ /^\s*$/){next();}
      chomp($_);
      my @line = split ("\t", $_);

      my $read_id = $line[0];
      my $subject_id = $line[1];
      my $percent_identity = $line[2];
      my $alignment_length = $line[3];
      my $subject_start = $line[8]; #Position on the subject sequence at which the probe alignment begins
      my $subject_end = $line[9]; #Position on the subject sequence at which the probe alignment ends
      my $bit_score = $line[11];

      #Get the base read ID which corresponds to the read pair.  Both reads of the pair will have the same base read ID
      my $read_num;
      if ($read_id =~ /([a-zA-Z0-9]+_\d+_\d+_\d+_\d+)_R(\d)/){
	$read_pair_id = $1;
	$read_num = $2;
      }else{
	print RED, "\nCould not extract read pair ID from the individual read ID: $read_id\n\n", RESET;
	exit();
      }

      #Watch out for reads that appear to align beyond the range of the subject sequence (caused by BWA alignment to a concatenated DB)
      my $subject_size = $seq_db_ref->{$subject_id}->{size};
      if ($subject_start < 1 || $subject_end < 1 || $subject_start > $subject_size || $subject_end > $subject_size){
        next();
      }


      #If this read_id is different from the working read pair ID, a new block of hits has started:
      #reset the working read id, store the top N hits for this read, reset the blast_results hash, and continue
      if ($read_pair_id ne $working_read_pair_id){

	my $reads_found = keys %temp_blast_results;
	$read_ids_processed += $reads_found;
	$gt_read_ids_processed += $reads_found;

        &processReadPairBlastHits('-temp_blast_results'=>\%temp_blast_results, '-read_pair_id'=>$working_read_pair_id, '-top_hits_db'=>$top_hits_db_ref, '-min_bit_score'=>$min_bit_score,
                                  '-unambig_readpair_reads'=>\$unambig_readpair_reads, '-gt_unambig_readpair_reads'=>\$gt_unambig_readpair_reads,
                                  '-ambig_readpair_reads'=>\$ambig_readpair_reads, '-gt_ambig_readpair_reads'=>\$gt_ambig_readpair_reads,
                                  '-unambig_individual_reads'=>\$unambig_individual_reads, '-gt_unambig_individual_reads'=>\$gt_unambig_individual_reads,
                                  '-ambig_individual_reads'=>\$ambig_individual_reads, '-gt_ambig_individual_reads'=>\$gt_ambig_individual_reads);

	%temp_blast_results = ();
	$working_read_pair_id = $read_pair_id;
      }

      #NOTE: The query sequence start and end positions are always written with the start less than the end.
      #However, the subject sequences can have the start position greater than the end, if the alignment is on the opposite DNA strand
      # from the way the subject sequence was originally put into the database.
      # - Note that although they are listed in reverse order they are still relative to the start of the subject sequence on the top strand!
      # - So you should swap them before doing coordinate tests below
      my $original_start;
      my $hit_strand = "+";
      if ($subject_start > $subject_end){
	$original_start = $subject_start;
	$subject_start = $subject_end;
	$subject_end = $original_start;
        $hit_strand = "-";
      }
      $blast_results_processed++;
      $gt_blast_results_processed++;

      #Note: using the following datastructure will not allow multiple hits by one read to the same subject to be recorded
      #Only the most significant hit (longest) will be stored when this this occurs
      #If there are multiple hits with the same length (the one with a higher percent identity will be chosen)

      #See if there are any hits recorded for this read already (one set for read1 the other for read2)
      if ($temp_blast_results{$read_id}){
	my $hits_ref = $temp_blast_results{$read_id}{blast_hits};

	#See if there has already been a hit for this subject (ie. multiple hits of one read to parts of the same subject sequence)
	if ($hits_ref->{$subject_id}){

	  #Replace the old hit record for this subject if this hit has a better bit score than the one already recorded
	  if ($bit_score > $hits_ref->{$subject_id}->{bs}){
	    $hits_ref->{$subject_id}->{pi} = $percent_identity;
	    $hits_ref->{$subject_id}->{al} = $alignment_length;
	    $hits_ref->{$subject_id}->{ss} = $subject_start;
	    $hits_ref->{$subject_id}->{se} = $subject_end;
	    $hits_ref->{$subject_id}->{bs} = $bit_score;
            $hits_ref->{$subject_id}->{hs} = $hit_strand;

	  }else{
	    $redundant_hits++;
	    $gt_redundant_hits++;
	    next(); #Skip less significant hits to the same sequence
	  }

	}else{
	  #This hit is to a new subject ID
	  $hits_ref->{$subject_id}->{pi} = $percent_identity;
	  $hits_ref->{$subject_id}->{al} = $alignment_length;
	  $hits_ref->{$subject_id}->{ss} = $subject_start;
	  $hits_ref->{$subject_id}->{se} = $subject_end;
	  $hits_ref->{$subject_id}->{bs} = $bit_score;
          $hits_ref->{$subject_id}->{hs} = $hit_strand;
	  $blast_results_stored++;
	  $gt_blast_results_stored++;
	}

      }else{
	#create the first hit record for this read
	my %hits;
	$hits{$subject_id}{pi} = $percent_identity;
	$hits{$subject_id}{al} = $alignment_length;
	$hits{$subject_id}{ss} = $subject_start;
	$hits{$subject_id}{se} = $subject_end;
	$hits{$subject_id}{bs} = $bit_score;
        $hits{$subject_id}{hs} = $hit_strand;
	$temp_blast_results{$read_id}{blast_hits} = \%hits;
	$temp_blast_results{$read_id}{read_num} = $read_num;
	$blast_results_stored++;
	$gt_blast_results_stored++;
      }
    }

    #Process last block in the file!!!!
    my $reads_found = keys %temp_blast_results;
    $read_ids_processed += $reads_found;
    $gt_read_ids_processed += $reads_found;

    unless ($reads_found == 0){
      &processReadPairBlastHits('-temp_blast_results'=>\%temp_blast_results, '-read_pair_id'=>$working_read_pair_id, '-top_hits_db'=>$top_hits_db_ref, '-min_bit_score'=>$min_bit_score,
                                '-unambig_readpair_reads'=>\$unambig_readpair_reads, '-gt_unambig_readpair_reads'=>\$gt_unambig_readpair_reads,
                                '-ambig_readpair_reads'=>\$ambig_readpair_reads, '-gt_ambig_readpair_reads'=>\$gt_ambig_readpair_reads,
                                '-unambig_individual_reads'=>\$unambig_individual_reads, '-gt_unambig_individual_reads'=>\$gt_unambig_individual_reads,
                                '-ambig_individual_reads'=>\$ambig_individual_reads, '-gt_ambig_individual_reads'=>\$gt_ambig_individual_reads);
    }
    %temp_blast_results = ();
    $working_read_pair_id = $read_pair_id;

    close (BLAST_RESULTS);

    #Total possible reads that could have been mapped in this blast results file is: $reads_analyzed_per_block
    my $percent_reads_with_hits = sprintf("%.2f", (($read_ids_processed / $reads_analyzed_per_block)*100));
    my $percent_reads_unambig_as_pair = sprintf("%.2f", (($unambig_readpair_reads / $reads_analyzed_per_block)*100));
    my $percent_reads_ambig_as_pair = sprintf("%.2f", (($ambig_readpair_reads / $reads_analyzed_per_block)*100));

    my $percent_unambig_individual_reads = sprintf("%.2f", (($unambig_individual_reads / $reads_analyzed_per_block)*100));
    my $percent_ambig_individual_reads = sprintf("%.2f", (($ambig_individual_reads / $reads_analyzed_per_block)*100));

    #Summarize what was found for this blast file before proceeding to the next one
    print BLUE, "\n\tProcessed $blast_results_processed blast results", RESET;
    print BLUE, "\n\tFound $redundant_hits of equal or lesser quality to single subject sequences", RESET;
    print BLUE, "\n\tStored a total of $blast_results_stored blast results corresponding to $read_ids_processed reads (only the top hit will be printed)\n", RESET;

    print BLUE, "\n\tFound hits for $read_ids_processed ($percent_reads_with_hits%) of $reads_analyzed_per_block reads in this block of reads (block size specified by user)", RESET;
    print BLUE, "\n\tReads which could be assigned to a $seq_type as a PAIR:", RESET;
    print BLUE, "\n\t\tFound a total of $unambig_readpair_reads ($percent_reads_unambig_as_pair%) reads belonging to unambiguous read-pair mappings", RESET;
    print BLUE, "\n\t\tFound a total of $ambig_readpair_reads ($percent_reads_ambig_as_pair%) reads belonging to ambiguous read-pair mappings", RESET;
    print BLUE, "\n\tReads which had to be assigned to a $seq_type INDIVIDUALLY:", RESET;
    print BLUE, "\n\t\tFound a total of $unambig_individual_reads ($percent_unambig_individual_reads%) individual reads mapping unambiguously", RESET;
    print BLUE, "\n\t\tFound a total of $ambig_individual_reads mapping ($percent_ambig_individual_reads%) ambiguously", RESET;

    print $log "\n\tProcessed $blast_results_processed blast results";
    print $log "\n\tFound $redundant_hits of equal or lesser quality to single subject sequences";
    print $log "\n\tStored a total of $blast_results_stored blast results corresponding to $read_ids_processed reads (only the top hit will be printed)\n";

    print $log "\n\tFound hits for $read_ids_processed ($percent_reads_with_hits%) of $reads_analyzed_per_block reads in this block of reads (block size specified by user)";
    print $log "\n\tReads which could be assigned to a $seq_type as a PAIR:";
    print $log "\n\t\tFound a total of $unambig_readpair_reads ($percent_reads_unambig_as_pair%) reads belonging to unambiguous read-pair mappings";
    print $log "\n\t\tFound a total of $ambig_readpair_reads ($percent_reads_ambig_as_pair%) reads belonging to ambiguous read-pair mappings";
    print $log "\n\tReads which had to be assigned to a $seq_type INDIVIDUALLY:";
    print $log "\n\t\tFound a total of $unambig_individual_reads ($percent_unambig_individual_reads%) individual reads mapping unambiguously";
    print $log "\n\t\tFound a total of $ambig_individual_reads mapping ($percent_ambig_individual_reads%) ambiguously";

  }
  my $total_reads_with_hits = keys %{$top_hits_db_ref};
  print BLUE, "\nFound at least one hit for a grand total of $total_reads_with_hits reads\n", RESET;
  print $log "\nFound at least one hit for a grand total of $total_reads_with_hits reads\n";

  #Summarize what was found for ALL BLAST FILES - Summarize relative to all possible reads that were blasted ($total_reads_analyzed specified by the user)
  my $gt_percent_reads_with_hits = sprintf("%.2f", (($gt_read_ids_processed / $total_reads_analyzed)*100));
  my $gt_percent_reads_unambig_as_pair = sprintf("%.2f", (($gt_unambig_readpair_reads / $total_reads_analyzed)*100));
  my $gt_percent_reads_ambig_as_pair = sprintf("%.2f", (($gt_ambig_readpair_reads / $total_reads_analyzed)*100));

  my $gt_percent_unambig_individual_reads = sprintf("%.2f", (($gt_unambig_individual_reads / $total_reads_analyzed)*100));
  my $gt_percent_ambig_individual_reads = sprintf("%.2f", (($gt_ambig_individual_reads / $total_reads_analyzed)*100));

  print BLUE, "\n\nGRAND TOTAL STATS", RESET;
  print BLUE, "\n\tProcessed $gt_blast_results_processed blast results", RESET;
  print BLUE, "\n\tFound $gt_redundant_hits of equal or lesser quality to single $seq_type sequences", RESET;
  print BLUE, "\n\tStored a total of $gt_blast_results_stored blast results corresponding to $gt_read_ids_processed reads (only the top hit will be printed)\n", RESET;

  print BLUE, "\n\tFound hits for $gt_read_ids_processed ($gt_percent_reads_with_hits%) of $total_reads_analyzed total reads analyzed (specified by user)", RESET;
  print BLUE, "\n\tReads which could be assigned to a $seq_type as a PAIR:", RESET;
  print BLUE, "\n\t\tFound a total of $gt_unambig_readpair_reads ($gt_percent_reads_unambig_as_pair%) reads belonging to unambiguous read-pair mappings", RESET;
  print BLUE, "\n\t\tFound a total of $gt_ambig_readpair_reads ($gt_percent_reads_ambig_as_pair%) reads belonging to ambiguous read-pair mappings", RESET;
  print BLUE, "\n\tReads which had to be assigned to a $seq_type INDIVIDUALLY:", RESET;
  print BLUE, "\n\t\tFound a total of $gt_unambig_individual_reads ($gt_percent_unambig_individual_reads%) individual reads mapping unambiguously", RESET;
  print BLUE, "\n\t\tFound a total of $gt_ambig_individual_reads ($gt_percent_ambig_individual_reads%) reads mapping ambiguously", RESET;
  print BLUE, "\n\tNOTE: Because of cases where either R1 or R2 has NO hits, these numbers will not correspond exactly to the number of lines in your file!", RESET;

  print $log "\n\nGRAND TOTAL STATS";
  print $log "\n\tProcessed $gt_blast_results_processed blast results";
  print $log "\n\tFound $gt_redundant_hits of equal or lesser quality to single $seq_type sequences";
  print $log "\n\tStored a total of $gt_blast_results_stored blast results corresponding to $gt_read_ids_processed reads (only the top hit will be printed)\n";

  print $log "\n\tFound hits for $gt_read_ids_processed ($gt_percent_reads_with_hits%) of $total_reads_analyzed total reads analyzed (specified by user)";
  print $log "\n\tReads which could be assigned to a $seq_type as a PAIR:";
  print $log "\n\t\tFound a total of $gt_unambig_readpair_reads ($gt_percent_reads_unambig_as_pair%) reads belonging to unambiguous read-pair mappings";
  print $log "\n\t\tFound a total of $gt_ambig_readpair_reads ($gt_percent_reads_ambig_as_pair%) reads belonging to ambiguous read-pair mappings";
  print $log "\n\tReads which had to be assigned to a $seq_type INDIVIDUALLY:";
  print $log "\n\t\tFound a total of $gt_unambig_individual_reads ($gt_percent_unambig_individual_reads%) individual reads mapping unambiguously";
  print $log "\n\t\tFound a total of $gt_ambig_individual_reads ($gt_percent_ambig_individual_reads%) reads mapping ambiguously";
  print $log "\n\tNOTE: Because of cases where either R1 or R2 has NO hits, these numbers will not correspond exactly to the number of lines in your file!";

  return($total_reads_with_hits);
}


=head2 processReadPairBlastHits

=over 3

=item Function:

Process a batch of Blast hits for a single PAIRED READ mapped to a database of subject sequences (e.g. introns or intergenic regions)

=item Return:

NULL

=item Args:

'-temp_blast_results'  - blast results for a read pair

'-read_pair_id' - current read pair ID

'-top_hits_db' - Reference to Berkeley DB object

'-min_bit_score' - Minimum bit score specified by user

'-unambig_readpair_reads' - Reference to counter

'-gt_unambig_readpair_reads' - Reference to counter

'-ambig_readpair_reads' - Reference to counter

'-gt_ambig_readpair_reads' - Reference to counter

'-unambig_individual_reads' - Reference to counter

'-gt_unambig_individual_reads' - Reference to counter

'-ambig_individual_reads' - Reference to counter

'-gt_ambig_individual_reads' - Reference to counter


=item Example(s):

&processReadPairBlastHits('-temp_blast_results'=>\%temp_blast_results, '-read_pair_id'=>$working_read_pair_id, '-top_hits_db'=>$top_hits_db_ref, '-min_bit_score'=>$min_bit_score,
                          '-unambig_readpair_reads'=>\$unambig_readpair_reads, '-gt_unambig_readpair_reads'=>\$gt_unambig_readpair_reads,
                          '-ambig_readpair_reads'=>\$ambig_readpair_reads, '-gt_ambig_readpair_reads'=>\$gt_ambig_readpair_reads,
                          '-unambig_individual_reads'=>\$unambig_individual_reads, '-gt_unambig_individual_reads'=>\$gt_unambig_individual_reads,
                          '-ambig_individual_reads'=>\$ambig_individual_reads, '-gt_ambig_individual_reads'=>\$gt_ambig_individual_reads);

=back

=cut

sub processReadPairBlastHits{
  my %args = @_;
  my %read_pair_blast_results = %{$args{'-temp_blast_results'}};
  my $read_pair_id = $args{'-read_pair_id'};
  my $top_hits_db_ref = $args{'-top_hits_db'};
  my $min_bit_score = $args{'-min_bit_score'};
  my $unambig_readpair_reads_ref = $args{'-unambig_readpair_reads'};
  my $gt_unambig_readpair_reads_ref = $args{'-gt_unambig_readpair_reads'};
  my $ambig_readpair_reads_ref = $args{'-ambig_readpair_reads'};
  my $gt_ambig_readpair_reads_ref = $args{'-gt_ambig_readpair_reads'};
  my $unambig_individual_reads_ref = $args{'-unambig_individual_reads'};
  my $gt_unambig_individual_reads_ref = $args{'-gt_unambig_individual_reads'};
  my $ambig_individual_reads_ref = $args{'-ambig_individual_reads'};
  my $gt_ambig_individual_reads_ref = $args{'-gt_ambig_individual_reads'};

  #Process the blast hits found for a pair of reads

  my $read1_id = "$read_pair_id"."_R1";
  my $read2_id = "$read_pair_id"."_R2";

  #First attempt to map both reads to a subject sequence - only if this fails will they be treated as individual reads

  #Generate a complete list of subjects hit by either read
  my %subject_hit_list;
  foreach my $read (keys %read_pair_blast_results){
    my $hits_ref = $read_pair_blast_results{$read}{blast_hits};
    foreach my $subject_id (keys %{$hits_ref}){
      if ($subject_hit_list{$subject_id}){
	$subject_hit_list{$subject_id}{combined_score} += $hits_ref->{$subject_id}->{bs};
	$subject_hit_list{$subject_id}{paired_hit} = "yes";
	if ($hits_ref->{$subject_id}->{bs} < $min_bit_score){
	  $subject_hit_list{$subject_id}{both_reads_pass} = "no";
	}
      }else{
	$subject_hit_list{$subject_id}{combined_score} = $hits_ref->{$subject_id}->{bs};
	$subject_hit_list{$subject_id}{paired_hit} = "no";
	if ($hits_ref->{$subject_id}->{bs} >= $min_bit_score){
	  $subject_hit_list{$subject_id}{both_reads_pass} = "yes";
	}else{
	  $subject_hit_list{$subject_id}{both_reads_pass} = "no";
	}
      }
    }
  }

  #Determine the rank of subject pair hits according to combined bit score - where the same subject was hit by both reads
  my $rank = 0;
  my %ranked_subject_pair_hits;
  foreach my $subject_id (sort {$subject_hit_list{$b}{combined_score}<=>$subject_hit_list{$a}{combined_score}} keys %subject_hit_list){
    unless ($subject_hit_list{$subject_id}{paired_hit} eq "yes" &&  $subject_hit_list{$subject_id}{both_reads_pass} eq "yes"){
      next();
    }
    $rank++;
    $ranked_subject_pair_hits{$rank}{subject_id} = $subject_id;
    $ranked_subject_pair_hits{$rank}{combined_score} = $subject_hit_list{$subject_id}{combined_score};
  }

  my $read_pair_subject_count = keys %ranked_subject_pair_hits;

  #Possible outcomes...
  #A.) There is a single subject where both reads map to the same subject, and both have bit scores greater than some cutoff ($min_bit_score)
  #    - Print out both reads of the pair using the hits to this single subject and return().

  #B.) There are multiple subject where both reads map to a single subject, and both reads have bit scores larger than some cutoff ($min_bit_score)
  #   - If one of these subjects has a clear 'best' combined bit score (simply the highest? or at least x larger than the next best?)
  #     - Print out both reads of the pair using the hits to this single subject and return().
  #   - If multiple subjects have the same combined score, then the hit is ambiguous.
  #     - Mark both reads as ambiguous, print them out and return().

  if ($read_pair_subject_count == 1){
    #Only one subject found with quality hits to both reads
    my $subject_id = $ranked_subject_pair_hits{1}{subject_id};

    my $read1_hits_ref = $read_pair_blast_results{$read1_id}{blast_hits};
    my $read2_hits_ref = $read_pair_blast_results{$read2_id}{blast_hits};

    unless ($read1_hits_ref->{$subject_id} && $read2_hits_ref->{$subject_id}){
      print RED, "\n$read_pair_id\tSingle Paired Read - Unambiguous.  Subject ID: $subject_id\tREAD1: $read1_id\tREAD2: $read2_id", RESET;
      print Dumper %read_pair_blast_results;
      print Dumper %subject_hit_list;
      print Dumper %ranked_subject_pair_hits;
      exit();
    }

    #    - Read_ID, Hit_Type (i.e. Ambiguous or Top Hit), Subject_ID, AlignmentLength, PercentIdentity, BitScore, SubjectStart, SubjectEnd, strand

    $top_hits_db_ref->{$read1_id} = "Top_Hit\t$subject_id\t$read1_hits_ref->{$subject_id}->{al}\t$read1_hits_ref->{$subject_id}->{pi}\t$read1_hits_ref->{$subject_id}->{bs}\t$read1_hits_ref->{$subject_id}->{ss}\t$read1_hits_ref->{$subject_id}->{se}\t$read1_hits_ref->{$subject_id}->{hs}";

    $top_hits_db_ref->{$read2_id} = "Top_Hit\t$subject_id\t$read2_hits_ref->{$subject_id}->{al}\t$read2_hits_ref->{$subject_id}->{pi}\t$read2_hits_ref->{$subject_id}->{bs}\t$read2_hits_ref->{$subject_id}->{ss}\t$read2_hits_ref->{$subject_id}->{se}\t$read2_hits_ref->{$subject_id}->{hs}";

    ${$unambig_readpair_reads_ref} += 2;
    ${$gt_unambig_readpair_reads_ref} += 2;

    return();

  }elsif($read_pair_subject_count > 1){

    #Multiple subjects found that have quality hits for both reads.  See if there is a clear winner according to combined score

    my $read1_hits_ref = $read_pair_blast_results{$read1_id}{blast_hits};
    my $read2_hits_ref = $read_pair_blast_results{$read2_id}{blast_hits};

    if ($ranked_subject_pair_hits{1}{combined_score} > $ranked_subject_pair_hits{2}{combined_score}){

      my $subject_id = $ranked_subject_pair_hits{1}{subject_id};

      unless ($read1_hits_ref->{$subject_id} && $read2_hits_ref->{$subject_id}){
	print RED, "\nMultiple Paired Reads - Unambiguous.  Subject ID: $subject_id\tREAD1: $read1_id\tREAD2: $read2_id", RESET;
	print Dumper $read1_hits_ref;
	print Dumper $read2_hits_ref;
	exit();
      }

      $top_hits_db_ref->{$read1_id} = "Top_Hit\t$subject_id\t$read1_hits_ref->{$subject_id}->{al}\t$read1_hits_ref->{$subject_id}->{pi}\t$read1_hits_ref->{$subject_id}->{bs}\t$read1_hits_ref->{$subject_id}->{ss}\t$read1_hits_ref->{$subject_id}->{se}\t$read1_hits_ref->{$subject_id}->{hs}";
      $top_hits_db_ref->{$read2_id} = "Top_Hit\t$subject_id\t$read2_hits_ref->{$subject_id}->{al}\t$read2_hits_ref->{$subject_id}->{pi}\t$read2_hits_ref->{$subject_id}->{bs}\t$read2_hits_ref->{$subject_id}->{ss}\t$read2_hits_ref->{$subject_id}->{se}\t$read2_hits_ref->{$subject_id}->{hs}";

      ${$unambig_readpair_reads_ref} += 2;
      ${$gt_unambig_readpair_reads_ref} += 2;

      return();
    }else{

      my $subject_id = $ranked_subject_pair_hits{1}{subject_id};

      unless ($read1_hits_ref->{$subject_id} && $read2_hits_ref->{$subject_id}){
	print RED, "\nPaired Read - Ambiguous.  Subject ID: $subject_id\tREAD1: $read1_id\tREAD2: $read2_id", RESET;
	print Dumper $read1_hits_ref;
	print Dumper $read2_hits_ref;
	exit();
      }

      $top_hits_db_ref->{$read1_id} = "Ambiguous\tNA\t$read1_hits_ref->{$subject_id}->{al}\t$read1_hits_ref->{$subject_id}->{pi}\t$read1_hits_ref->{$subject_id}->{bs}\tNA\tNA\tNA";
      $top_hits_db_ref->{$read2_id} = "Ambiguous\tNA\t$read2_hits_ref->{$subject_id}->{al}\t$read2_hits_ref->{$subject_id}->{pi}\t$read2_hits_ref->{$subject_id}->{bs}\tNA\tNA\tNA";

      ${$ambig_readpair_reads_ref} += 2;
      ${$gt_ambig_readpair_reads_ref} += 2;

      return();
    }
  }

  #C.) There are no subjects where both reads map to a single subject.  OR the BIT SCORE of one or both of these reads is below the cutoff.
  #   - Treat these as individual reads.  Use the bit score to identify the 'top hit'.  If there are multiple equal top hits, mark the read as Ambiguous
  #   - Once both reads have been processed, return().

  #print YELLOW, "\n\tNo Single Subject for Read Pair - Treating reads individually", RESET;

  foreach my $working_read_id (sort {$read_pair_blast_results{$a}{read_num} <=> $read_pair_blast_results{$b}{read_num}} keys %read_pair_blast_results){

    my $hits_ref = $read_pair_blast_results{$working_read_id}{blast_hits};
    my $hit_count = keys %{$hits_ref};

    #Create a sorted list to help determine the ranking of hits to each subject according to Bit Scores
    my %sorted_list;
    foreach my $subject_id (sort keys %{$hits_ref}){
      my $bs = $hits_ref->{$subject_id}->{bs};

      if ($sorted_list{$bs}){

	my $list_ref = $sorted_list{$bs}{subject_list};
	push (@{$list_ref}, $subject_id);
      }else{
	my @subject_list;
	push (@subject_list, $subject_id);
	$sorted_list{$bs}{subject_list} = \@subject_list;
      }
    }

  BS:foreach my $bs (sort {$b <=> $a} keys %sorted_list){
      my @subject_list = @{$sorted_list{$bs}{subject_list}};

      #D.) Print the best hits for individual reads to the 'top_hits' file

      my $top_subject_count = scalar(@subject_list);
      my $subject_id = $subject_list[0];

      if ($top_subject_count == 1){

	#DEBUG
	unless ($hits_ref->{$subject_id}){
	  print RED, "\nSingle Read - Unambiguous.  Subject ID: $subject_id\tREAD: $working_read_id", RESET;
	  print Dumper $hits_ref;
	  exit();
	}

	$top_hits_db_ref->{$working_read_id} = "Top_Hit\t$subject_id\t$hits_ref->{$subject_id}->{al}\t$hits_ref->{$subject_id}->{pi}\t$hits_ref->{$subject_id}->{bs}\t$hits_ref->{$subject_id}->{ss}\t$hits_ref->{$subject_id}->{se}\t$hits_ref->{$subject_id}->{hs}";

	${$unambig_individual_reads_ref}++;
	${$gt_unambig_individual_reads_ref}++;
	last (BS);

      }else{

	unless ($hits_ref->{$subject_id}){
	  print RED, "\nSingle Read - Ambiguous.  Subject ID: $subject_id\tREAD: $working_read_id", RESET;
	  print Dumper $hits_ref;
	  exit();
	}

	$top_hits_db_ref->{$working_read_id} = "Ambiguous\tNA\t$hits_ref->{$subject_id}->{al}\t$hits_ref->{$subject_id}->{pi}\t$hits_ref->{$subject_id}->{bs}\tNA\tNA\tNA";

	${$ambig_individual_reads_ref}++;
	${$gt_ambig_individual_reads_ref}++;
	last (BS);

      }
    }
  }
  return();
}



=head2 joinRecords

=over 3

=item Function:

Join top_hit records to read record summary files for a particular type of sequence mapping (e.g. introns) and update the master read record file to assign reads with quality hits to that type

=item Return:

NULL

=item Args:

'-subjects'   - Reference to map of subject sequences (e.g. introns)

'-seq_type'   - Sequence type of subject sequences (e.g. 'Intron')

'-upper_bit_score'  - Bit Score used to assign reads to the specified seq type

'-top_hits_db' - Reference to Berkeley DB object

'-hit_u_count' - Grand count of unambiguous hits of this type

'-hit_a_count' - Grand count of ambiguous hits of this type

'-read_records_infile'  - top_hits temp file

'-summary_outfile'   - summary read records file specific to the type of alignment parsed (e.g. intron alignments)

'-library_name'   - name of library currently being processed

'-working_dir'   - directory containing temp files

'-log_file_handle' - Handle for LOG file

=item Example(s):

&joinRecords('-subjects'=>\%subjects, '-seq_type'=>$seq_type, '-upper_bit_score'=>$upper_bit_score, '-top_hits_db'=>\%top_hits_db, '-hit_u_count'=>\$hit_u_count, '-hit_a_count'=>\$hit_a_count,
             '-read_records_file'=>$read_records_file, '-summary_outfile'=>$summary_outfile, '-library_name'=>$library_name, '-working_dir'=>$working_dir, '-min_align_length'=>$min_align_length,
             '-log_file_handle'=>$log_fh);


=back

=cut

sub joinRecords{
  my %args = @_;
  my $subjects_ref = $args{'-subjects'};
  my $seq_type = $args{'-seq_type'};
  my $upper_bit_score = $args{'-upper_bit_score'};
  my $min_align_length = $args{'-min_align_length'};
  my $top_hits_db_ref = $args{'-top_hits_db'};
  my $hit_u_count = $args{'-hit_u_count'};
  my $hit_a_count = $args{'-hit_a_count'};
  my $read_records_file = $args{'-read_records_file'};
  my $summary_outfile = $args{'-summary_outfile'};
  my $library_name = $args{'-library_name'};
  my $working_dir = $args{'-working_dir'};
  my $log_fh = $args{'-log_file_handle'};

  my $base;
  if ($read_records_file =~ /(.*)\.gz$/){
    $base = $1;
  }else{
    print RED, "\nFormat of read records infile name not understood: $read_records_file - not compressed?\n\n", RESET;
    exit();
  }

  my %read_pairs_db;
  my $read_pairs_db_file = "$working_dir"."$library_name"."_ReadPairs.btree";
  system ("rm -f $read_pairs_db_file");
  tie(%read_pairs_db, 'BerkeleyDB::Btree', -Cachesize => 256000000, -Filename=> $read_pairs_db_file, -Flags => DB_CREATE) or die "can't open file $read_pairs_db_file: $! $BerkeleyDB::Error\n";

  print BLUE, "\n\nCreating Berkley DB of read records file.\n", RESET;
  print $log_fh "\n\nCreating Berkley DB of read records file.\n";

  my $r1_subject_name = "R1_"."$seq_type"."Name";
  my $r2_subject_name = "R2_"."$seq_type"."Name";

  my $seq_type_uc = uc($seq_type);
  my $status_u = "$seq_type_uc"."_U";
  my $status_a = "$seq_type_uc"."_A";

  #first create a berkley DB hash of the read records file
  my $counter = 0;
  my $header = 1;
  open (READS, "zcat $read_records_file |") || die "\nCould not open read records file: $read_records_file\n\n";
  while(<READS>){
    $counter++;
    chomp($_);
    if ($header == 1){
      $header = 0;
      next();
    }

    if ($counter == 100000){
      $| = 1; print BLUE, ".", RESET;  $| = 0;
      $counter = 0;
    }

    my @line = split("\t", $_);
    my $read_pair_id = $line[0];

    my $read1_status = $line[3];
    my $read2_status = $line[4];

    #Format: key{read_pair_id} value{read1_status\tread2_status}
    $read_pairs_db{$read_pair_id} = "$read1_status\t$read2_status";
  }
  close(READS);
  my $message = &memoryUsage();
  print YELLOW, "\n\n$message", RESET;
  print $log_fh "\n\n$message";

  #Now go through every read_pair_id, join with read1 and read2 BLAST hit info, and print result to outfile
  print BLUE, "\n\nCreate summary file: $summary_outfile\n", RESET;
  print $log_fh "\n\nCreate summary file: $summary_outfile\n";

  open (SUMMARY, ">$summary_outfile") || die "\nCould not open summary outfile: $summary_outfile\n\n";

  print SUMMARY "Read_ID\tDistanceBetweenReads\tR1_ID\tR1_HitType\t$r1_subject_name\tR1_Strand\tR1_AlignmentLength\tR1_PercentIdentity\tR1_BitScore\tR1_Chr\tR1_ChrStart\tR1_ChrEnd\tR2_ID\tR2_HitType\t$r2_subject_name\tR2_Strand\tR2_AlignmentLength\tR2_PercentIdentity\tR2_BitScore\tR2_Chr\tR2_ChrStart\tR2_ChrEnd\n";

  $counter = 0;
  my $read_count = 0;
  my $read_pairs_to_be_printed = 0;

  while (my ($read_pair_id, $read_pair_string) = each %read_pairs_db){
    $read_count+=2;
    $counter++;
    if ($counter == 100000){
       $| = 1; print BLUE, ".", RESET;  $| = 0;
       $counter = 0;
    }

    my $read1_id = "$read_pair_id"."_R1";
    my $read2_id = "$read_pair_id"."_R2";

    my @read_pair_string = split("\t", $read_pair_string);
    my $read1_status = $read_pair_string[0];
    my $read2_status = $read_pair_string[1];

    #Watch for undefined values that might indicate data corruption of simple mis-firing of the script for some I/O related insanity
    unless(defined($read1_status) && defined($read2_status)){
      print RED, "\n\nA read status was undefined - aborting\n\n", RESET;
      exit();
    }

    #make sure there was a hit for at least one of the two reads
    unless ($top_hits_db_ref->{$read1_id} || $top_hits_db_ref->{$read2_id}){
      next();
    }
    $read_pairs_to_be_printed++;

    #Format of blast hit record stored:
    #Top_Hit\t$subject_id\t$hits_ref->{$subject_id}->{al}\t$hits_ref->{$subject_id}->{pi}\t$hits_ref->{$subject_id}->{bs}\t$hits_ref->{$subject_id}->{ss}\t$hits_ref->{$subject_id}->{se}\t$hits_ref->{$subject_id}->{hs}


    #Determine the fragment size based on the alignment of the two reads.  Only consider cases where both reads have 'Top_Hits'
    #Use the outer distance of the coords
    my $distance = "NA";

    if ($top_hits_db_ref->{$read1_id} && $top_hits_db_ref->{$read2_id}){
      my $string1 = $top_hits_db_ref->{$read1_id};
      my @string1 = split("\t", $string1);
      my $string2 = $top_hits_db_ref->{$read2_id};
      my @string2 = split("\t", $string2);

      my $type1 = $string1[0];
      my $type2 = $string2[0];
      my $subject_id1 = $string1[1];
      my $subject_id2 = $string2[1];
      my $start1 = $string1[5];
      my $start2 = $string2[5];
      my $end1 = $string1[6];
      my $end2 = $string2[6];
      my @coords = ($start1, $start2, $end1, $end2);

      if ($type1 eq "Top_Hit" && $type2 eq "Top_Hit"){
        my @coords_sort = sort{$a <=> $b}(@coords);

        #Calculate the distance between Read1 and Read2 in the subject sequence
	if ($subject_id1 eq $subject_id2){
	  $distance = $coords_sort[3] - $coords_sort[0];
	}else{
	  $distance = "$seq_type"."_Mismatch";
	}
      }
    }

    print SUMMARY "$read_pair_id\t$distance";

    #READ 1#
    if ($top_hits_db_ref->{$read1_id}){
      my $string = $top_hits_db_ref->{$read1_id};
      my @string = split("\t", $string);

      my $hit_type = $string[0];
      my $subject_id = $string[1];
      my $align_length = $string[2];
      my $percent_id = $string[3];
      my $bit_score = $string[4];
      my $start = $string[5];
      my $end = $string[6];
      my $strand = $string[7];

      #Convert coordinates
      my $chr_start = "NA";
      my $chr_end = "NA";
      my $chr = "NA";
      if ($hit_type eq "Top_Hit"){      
        my $source_strand = $subjects_ref->{$subject_id}->{strand};
        my $subject_start = 1;
        my $subject_end = ($subjects_ref->{$subject_id}->{end} - $subjects_ref->{$subject_id}->{start});
        my $subject_start_chr = $subjects_ref->{$subject_id}->{start_chr};
        my $subject_end_chr = $subjects_ref->{$subject_id}->{end_chr};
        my $coords_ref = &convertGeneCoordinatesGeneric ('-subject_id'=>$subject_id, '-source_strand'=>$source_strand, 
                                                         '-subject_start'=>$subject_start, '-subject_end'=>$subject_end, '-subject_start_chr'=>$subject_start_chr, '-subject_end_chr'=>$subject_end_chr,
                                                         '-query_start'=>$start, '-query_end'=>$end, '-ordered'=>"yes");

        $chr = $subjects_ref->{$subject_id}->{chromosome};
        $chr_start = $coords_ref->{$subject_id}->{chr_start};
        $chr_end = $coords_ref->{$subject_id}->{chr_end};
      }

      #print SUMMARY "\tR1_ID\tR1_HitType\tR1_SubjectID\tR1_Strand\tR1_AlignmentLength\tR1_PercentIdentity\tR1_BitScore\tR1_ChrStart\R1_ChrEnd";
      print SUMMARY "\t$read1_id\t$hit_type\t$subject_id\t$strand\t$align_length\t$percent_id\t$bit_score\t$chr\t$chr_start\t$chr_end";

      #If the upper bit score criteria is met and the status of this read is still 'Unassigned' in the read record file, change the status
      #Possible status is '$status_u' or '$status_a'

      #If a read has already been assigned to '$status_u' or '$status_a', start by resetting these to 'Unassigned'
      if ($read1_status =~ /$status_u/i || $read1_status =~ /$status_a/i){
	$read1_status = "Unassigned";
      }

      #Make sure that the status of this reads is still 'Unassigned'
      if ($read1_status eq "Unassigned"){

	#Now, if the top hit for this read meets the user specified threshold, change its status
	#- The resulting status will depend on whether the hit was ambiguous or not ...
	if (($hit_type eq "Top_Hit") && ($bit_score >= $upper_bit_score) && ($align_length >= $min_align_length)){
	  ${$hit_u_count}++;
	  $read1_status = $status_u;
	}elsif(($hit_type eq "Ambiguous") && ($bit_score >= $upper_bit_score) && ($align_length >= $min_align_length)){
	  ${$hit_a_count}++;
	  $read1_status = $status_a;
	}
      }
    }else{
      #Read 1 is not defined but READ 2 was
      print SUMMARY "\t$read1_id\tNone\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA";
    }

    #READ 2#
    if ($top_hits_db_ref->{$read2_id}){
      my $string = $top_hits_db_ref->{$read2_id};
      my @string = split("\t", $string);

      my $hit_type = $string[0];
      my $subject_id = $string[1];
      my $align_length = $string[2];
      my $percent_id = $string[3];
      my $bit_score = $string[4];
      my $start = $string[5];
      my $end = $string[6];
      my $strand = $string[7];


      #Convert coordinates
      my $chr_start = "NA";
      my $chr_end = "NA";
      my $chr= "NA";
      if ($hit_type eq "Top_Hit"){      
        my $source_strand = $subjects_ref->{$subject_id}->{strand};
        my $subject_start = 1;
        my $subject_end = ($subjects_ref->{$subject_id}->{end} - $subjects_ref->{$subject_id}->{start});
        my $subject_start_chr = $subjects_ref->{$subject_id}->{start_chr};
        my $subject_end_chr = $subjects_ref->{$subject_id}->{end_chr};
        my $coords_ref = &convertGeneCoordinatesGeneric ('-subject_id'=>$subject_id, '-source_strand'=>$source_strand, 
                                                         '-subject_start'=>$subject_start, '-subject_end'=>$subject_end, '-subject_start_chr'=>$subject_start_chr, '-subject_end_chr'=>$subject_end_chr,
                                                         '-query_start'=>$start, '-query_end'=>$end, '-ordered'=>"yes");

        $chr = $subjects_ref->{$subject_id}->{chromosome};
        $chr_start = $coords_ref->{$subject_id}->{chr_start};
        $chr_end = $coords_ref->{$subject_id}->{chr_end};
      }

      #print SUMMARY "\tR2_ID\tR2_HitType\tR2_SubjectID\tR2_AlignmentLength\tR2_PercentIdentity\tR2_BitScore\tR2_ChrStart\tR2_Chr_End\n";
      print SUMMARY "\t$read2_id\t$hit_type\t$subject_id\t$strand\t$align_length\t$percent_id\t$bit_score\t$chr\t$chr_start\t$chr_end\n";

      #If the upper bit score criteria is met and the status of this read is still 'Unassigned' in the read record file, change the status
      #Possible status is '$status_u' or '$status_a'

      #If a read has already been assigned to '$status_u' or '$status_a', start by resetting these to 'Unassigned'
      if ($read2_status =~ /$status_u/i || $read2_status =~ /$status_a/i){
	$read2_status = "Unassigned";
      }

      #Make sure that the status of this reads is still 'Unassigned'
      if ($read2_status eq "Unassigned"){

	#Now, if the top hit for this read meets the user specified threshold, change its status
	#- The resulting status will depend on whether the hit was ambiguous or not ...
	if (($hit_type eq "Top_Hit") && ($bit_score >= $upper_bit_score) && ($align_length >= $min_align_length)){
	  ${$hit_u_count}++;
	  $read2_status = $status_u;
	}elsif(($hit_type eq "Ambiguous") && ($bit_score >= $upper_bit_score) && ($align_length >= $min_align_length)){
	  ${$hit_a_count}++;
	  $read2_status = $status_a;
	}
      }
    }else{
      #Read 2 is not defined but Read 1 was
      print SUMMARY "\t$read2_id\tNone\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n";

    }
    $read_pairs_db{$read_pair_id} = "$read1_status\t$read2_status";
  }
  close(SUMMARY);
  $message = &memoryUsage();
  print YELLOW, "\n\n$message", RESET;
  print $log_fh "\n\n$message";
  
  print BLUE, "\n\nFound hits for one or both reads of $read_pairs_to_be_printed read pairs\n", RESET;

  print BLUE, "\n\nCompress summary outfile: $summary_outfile\n", RESET;
  my $compress_cmd = "gzip -f $summary_outfile";
  system($compress_cmd);


  #Finally go through the original read records file and update the read status values where appropriate
  print BLUE, "\n\nUpdate read status values in read record file: $read_records_file\n", RESET;
  print $log_fh "\n\nUpdate read status values in read record file: $read_records_file\n";

  my $new_read_record_file = "$base"."2";
  $header = 0;
  open (READS, "zcat $read_records_file |") || die "\nCould not open read records file: $read_records_file\n\n";
  open (NEW_READS, ">$new_read_record_file") || die "\nCould not open new read records file: $new_read_record_file\n\n";

  $counter= 0;
  while(<READS>){
    $counter++;
    chomp($_);

    if ($header == 0){
      $header = 1;
      print NEW_READS "$_\n";
      next();
    }

    if ($counter == 100000){
       $| = 1; print BLUE, ".", RESET;  $| = 0;
       $counter = 0;
    }

    #Values from $line[0] to $line[10];
    my @line = split("\t", $_);

    #Make sure all values are defined!!  Missing values might indicate file corruption at some point
    unless ($line[0] && $line[1] && $line[2] && $line[3] && $line[4] && $line[5] && $line[6] && ($line[7] =~ /\d+/) && ($line[8] =~ /\d+/) && ($line[9] =~ /\d+/) && ($line[10] =~ /\d+/)){
      print RED, "\n\nFound an undefined value in the input read records file!!!\n\n", RESET;
      print RED, "RECORD: $_\n\n", RESET;
      exit();
    }

    my $read_id = $line[0];

    my $string = $read_pairs_db{$read_id};
    my @string = split("\t", $string);

    my $read1_status = $string[0];
    my $read2_status = $string[1];

    print NEW_READS "$line[0]\t$line[1]\t$line[2]\t$read1_status\t$read2_status\t$line[5]\t$line[6]\t$line[7]\t$line[8]\t$line[9]\t$line[10]\n";
  }
  close(READS);
  close(NEW_READS);
  $message = &memoryUsage();
  print YELLOW, "\n\n$message", RESET;
  print $log_fh "\n\n$message";

  untie(%read_pairs_db);

  my $rm_cmd2 = "rm -f $read_pairs_db_file";
  system($rm_cmd2);

  my $percent_u = sprintf("%.2f", ((${$hit_u_count}/$read_count)*100));
  my $percent_a = sprintf("%.2f", ((${$hit_a_count}/$read_count)*100));

  print BLUE, "\n\nAssigned ${$hit_u_count} reads ($percent_u%) to the status '$status_u', and ${$hit_a_count} ($percent_a%) reads to the status '$status_a'\n\n", RESET;
  print $log_fh "\n\nAssigned ${$hit_u_count} reads ($percent_u%) to the status '$status_u', and ${$hit_a_count} ($percent_a%) reads to the status '$status_a'\n\n";

  my $cmd_gzip = "gzip -f $new_read_record_file";
  system ($cmd_gzip);

  my $cmd_mv = "mv $new_read_record_file".".gz"." $read_records_file";
  print BLUE, "\n\nOverwriting existing read records file:\n\t$cmd_mv\n\n", RESET;
  print $log_fh "\n\nOverwriting existing read records file:\n\t$cmd_mv\n\n";

  system($cmd_mv);

  return();
}


=head2 convertGeneCoordinatesGeneric()

=over 3

=item Function:

Convert coordinates that are relative to the a subsequence of a chromosome (could be an intron, exon, or arbitrary region) to chromosome coordinates

Say you have mapped reads within a sub-sequence of a chromosome (e.g. an intron).  If you know the size of this region and its chromosome coordinates ...

You can calculate the chromosome coordinates of any sub-sequence within this sub-sequence

NOTE: this script assumes that the sequence objects in use have been oriented to the +ve strand (i.e. things on the -ve strand have been reverse-complemented.) - so source strand has to be accounted for!

=item Return:

Coordinate object (keyed on subject id) as a hash - chromosome start pos, and chromosome end pos

=item Args:

'-subject_id' => ID of intron, etc. for which a coordinate conversion is required

'-source_strand'  => Source chromosomal strand of the object

'-subject_start' => Subject based start position (should always be 1)

'-subject_end' => Subject based end position (should be equal to the length of the subject)

'-subject_start_chr' => Chromosome coordinate of start position of subject sequence

'-subject_end_chr' => Chromosome coordinate of end position of subject sequence

'-query_start' => Subject based start coordinate to be converted to chromosome coordinates

'-query_end' => Subject based end coordinate to be converted to chromosome coordinates

'-ordered' => "yes|no"  #If 'yes', ensure that start coord is always smaller than end regardless of source strand

=item Example(s):

my %coords = %{&convertGeneCoordinatesGeneric ('-subject_id'=>$intron_id, '-source_strand'=>$source_strand, 
                                               '-subject_start'=>$intron_start, '-subject_end'=>$intron_end, '-subject_start_chr'=>$intron_start_chr, '-subject_end_chr'=>$intron_end_chr,
                                               '-query_start'=>$start, '-query_end'=>$end, '-ordered'=>"yes")};

my $chr_start = $coords{$subject_id}{chr_start};

my $chr_end = $coords{$subject_id}{chr_end};

=back

=cut


###################################
#convertGeneCoordinatesGeneric    #
###################################
sub convertGeneCoordinatesGeneric{

  my %args = @_;
  my $subject_id = $args{'-subject_id'};
  my $chr_strand = $args{'-source_strand'};
  my $subject_start = $args{'-subject_start'};
  my $subject_end = $args{'-subject_end'};
  my $subject_start_chr = $args{'-subject_start_chr'};
  my $subject_end_chr = $args{'-subject_end_chr'};
  my $query_start = $args{'-query_start'};
  my $query_end = $args{'-query_end'};
  my $ordered = $args{'-ordered'};

  my %coords;

  #Check input
  unless ($subject_id && ($chr_strand eq "1" || $chr_strand eq "-1") && ($subject_start =~ /\d+/) && ($subject_end =~ /\d+/) && ($subject_start_chr =~ /\d+/) && ($subject_end_chr =~ /\d+/) && ($query_start =~ /\d+/) && ($query_end =~ /\d+/) && ($ordered =~ /yes|y|no|n/i)){
    print Dumper %args;
    print YELLOW, "subject_id: $subject_id\tchr_strand: $chr_strand\tsubject_start: $subject_start\tsubject_end: $subject_end\tsubject_start_chr: $subject_start_chr\tsubject_end_chr: $subject_end_chr\tquery_start: $query_start\tquery_end: $query_end\tordered: $ordered", RESET;
    print RED, "One of the parameters of convertGeneCoordinatesGeneric() is missing or incorrect format!", RESET;
    exit();
  }

  #Note: All gene/exon/intron coordinates stored in ALEXA are relative to the coding strand (i.e. start codon near beginning, end codon near the end)
  #This means that when coverting exon, or other coordinates back to the chromosome context, the original 'strand' of the gene on the chromosome
  #needs to be considered.

  #Make sure the supplied coordinates are actually within the specified gene
  unless ($query_start >= $subject_start-1 && $query_start <= $subject_end+1){
    print RED, "\nQuery Start coordinate ($query_start) provided to convertGeneCoordinatesGeneric() is not within the specified Subject coordinates\n\n", RESET;
    exit();
  }
  unless ($query_end >= $subject_start-1 && $query_end <= $subject_end+1){
    print RED, "\nQuery End coordinate ($query_end) provided to convertGeneCoordinatesGeneric() is not within the specified Subject coordinates\n\n", RESET;
    exit();
  }

  #Convert provided subject coordinates to coordinates relative to the chromosome
  if ($chr_strand == 1){
    my $query_chr_start = $subject_start_chr + $query_start - 1;
    my $query_chr_end = $subject_start_chr + $query_end - 1;

    if ($ordered =~ /yes|y/i){
      #Make sure the start and end are reported such that start is always smaller than end
      my $temp;
      if ($query_chr_start > $query_chr_end){
	$temp = $query_chr_start;
	$query_chr_start = $query_chr_end;
	$query_chr_end = $temp;
      }
    }
    $coords{$subject_id}{chr_start} = $query_chr_start;
    $coords{$subject_id}{chr_end} = $query_chr_end;

  }elsif ($chr_strand == -1){

    my $query_chr_start = $subject_end_chr - $query_end + 1;
    my $query_chr_end = $subject_end_chr - $query_start + 1;

    if ($ordered =~ /yes|y/i){
      #Make sure the start and end are reported such that start is always smaller than end
      my $temp;
      if ($query_chr_start > $query_chr_end){
	$temp = $query_chr_start;
	$query_chr_start = $query_chr_end;
	$query_chr_end = $temp;
      }
    }

    $coords{$subject_id}{chr_start} = $query_chr_start;
    $coords{$subject_id}{chr_end} = $query_chr_end;

  }else{
    print "\nStrand format: $chr_strand not understood by convertGeneCoordinates()!\n\n";
    exit();
  }

  return(\%coords);
}


=head2 removeFailingReads()

=over 3

=item Function:

Given a directory containing Read Records file, populate a hash with only those Read IDs that pass user specified criteria

Read record files are assumed to compressed

The purpose of this function is to select only those reads that have been assigned to a particular read class after parsing...

Ambiguous hits will be excluded, as well as hits that did not reach the bit_score cutoff as well as reads

Similarly a read that maps to an Intron will not be allowed if it has been already been assigned to an EnsEMBL transcript (prevent a read from being used more than once for expression estimates)

=item Return:

Total number of unambiguous reads mapped in the library (all classes combined)

=item Args:

'-read_record_dir' => Path to directory containing read record files

'-reads_object' => reference to hash to be populated (may be a hash tied to a berkeley DB file)

'-read_class' => class of read which will allow it to be added to the passing read list (e.g. ENST_U, NOVEL_JUNCTION_U, NOVEL_BOUNDARY_U, INTRON_U, INTERGENIC_U)

=item Example(s):

my $pass_read_count = &importPassingReads ('-read_record_dir'=>$read_record_dir, '-reads_object'=>\%passing_reads, '-read_class'=>$read_class);

=back

=cut


###################################
#removeFailingReads               #
###################################
sub removeFailingReads{
  my %args = @_;
  my $read_record_dir = $args{'-read_record_dir'};
  my $reads_ref = $args{'-reads_object'};
  my $read_class = $args{'-read_class'};

  my $library_size = 0; #All reads - any class
  my $library_bases = 0; #All bases - any class

  my $u_class = "ENST_U INTRON_U INTERGENIC_U NOVEL_JUNCTION_U NOVEL_BOUNDARY_U"; #Classes that qualify as a mapped read
  #my $mapped_reads = 0; #Reads that have been assigned to one of these classes
  my $mapped_bases = 0;

  #NOTE:  Instead of mapped reads, this code should probably be changed to use mapped bases (to account for libraries that use different reads lengths)
  #       - Now if one library has longer reads than a comparison library, it will tend to have higher expression values and this will Skew the DE distribution
  #       - Since reads are required to map over most of their length, perhaps use read length for assigned reads to determined the number of mapped bases

  my $start_pass_count = keys %{$reads_ref};
  my $end_pass_count = 0;
  print BLUE, "\n\nStarting with a total of $start_pass_count reads - checking their status\n", RESET;

  #If there are no starting reads there is no point in continuing
  if ($start_pass_count == 0){
    print YELLOW, "\n\nThere are no reads that need to be checked for status - returning\n\n", RESET;
    return($end_pass_count);
  }

  #Get files from this directory
  print BLUE, "\n\nSearching $read_record_dir for read records files", RESET;

  my %files;
  opendir(DIRHANDLE, "$read_record_dir") || die "\nCannot open directory: $read_record_dir\n\n";
  my @test_files = readdir(DIRHANDLE);
  my $file_count = 0;

  foreach my $test_file (sort @test_files){
    my $file_path = "$read_record_dir"."$test_file";
    my $first_read_pair_id;
    my $test = 1;

    #Skip directories within the specified directory
    if (-e $file_path && -d $file_path){
      print YELLOW, "\n\t$file_path  is a directory - skipping", RESET;
      next();
    }

    #If the results file is compressed uncompress it
    unless ($file_path =~ /(.*)\.gz$/){
      print RED, "\nFound an uncompressed file: $file_path - make sure no files in progress are in this directory!!\n\n", RESET;
      exit();
    }
    $file_count++;
    print BLUE, "\n\t$file_path was added to the list of files to be processed", RESET;
    $files{$file_count}{path} = $file_path;
  }
 
  print BLUE, "\n\nProcessing these files to find reads with a status of: $read_class (all chromosomes)", RESET;
  foreach my $fc (sort {$a <=> $b} keys %files){
    my $file_path = $files{$fc}{path};
    my %columns;

    print BLUE, "\n\tProcessing: $file_path", RESET;

    open (READ, "zcat $file_path |") || die "\nCould not open file: $file_path\n\n";

    my $header = 1;
    while(<READ>){
      chomp($_);
      my @line = split("\t", $_);

      if ($header == 1){
        my $column_count = 0;
        foreach my $column (@line){
          $columns{$column}{column_pos} = $column_count;
          $column_count++;
        }
        $header = 0;
        next();
      }
      my $r1_id = $line[$columns{'Read1_ID'}{column_pos}];
      my $r2_id = $line[$columns{'Read2_ID'}{column_pos}];
      my $r1_status = $line[$columns{'Read1_Status'}{column_pos}];
      my $r2_status = $line[$columns{'Read2_Status'}{column_pos}];
      my $r1_seq_length = length($line[$columns{'Read1_Seq'}{column_pos}]);
      my $r2_seq_length = length($line[$columns{'Read2_seq'}{column_pos}]);
      $library_size++;
      $library_bases += $r1_seq_length;
      $library_bases += $r2_seq_length;

      #If read was mapped unambiguously to any class (except repeat), add it to the grand library
      if ($u_class =~ /($r1_status)/){
        #$mapped_reads++;
        $mapped_bases += $r1_seq_length;
      }
      if ($u_class =~ /($r2_status)/){
        #$mapped_reads++;
        $mapped_bases += $r2_seq_length;
      }

      if ($reads_ref->{$r1_id}){
        if ($read_class =~ /($r1_status)/){
          #Do nothing
        }else{
          delete($reads_ref->{$r1_id});
        }
      }
      if ($reads_ref->{$r2_id}){
        if ($read_class =~ /($r2_status)/){
          #Do nothing
        }else{
          delete($reads_ref->{$r2_id});
        }
      }
    }
    close(READ);
  }
  
  $end_pass_count = keys %{$reads_ref};
  my $pass_count_p = sprintf("%.2f",(($end_pass_count/$start_pass_count)*100));
  #print BLUE, "\n\nFound a grand total of $library_size reads for this library - A total of $mapped_reads reads could be mapped to a transcript, junction, boundary, intron or intergenic region", RESET;
  print BLUE, "\n\nFound a grand total of $library_size reads and $library_bases bases for this library - A total of $mapped_bases bases could be mapped to a transcript, junction, boundary, intron or intergenic region", RESET;
  print BLUE, "\nAfter deleting those without the specified class: $read_class a total of $end_pass_count ($pass_count_p%) reads remain\n", RESET;
  
  return($mapped_bases);
}


=head2 importExpressionCutoffs()

=over 3

=item Function:

Given a file containing pre-calculated gene-by-gene expression cutoffs, create a hash to represent this data

Should find one record for every gene, plus a false gene id '0' which represents the cutoffs for intergenic regions

=item Return:

Hash keyed on ALEXA gene id with expression cutoff value for each gene

=item Args:

'-cutoffs_file' => Path to file containining cutoff data for each gene

=item Example(s):

my $gene_cutoffs_ref = &importExpressionCutoffs ('-cutoffs_file'=>$cutoffs_file);

=back

=cut

###################################
#importExpressionCutoffs          #
###################################
sub importExpressionCutoffs{
  my %args = @_;
  my $cutoffs_file = $args{'-cutoffs_file'};

  unless(-e $cutoffs_file){
    print RED, "\nCutoffs file supplied to importExpressionCutoffs() does not appear valid: $cutoffs_file\n\n", RESET;
    exit();
  }

  my %gene_cutoffs;

  open(CUTOFFS, "$cutoffs_file") || die "\nCould not open cutoffs file: $cutoffs_file\n\n"; 

  my $header = 1;
  while(<CUTOFFS>){
    chomp($_);
    my @line = split("\t", $_);
    if ($header == 1){
      $header = 0;
      next();
    }
    my $gene_id = $line[0];
    my $average_coverage = $line[1];
    my $cutoff = $line[3];
    $gene_cutoffs{$gene_id}{cutoff} = $cutoff;
    $gene_cutoffs{$gene_id}{average_coverage} = $average_coverage;
  }
  close(CUTOFFS);

  my $cutoffs_found = keys %gene_cutoffs;

  print BLUE, "\n\nFound $cutoffs_found gene-by-gene cutoffs (including 1 intergenic cutoff)\n\n", RESET;

  return(\%gene_cutoffs);
}


=head2 testExpression()

=over 3

=item Function:

Given a gene_cutoffs object, a gene id (or list of gene IDs) and the observed expression level, decide whether an element (gene, transcript, exon, etc.) is expressed

=item Return:

'1' for expressed, '0' for not expressed

=item Args:

'-cutoffs_file' => Path to file containining cutoff data for each gene

=item Example(s):

my $cutoff_test = &testExpression('-cutoffs_ref'=>$gene_cutoffs_ref, '-gene_id'=>$gene_id, '-expression_value'=>$normalized_expression_level)

=back

=cut

###################################
#testExpression                   #
###################################
sub testExpression{
  my %args = @_;
  my $cutoffs_ref = $args{'-cutoffs_ref'};
  my $gene_id = $args{'-gene_id'};
  my $norm_expression_value = $args{'-norm_expression_value'};
  my $raw_expression_value = $args{'-raw_expression_value'}; 
  my $percent_gene_expression_cutoff = $args{'-percent_gene_expression_cutoff'};

  my $expressed = 0;
  my $percent_gene_expression = 0;

  my @gene_list = split(" ", $gene_id);

  #NOTE:  We want to use the normalized (and corrected) expression values to test for expressed status

  if (scalar(@gene_list) == 1){

    #Make sure specified gene is defined in the gene_cutoffs object
    my $gene = $gene_list[0];
    unless($cutoffs_ref->{$gene}){
      print RED, "\nGene ID: $gene is not defined in the gene_cutoffs_ref supplied to testExpression()\n\n", RESET;
      exit();
    }

    if ($norm_expression_value > $cutoffs_ref->{$gene}->{cutoff}){
      $expressed = 1;
    }
    if ($raw_expression_value > 0){
      if ($cutoffs_ref->{$gene}->{average_coverage} =~ /\d+/){
        if ($cutoffs_ref->{$gene}->{average_coverage} > 0){
          my $x = ($raw_expression_value/$cutoffs_ref->{$gene}->{average_coverage})*100;
          $percent_gene_expression = sprintf("%.2f", $x);
        }
      }
    }
    #If a non-zero percent gene expression cutoff was supplied and this feature has a percent gene expression less than this cutoff, then the feature will not be expressed even if it's expression level is high enough
    if ($percent_gene_expression < $percent_gene_expression_cutoff){
      $expressed = 0;
    }

  }elsif(scalar(@gene_list > 1)){
    #Feature is associated with multiple genes - (i.e. introns of overlaping genes on opposite strands mostly)
    $expressed = 1;
    my $temp = 100000000000000000000000000000000000000000000000000;
    foreach my $gene (@gene_list){
      #Make sure specified gene is defined in the gene_cutoffs object
      unless($cutoffs_ref->{$gene}){
        print RED, "\nGene ID: $gene_id is not defined in the gene_cutoffs_ref supplied to testExpression()\n\n", RESET;
        exit();
      }
      #If the expression value is lower than the expression cutoff for any gene, this element will not be expressed
      #Similarly if the percent_gene_expression is lower than the specified precent_gene_expression_cutoff for any gene, this element will not be expressed
      if ($norm_expression_value <= $cutoffs_ref->{$gene}->{cutoff}){
        $expressed = 0;
      }
      if ($raw_expression_value > 0){
        if ($cutoffs_ref->{$gene}->{average_coverage} =~ /\d+/){
          if ($cutoffs_ref->{$gene}->{average_coverage} > 0){
            my $x = ($raw_expression_value/$cutoffs_ref->{$gene}->{average_coverage})*100;
            my $y = sprintf("%.2f", $x);
            if ($y < $temp){
              $temp = $y;
            }
          }
        }
      }
    }
    if ($temp <= 10000000000000000000000000000000000000000000000000){
      $percent_gene_expression = $temp;
    }
    if ($percent_gene_expression < $percent_gene_expression_cutoff){
      $expressed = 0;
    }
  }else{
    print RED, "\nGene ID/list supplied to testExpression() was not understood\n\n", RESET;
    exit();
  }

  my @results;
  push(@results, $expressed);
  push(@results, $percent_gene_expression);
  return(\@results);
}


=head2 importLibrarySize()

=over 3

=item Function:

Parse a MappedReads.txt file in the Summary dir of a Read Records Dir and get the library size (total mapped bases)

=item Return:

Integer value

=item Args:

'-read_records_dir' => Path to read records file for a particular library

=item Example(s):

my $library_size = &importLibrarySize ('-read_records_dir'=>$read_records_dir);

=back

=cut

###################################
#importLibrarySize                #
###################################
sub importLibrarySize{
  my %args = @_;
  my $read_records_dir = $args{'-read_records_dir'};

  my $mapped_reads_file = "$read_records_dir"."Summary/MappedReads.txt";
  unless(-e $mapped_reads_file){
    print RED, "\nSummary mapped_reads_file supplied to importLibrarySize() could noe be found!!  Have you calculated library sizes?  - Aborting\n\n", RESET;
    exit();
  }
  my $library_size;

  my %columns;
  open (IN, "$mapped_reads_file") || die "\n\nCould not open input file: $mapped_reads_file\n\n";
  my $header = 1;
  my $col_count = 0;
  while(<IN>){
    chomp($_);
    my @line = split("\t", $_);
    if ($header == 1){
      $header = 0;
      foreach my $value (@line){
        $columns{$value}{column_pos} = $col_count;
        $col_count++;
      }
      next();
    }

    #Watch for the grand total line in the file
    if ($line[0] eq "TOTAL"){
      $library_size = $line[$columns{'LIBRARY_SIZE_Bases'}{column_pos}];
    }
  }
  close(IN);

  return($library_size);
}



1;







