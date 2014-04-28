#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#The purpose of this script is to parse blast results files for Solexa Paired Reads that were BLASTed against a database of human Repeat Elements

#1.)  Parse blast results files one at a time
#1a.) Go through the specified blast results files (in the specified directory)
#    - For each blast results file go through the hits for a single read
#    - Store the best Repeat hit
#    - If the best and second best hit have the same bit score, consider this a tie and make note of it as an ambiguous repeat hit

#2.) Print the best hit to the 'top_hits' file
#   - Read_ID, Hit_Type (i.e. Ambiguous or Top Hit), Repeat_ID,  AlignLength, Percent Identity, Subject Start, Subject End

#3.) Join top hits results for each read of a pair to a single line of data and print as a summary of the blast results
#    - Count entries in the read records file.  Store all read IDs as an array.
#    - Note that not all read records will have a top hits record.  There may not have been any blast hits to either read of the pair
#    - Parse a BLOCK of full read records.  Create a hash and key on Read_ID
#    - Now parse the top hits file and build a second hash for all individual reads that correspond to those in the BLOCK
#    - Once these two hashes are built go through the read records hash and print out corresponding top hits to a new summary file
#      - If a particular read of a read pair in the read records hash does not have a top hits record, fill in values with 'NA'
#      - Each read record will have HitType value which is either the 'Top_Hit', 'Ambiguous', or 'None'
#    - Initialize both hashes and start processing the next block

#    - Columns to be added to the Repeats blast results summary file:
#    - Read_ID, DistanceBetweenReads, R1_ID, R1_HitType, R1_RepeatName[repeat|family|source], R1_AlignmentLength, R1_PercentIdentity, R1_BitScore, R2_ID, R2_HitType, R2_RepeatName[repeat|family|source], R2_AlignmentLength, R2_PercentIdentity, R2_BitScore
#    - Reorganize these columns so that all the Read1 data is together, followed by all the Read2 data


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

my $read_records_infile = '';
my $blast_results_dir = '';
my $reads_analyzed_per_block = '';
my $working_dir = '';
my $total_reads_analyzed = '';
my $summary_outfile = '';
my $min_bit_score = '';    #Hits below the minimum Bit Score will not be considered reliable ...
my $upper_bit_score = '';  #The Bit Score used to assign a read to this source database (should reflect a perfect or near perfect alignment)
my $min_align_length = '';
my $logfile = '';

GetOptions ('read_records_infile=s'=>\$read_records_infile, 'blast_results_dir=s'=>\$blast_results_dir,
	    'reads_analyzed_per_block=i'=>\$reads_analyzed_per_block, 'total_reads_analyzed=i'=>\$total_reads_analyzed,
	    'working_dir=s'=>\$working_dir, 'summary_outfile=s'=>\$summary_outfile,
	    'min_bit_score=f'=>\$min_bit_score, 'upper_bit_score=f'=>\$upper_bit_score, 'min_align_length=i'=>\$min_align_length,
	    'logfile=s'=>\$logfile);

#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tThis parses blast results for Solexa paired reads BLASTED against a database of Human Repeat Elements", RESET;
print GREEN, "\n\tThis script takes a read record file and blast results files as input, parses the results and creates a summary file", RESET;
print GREEN, "\n\tIt also updates the master read record file to assign perfect matching reads to this source (if they have not already been assigned)", RESET;
print GREEN, "\n\tSpecify a tab-delimited master read records file using: --read_records_infile", RESET;
print GREEN, "\n\tSpecify a directory containing tab-delimited blast results using: --blast_results_dir", RESET;
print GREEN, "\n\tSpecify the number of reads that were blasted to create each blast results file using: --reads_analyzed_per_block", RESET;
print GREEN, "\n\t\tThis number will be used to get an initial sense of the mapping success rate", RESET;
print GREEN, "\n\t\tCheck the corresponding fasta file directory to confirm the number of reads in each block", RESET;
print GREEN, "\n\tSimilarly specify the grand total number of reads blasted using: --total_reads_analyzed", RESET;
print GREEN, "\n\tSpecify a working directory for temp files using: --working_dir", RESET;
print GREEN, "\n\tSpecify the name of the resulting summary of blast hits for each read pair using: --summary_outfile", RESET;
print GREEN, "\n\tSpecify the minimum bit score for a blast hit to be considered reliable for resolving ambiguous hits with pairing info using: --min_bit_score", RESET;
print GREEN, "\n\t\tNOTE: that this min_bit_score is only used in this very limited sense - actual quality filtering will be done later!", RESET;
print GREEN, "\n\tSpecify the minimum bit score required for a read to be assigned to this source database using: --upper_bit_score", RESET;
print GREEN, "\n\t\tThis cutoff should represent a perfect or near perfect BLAST hit - if the cutoff is met the read status will be changed in the master read record file", RESET;
print GREEN, "\n\t\tYou can try different cutoffs with the same BLAST results data and the read record file will be updated", RESET;
print GREEN, "\n\tSpecify the minimum alignment length required for a read to be a assigned to this source database using: --min_align_length", RESET;
print GREEN, "\n\tSpecify a log file to use for script output using: --logfile", RESET;

print GREEN, "\n\nExample: parseRepeatBlastResults.pl  --read_records_infile=/projects/malachig/solexa/read_records/HS04391/20821AAXX_Lane1.txt.gz  --blast_results_dir=/projects/malachig/solexa/blast_results/HS04391/20821AAXX_Lane1/human_and_ancestral_repeats_RepBase13.05/  --reads_analyzed_per_block=125000  --total_reads_analyzed=5738676  --working_dir=/projects/malachig/solexa/blast_results/HS04391/20821AAXX_Lane1/temp/  --summary_outfile=/projects/malachig/solexa/read_records/HS04391/Repeats/20821AAXX_Lane1_Repeats.txt  --min_bit_score=48.1  --upper_bit_score=60.0  --min_align_length=32  --logfile=/projects/malachig/solexa/logs/HS04391/20821AAXX_Lane1/parseRepeatBlastResults_LOG.txt\n\n", RESET;

unless ($read_records_infile && $blast_results_dir && $reads_analyzed_per_block && $total_reads_analyzed && $working_dir && $summary_outfile && $min_bit_score && $upper_bit_score && $min_align_length && $logfile){
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  exit();
}

#Check dirs before getting started
$working_dir = &checkDir('-dir'=>$working_dir, '-clear'=>"no");
$blast_results_dir = &checkDir('-dir'=>$blast_results_dir, '-clear'=>"no");

my $base;
if ($read_records_infile =~ /(.*)\.gz$/){
  $base = $1;  
}else{
  print RED, "\nFormat of read records infile name not understood: $read_records_infile - not compressed?\n\n", RESET;
  exit();
}


open (LOG, ">$logfile") || die "\nCould not open logfile: $logfile\n\n";

#1.)  Parse blast results files one at a time
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

my $repeat_u_count = 0;
my $repeat_a_count = 0;

my $library_name;
if ($read_records_infile =~ /\/(\w+_Lane\d+)\.txt/){
  $library_name = $1;
}else{
  print RED, "\nCould not determine library name from read records infile: $read_records_infile\n\n", RESET;
  exit();
}

#Open global Berkley DB files for reading or writing
my $cache_size = 25600000;
my %top_hits_db;
my $top_hits_db_file = "$working_dir"."$library_name"."_TopHits.btree";
system ("rm -f $top_hits_db_file");

tie(%top_hits_db, 'BerkeleyDB::Btree', -Cachesize => $cache_size, -Filename=> $top_hits_db_file, -Flags => DB_CREATE) or die "can't open file $top_hits_db_file: $! $BerkeleyDB::Error\n";

&parseBlastFiles('-blast_results_dir'=>$blast_results_dir);
my $message = &memoryUsage();
print YELLOW, "\n\n$message", RESET;
print LOG "\n\n$message";

#Cleanly close the berkeley DB files
untie(%top_hits_db);

#Retie top hits file as readonly
tie(%top_hits_db, 'BerkeleyDB::Btree', -Cachesize => $cache_size, -Filename=> $top_hits_db_file, -Flags => DB_RDONLY ) or die "can't open file $top_hits_db_file: $! $BerkeleyDB::Error\n";


#3.) Create a summary of top hits and associated info with both reads of a pair on a single line
#    - Count entries in the read records file.  Store all read IDs as an array.
#    - Note that not all read records will have a top hits record.  There may not have been any blast hits to either read of the pair
#    - Parse a BLOCK of full read records.  Create a hash and key on Read_ID
#    - Now parse the top hits file and build a second hash for all individual reads that correspond to those in the BLOCK
#    - Once these two hashes are built go through the read records hash and print out corresponding top hits to a new summary file
#      - If a particular read of a read pair in the read records hash does not have a top hits record, fill in values with 'NA'
#      - Each read record will have HitType value which is either 'Top_Hit', 'Ambiguous', or 'None'
#    - Initialize both hashes and start processing the next block

#    - Columns that will be in the summary file:
#    - Read_ID, DistanceBetweenReads, R1_ID, R1_HitType, R1_RepeatName[repeat|family|source], R1_AlignmentLength, R1_PercentIdentity, R1_BitScore, R2_ID, R2_HitType, R2_RepeatName[repeat|family|source], R2_AlignmentLength, R2_PercentIdentity, R2_BitScore

&joinRecords('-read_records_infile'=>$read_records_infile, '-summary_outfile'=>$summary_outfile, '-library_name'=>$library_name, '-working_dir'=>$working_dir);
$message = &memoryUsage();
print YELLOW, "\n\n$message", RESET;
print LOG "\n\n$message";

print "\n\n";
print LOG "\n\n";

untie(%top_hits_db);
print BLUE, "\nDeleting top hits db file: $top_hits_db_file\n\n", RESET;
print LOG "\nDeleting top hits db file: $top_hits_db_file\n\n";
my $cmd = "rm -f $top_hits_db_file";
system($cmd);

$message = &memoryUsage();
print YELLOW, "\n\n$message", RESET;
print LOG "\n\n$message";

print LOG "\n\nSCRIPT COMPLETE\n\n";
close (LOG);
print "\n\nSCRIPT COMPLETE\n\n";

exit();



#############################################################################################################################
#1.) Parse blast results files one at a time                                                                                #
#2.)  Print the best hit to the 'top_hits' file                                                                             #
#############################################################################################################################
sub parseBlastFiles{
  my %args = @_;
  my $blast_results_dir = $args{'-blast_results_dir'};

  #Get files from this directory
  print BLUE, "\n\nSearching $blast_results_dir for blast results files", RESET;
  print LOG "\n\nSearching $blast_results_dir for blast results files";

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
	print LOG "\n\t$file_path  is a directory - skipping";
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
	  print LOG "\n\t$file_path does not appear to be a blast tabular output file!!\n\n";
	  close(TEST);
	  next(FILE_TEST);
	}
	$line_count++;
      }
      close(TEST);
      $file_count++;
      print BLUE, "\n\t$file_path was added to the list of files to be processed", RESET;
      print LOG "\n\t$file_path was added to the list of files to be processed";

      $blast_files{$file_count}{path} = $file_path;
      $blast_files{$file_count}{first_read_pair_id} = $first_read_pair_id;
    }

  #Parse blast files one at a time and write top hits to berkley DB file as each read is processed
  my %temp_blast_results;
  my %stored_blast_results;

  my $num_files = keys %blast_files;
  print BLUE, "\n\nBegin parsing $num_files blast results files", RESET;
  print LOG "\n\nBegin parsing $num_files blast results files";

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
    print LOG "\n\n\tParsing $results_file for blast results";

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

      #If this read_id is different from the working read pair ID, a new block of hits has started:
      #reset the working read id, store the top N hits for this read, reset the blast_results hash, and continue
      if ($read_pair_id ne $working_read_pair_id){

	my $reads_found = keys %temp_blast_results;
	$read_ids_processed += $reads_found;
	$gt_read_ids_processed += $reads_found;

	&processReadPairBlastHits('-temp_blast_results'=>\%temp_blast_results, '-read_pair_id'=>$working_read_pair_id);

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
      &processReadPairBlastHits('-temp_blast_results'=>\%temp_blast_results, '-read_pair_id'=>$working_read_pair_id);
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
    print BLUE, "\n\tReads which could be assigned to a Repeat as a PAIR:", RESET;
    print BLUE, "\n\t\tFound a total of $unambig_readpair_reads ($percent_reads_unambig_as_pair%) reads belonging to unambiguous read-pair mappings", RESET;
    print BLUE, "\n\t\tFound a total of $ambig_readpair_reads ($percent_reads_ambig_as_pair%) reads belonging to ambiguous read-pair mappings", RESET;
    print BLUE, "\n\tReads which had to be assigned to a Repeat INDIVIDUALLY:", RESET;
    print BLUE, "\n\t\tFound a total of $unambig_individual_reads ($percent_unambig_individual_reads%) individual reads mapping unambiguously", RESET;
    print BLUE, "\n\t\tFound a total of $ambig_individual_reads mapping ($percent_ambig_individual_reads%) ambiguously", RESET;

    print LOG "\n\tProcessed $blast_results_processed blast results";
    print LOG "\n\tFound $redundant_hits of equal or lesser quality to single subject sequences";
    print LOG "\n\tStored a total of $blast_results_stored blast results corresponding to $read_ids_processed reads (only the top hit will be printed)\n";

    print LOG "\n\tFound hits for $read_ids_processed ($percent_reads_with_hits%) of $reads_analyzed_per_block reads in this block of reads (block size specified by user)";
    print LOG "\n\tReads which could be assigned to a Repeat as a PAIR:";
    print LOG "\n\t\tFound a total of $unambig_readpair_reads ($percent_reads_unambig_as_pair%) reads belonging to unambiguous read-pair mappings";
    print LOG "\n\t\tFound a total of $ambig_readpair_reads ($percent_reads_ambig_as_pair%) reads belonging to ambiguous read-pair mappings";
    print LOG "\n\tReads which had to be assigned to a Repeat INDIVIDUALLY:";
    print LOG "\n\t\tFound a total of $unambig_individual_reads ($percent_unambig_individual_reads%) individual reads mapping unambiguously";
    print LOG "\n\t\tFound a total of $ambig_individual_reads mapping ($percent_ambig_individual_reads%) ambiguously";

  }
  my $total_reads_with_hits = keys %top_hits_db;
  print BLUE, "\nFound at least one hit for a grand total of $total_reads_with_hits reads\n", RESET;
  print LOG "\nFound at least one hit for a grand total of $total_reads_with_hits reads\n";

  #Summarize what was found for ALL BLAST FILES - Summarize relative to all possible reads that were blasted ($total_reads_analyzed specified by the user)
  my $gt_percent_reads_with_hits = sprintf("%.2f", (($gt_read_ids_processed / $total_reads_analyzed)*100));
  my $gt_percent_reads_unambig_as_pair = sprintf("%.2f", (($gt_unambig_readpair_reads / $total_reads_analyzed)*100));
  my $gt_percent_reads_ambig_as_pair = sprintf("%.2f", (($gt_ambig_readpair_reads / $total_reads_analyzed)*100));

  my $gt_percent_unambig_individual_reads = sprintf("%.2f", (($gt_unambig_individual_reads / $total_reads_analyzed)*100));
  my $gt_percent_ambig_individual_reads = sprintf("%.2f", (($gt_ambig_individual_reads / $total_reads_analyzed)*100));

  print BLUE, "\n\nGRAND TOTAL STATS", RESET;
  print BLUE, "\n\tProcessed $gt_blast_results_processed blast results", RESET;
  print BLUE, "\n\tFound $gt_redundant_hits of equal or lesser quality to single Repeat sequences", RESET;
  print BLUE, "\n\tStored a total of $gt_blast_results_stored blast results corresponding to $gt_read_ids_processed reads (only the top hit will be printed)\n", RESET;

  print BLUE, "\n\tFound hits for $gt_read_ids_processed ($gt_percent_reads_with_hits%) of $total_reads_analyzed total reads analyzed (specified by user)", RESET;
  print BLUE, "\n\tReads which could be assigned to a Repeat as a PAIR:", RESET;
  print BLUE, "\n\t\tFound a total of $gt_unambig_readpair_reads ($gt_percent_reads_unambig_as_pair%) reads belonging to unambiguous read-pair mappings", RESET;
  print BLUE, "\n\t\tFound a total of $gt_ambig_readpair_reads ($gt_percent_reads_ambig_as_pair%) reads belonging to ambiguous read-pair mappings", RESET;
  print BLUE, "\n\tReads which had to be assigned to a Repeat INDIVIDUALLY:", RESET;
  print BLUE, "\n\t\tFound a total of $gt_unambig_individual_reads ($gt_percent_unambig_individual_reads%) individual reads mapping unambiguously", RESET;
  print BLUE, "\n\t\tFound a total of $gt_ambig_individual_reads ($gt_percent_ambig_individual_reads%) reads mapping ambiguously", RESET;
  print BLUE, "\n\tNOTE: Because of cases were either R1 or R2 has NO hits, these numbers will not correspond exactly to the number of lines in your file!", RESET;

  print LOG "\n\nGRAND TOTAL STATS";
  print LOG "\n\tProcessed $gt_blast_results_processed blast results";
  print LOG "\n\tFound $gt_redundant_hits of equal or lesser quality to single Repeat sequences";
  print LOG "\n\tStored a total of $gt_blast_results_stored blast results corresponding to $gt_read_ids_processed reads (only the top hit will be printed)\n";

  print LOG "\n\tFound hits for $gt_read_ids_processed ($gt_percent_reads_with_hits%) of $total_reads_analyzed total reads analyzed (specified by user)";
  print LOG "\n\tReads which could be assigned to a Repeat as a PAIR:";
  print LOG "\n\t\tFound a total of $gt_unambig_readpair_reads ($gt_percent_reads_unambig_as_pair%) reads belonging to unambiguous read-pair mappings";
  print LOG "\n\t\tFound a total of $gt_ambig_readpair_reads ($gt_percent_reads_ambig_as_pair%) reads belonging to ambiguous read-pair mappings";
  print LOG "\n\tReads which had to be assigned to a Repeat INDIVIDUALLY:";
  print LOG "\n\t\tFound a total of $gt_unambig_individual_reads ($gt_percent_unambig_individual_reads%) individual reads mapping unambiguously";
  print LOG "\n\t\tFound a total of $gt_ambig_individual_reads ($gt_percent_ambig_individual_reads%) reads mapping ambiguously";
  print LOG "\n\tNOTE: Because of cases were either R1 or R2 has NO hits, these numbers will not correspond exactly to the number of lines in your file!";

  return();
}


#############################################################################################################################
#Process the blast hits found for a pair of reads
#############################################################################################################################
sub processReadPairBlastHits{
  my %args = @_;
  my %read_pair_blast_results = %{$args{'-temp_blast_results'}};
  my $read_pair_id = $args{'-read_pair_id'};

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

    #DEBUG
    unless ($read1_hits_ref->{$subject_id} && $read2_hits_ref->{$subject_id}){
      print RED, "\n$read_pair_id\tSingle Paired Read - Unambiguous.  Subject ID: $subject_id\tREAD1: $read1_id\tREAD2: $read2_id", RESET;
      #print Dumper $read1_hits_ref;
      #print Dumper $read2_hits_ref;
      print Dumper %read_pair_blast_results;
      print Dumper %subject_hit_list;
      print Dumper %ranked_subject_pair_hits;
      exit();
    }

    #    - Read_ID, Hit_Type (i.e. Ambiguous or Top Hit), Subject_ID, AlignmentLength, PercentIdentity, BitScore, SubjectStart, SubjectEnd, strand

    $top_hits_db{$read1_id} = "Top_Hit\t$subject_id\t$read1_hits_ref->{$subject_id}->{al}\t$read1_hits_ref->{$subject_id}->{pi}\t$read1_hits_ref->{$subject_id}->{bs}\t$read1_hits_ref->{$subject_id}->{ss}\t$read1_hits_ref->{$subject_id}->{se}\t$read1_hits_ref->{$subject_id}->{hs}";

    $top_hits_db{$read2_id} = "Top_Hit\t$subject_id\t$read2_hits_ref->{$subject_id}->{al}\t$read2_hits_ref->{$subject_id}->{pi}\t$read2_hits_ref->{$subject_id}->{bs}\t$read2_hits_ref->{$subject_id}->{ss}\t$read2_hits_ref->{$subject_id}->{se}\t$read2_hits_ref->{$subject_id}->{hs}";

    $unambig_readpair_reads += 2;
    $gt_unambig_readpair_reads += 2;

    #print YELLOW, "\n\tPaired Read Mapping. Single subject - two quality hits", RESET;
    return();

  }elsif($read_pair_subject_count > 1){

    #Multiple subjects found that have quality hits for both reads.  See if there is a clear winner according to combined score

    my $read1_hits_ref = $read_pair_blast_results{$read1_id}{blast_hits};
    my $read2_hits_ref = $read_pair_blast_results{$read2_id}{blast_hits};

    if ($ranked_subject_pair_hits{1}{combined_score} > $ranked_subject_pair_hits{2}{combined_score}){

      my $subject_id = $ranked_subject_pair_hits{1}{subject_id};

      #DEBUG
      unless ($read1_hits_ref->{$subject_id} && $read2_hits_ref->{$subject_id}){
	print RED, "\nMultiple Paired Reads - Unambiguous.  Subject ID: $subject_id\tREAD1: $read1_id\tREAD2: $read2_id", RESET;
	print Dumper $read1_hits_ref;
	print Dumper $read2_hits_ref;
	exit();
      }

      $top_hits_db{$read1_id} = "Top_Hit\t$subject_id\t$read1_hits_ref->{$subject_id}->{al}\t$read1_hits_ref->{$subject_id}->{pi}\t$read1_hits_ref->{$subject_id}->{bs}\t$read1_hits_ref->{$subject_id}->{ss}\t$read1_hits_ref->{$subject_id}->{se}\t$read1_hits_ref->{$subject_id}->{hs}";
      $top_hits_db{$read2_id} = "Top_Hit\t$subject_id\t$read2_hits_ref->{$subject_id}->{al}\t$read2_hits_ref->{$subject_id}->{pi}\t$read2_hits_ref->{$subject_id}->{bs}\t$read2_hits_ref->{$subject_id}->{ss}\t$read2_hits_ref->{$subject_id}->{se}\t$read2_hits_ref->{$subject_id}->{hs}";

      $unambig_readpair_reads += 2;
      $gt_unambig_readpair_reads += 2;

      #print YELLOW, "\n\tPaired Read Mapping. Multiple subjects - two quality hits - but a 'best' one was found", RESET;
      return();
    }else{

      my $subject_id = $ranked_subject_pair_hits{1}{subject_id};

      #DEBUG
      unless ($read1_hits_ref->{$subject_id} && $read2_hits_ref->{$subject_id}){
	print RED, "\nPaired Read - Ambiguous.  Subject ID: $subject_id\tREAD1: $read1_id\tREAD2: $read2_id", RESET;
	print Dumper $read1_hits_ref;
	print Dumper $read2_hits_ref;
	exit();
      }

      $top_hits_db{$read1_id} = "Ambiguous\tNA\t$read1_hits_ref->{$subject_id}->{al}\t$read1_hits_ref->{$subject_id}->{pi}\t$read1_hits_ref->{$subject_id}->{bs}\tNA\tNA\tNA";
      $top_hits_db{$read2_id} = "Ambiguous\tNA\t$read2_hits_ref->{$subject_id}->{al}\t$read2_hits_ref->{$subject_id}->{pi}\t$read2_hits_ref->{$subject_id}->{bs}\tNA\tNA\tNA";

      $ambig_readpair_reads += 2;
      $gt_ambig_readpair_reads += 2;

      #print YELLOW, "\n\tPaired Read Mapping. Multiple subjects - two quality hits - ambiguous", RESET;
      return();
    }
  }

  #C.) There are no subjects where both reads map to a single subject.  OR the BIT SCORE of one or both of these reads is below the cutoff.
  #   - Treat these as individual reads.  Use the bit score to identify the 'top hit'.  If there are multiple equal top hits, mark the read as Ambiguous
  #   - Once both reads have been processed, return().

  #print YELLOW, "\n\tNo Single Subject for Read Pair - Treating reads individually", RESET;

  foreach my $working_read_id (sort {$read_pair_blast_results{$a}{read_num} <=> $read_pair_blast_results{$b}{read_num}} keys %read_pair_blast_results){

    #DEBUG
    #unless ($read_pair_blast_results{$working_read_id}{read_num} =~ /\d+/){
    #  print Dumper %read_pair_blast_results;
    #  exit();
    #}

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

	$top_hits_db{$working_read_id} = "Top_Hit\t$subject_id\t$hits_ref->{$subject_id}->{al}\t$hits_ref->{$subject_id}->{pi}\t$hits_ref->{$subject_id}->{bs}\t$hits_ref->{$subject_id}->{ss}\t$hits_ref->{$subject_id}->{se}\t$hits_ref->{$subject_id}->{hs}";
	#print YELLOW, "\n\tSingle Read Mapping. Either there was only one hit or a 'best' one was found", RESET;

	$unambig_individual_reads++;
	$gt_unambig_individual_reads++;
	last (BS);

      }else{

	#DEBUG
	unless ($hits_ref->{$subject_id}){
	  print RED, "\nSingle Read - Ambiguous.  Subject ID: $subject_id\tREAD: $working_read_id", RESET;
	  print Dumper $hits_ref;
	  exit();
	}

	$top_hits_db{$working_read_id} = "Ambiguous\tNA\t$hits_ref->{$subject_id}->{al}\t$hits_ref->{$subject_id}->{pi}\t$hits_ref->{$subject_id}->{bs}\tNA\tNA\tNA";
	#print YELLOW, "\n\tSingle Read Mapping. Multiple subjects - two quality hits - ambiguous", RESET;

	$ambig_individual_reads++;
	$gt_ambig_individual_reads++;
	last (BS);

      }
    }
  }
  return();
}


#############################################################################################################################
#3.) Create a summary of top hits and associated info with both reads of a pair on a single line
#############################################################################################################################
sub joinRecords{
  my %args = @_;
  my $read_records_file = $args{'-read_records_infile'};
  my $summary_outfile = $args{'-summary_outfile'};
  my $library_name = $args{'-library_name'};
  my $working_dir = $args{'-working_dir'};

  my %read_pairs_db;
  my $read_pairs_db_file = "$working_dir"."$library_name"."_ReadPairs.btree";
  system ("rm -f $read_pairs_db_file");
  tie(%read_pairs_db, 'BerkeleyDB::Btree', -Cachesize => $cache_size, -Filename=> $read_pairs_db_file, -Flags => DB_CREATE) or die "can't open file $read_pairs_db_file: $! $BerkeleyDB::Error\n";

  print BLUE, "\n\nCreating Berkley DB of read records file.\n", RESET;
  print LOG "\n\nCreating Berkley DB of read records file.\n";


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
  $message = &memoryUsage();
  print YELLOW, "\n\n$message", RESET;
  print LOG "\n\n$message";

  #Now go through every read_pair_id, join with read1 and read2 BLAST hit info, and print result to outfile
  print BLUE, "\n\nCreate summary file: $summary_outfile\n", RESET;
  print LOG "\n\nCreate summary file: $summary_outfile\n";

  open (SUMMARY, ">$summary_outfile") || die "\nCould not open summary outfile: $summary_outfile\n\n";

  print SUMMARY "Read_ID\tDistanceBetweenReads\tR1_ID\tR1_HitType\tR1_RepeatName\tR1_Strand\tR1_AlignmentLength\tR1_PercentIdentity\tR1_BitScore\tR2_ID\tR2_HitType\tR2_RepeatName\tR2_Strand\tR2_AlignmentLength\tR2_PercentIdentity\tR2_BitScore\n";

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
    unless ($top_hits_db{$read1_id} || $top_hits_db{$read2_id}){
      next();
    }
    $read_pairs_to_be_printed++;

    #Format of blast hit record stored:
    #Top_Hit\t$subject_id\t$hits_ref->{$subject_id}->{al}\t$hits_ref->{$subject_id}->{pi}\t$hits_ref->{$subject_id}->{bs}\t$hits_ref->{$subject_id}->{ss}\t$hits_ref->{$subject_id}->{se}\t$hits_ref->{$subject_id}->{hs}


    #Determine the fragment size based on the alignment of the two reads.  Only consider cases where both reads have 'Top_Hits'
    #Use the outer distance of the coords
    my $distance = "NA";

    if ($top_hits_db{$read1_id} && $top_hits_db{$read2_id}){
      my $string1 = $top_hits_db{$read1_id};
      my @string1 = split("\t", $string1);
      my $string2 = $top_hits_db{$read2_id};
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
	  $distance = "Repeat_Mismatch";
	}
      }
    }

    print SUMMARY "$read_pair_id\t$distance";

    #READ 1#
    if ($top_hits_db{$read1_id}){
      my $string = $top_hits_db{$read1_id};
      my @string = split("\t", $string);

      my $hit_type = $string[0];
      my $subject_id = $string[1];
      my $align_length = $string[2];
      my $percent_id = $string[3];
      my $bit_score = $string[4];
      my $strand = $string[7];

      #print SUMMARY "\tR1_ID\tR1_HitType\tR1_RepeatName\tR1_Strand\tR1_AlignmentLength\tR1_PercentIdentity\tR1_BitScore";
      print SUMMARY "\t$read1_id\t$hit_type\t$subject_id\t$strand\t$align_length\t$percent_id\t$bit_score";

      #If the upper bit score criteria is met and the status of this read is still 'Unassigned' in the read record file, change the status
      #Possible status is 'REPEAT_U' or 'REPEAT_A'

      #If a read has already been assigned to 'Repeat_U' or 'Repeat_A', start by resetting these to 'Unassigned'
      if ($read1_status eq "Repeat_U" || $read1_status eq "Repeat_A"){
	$read1_status = "Unassigned";
      }

      #Make sure that the status of this reads is still 'Unassigned'
      if ($read1_status eq "Unassigned"){

	#Now, if the top hit for this read meets the user specified threshold, change its status
	#- The resulting status will depend on whether the hit was ambiguous or not ...
	if (($hit_type eq "Top_Hit") && ($bit_score >= $upper_bit_score) && ($align_length >= $min_align_length)){
	  $repeat_u_count++;
	  $read1_status = "Repeat_U";
	}elsif(($hit_type eq "Ambiguous") && ($bit_score >= $upper_bit_score) && ($align_length >= $min_align_length)){
	  $repeat_a_count++;
	  $read1_status = "Repeat_A";
	}
      }
    }else{
      #Read 1 is not defined but READ 2 was
      print SUMMARY "\t$read1_id\tNone\tNA\tNA\tNA\tNA\tNA";
    }

    #READ 2#
    if ($top_hits_db{$read2_id}){
      my $string = $top_hits_db{$read2_id};
      my @string = split("\t", $string);

      my $hit_type = $string[0];
      my $subject_id = $string[1];
      my $align_length = $string[2];
      my $percent_id = $string[3];
      my $bit_score = $string[4];
      my $strand = $string[7];

      #print SUMMARY "\tR2_ID\tR2_HitType\tR2_RepeatName\tR2_AlignmentLength\tR2_PercentIdentity\tR2_BitScore\n";
      print SUMMARY "\t$read2_id\t$hit_type\t$subject_id\t$strand\t$align_length\t$percent_id\t$bit_score\n";

      #If the upper bit score criteria is met and the status of this read is still 'Unassigned' in the read record file, change the status
      #Possible status is 'Repeat_U' or 'Repeat_A'

      #If a read has already been assigned to 'Repeat_U' or 'Repeat_A', start by resetting these to 'Unassigned'
      if ($read2_status eq "Repeat_U" || $read2_status eq "Repeat_A"){
	$read2_status = "Unassigned";
      }

      #Make sure that the status of this reads is still 'Unassigned'
      if ($read2_status eq "Unassigned"){

	#Now, if the top hit for this read meets the user specified threshold, change its status
	#- The resulting status will depend on whether the hit was ambiguous or not ...
	if (($hit_type eq "Top_Hit") && ($bit_score >= $upper_bit_score) && ($align_length >= $min_align_length)){
	  $repeat_u_count++;
	  $read2_status = "Repeat_U";
	}elsif(($hit_type eq "Ambiguous") && ($bit_score >= $upper_bit_score) && ($align_length >= $min_align_length)){
	  $repeat_a_count++;
	  $read2_status = "Repeat_A";
	}
      }
    }else{
      #Read 2 is not defined but Read 1 was
      print SUMMARY "\t$read2_id\tNone\tNA\tNA\tNA\tNA\tNA\n";

    }
    $read_pairs_db{$read_pair_id} = "$read1_status\t$read2_status";
  }
  close(SUMMARY);
  $message = &memoryUsage();
  print YELLOW, "\n\n$message", RESET;
  print LOG "\n\n$message";
  print BLUE, "\n\nFound hits for one or both reads of $read_pairs_to_be_printed read pairs\n", RESET;

  print BLUE, "\n\nCompress summary outfile: $summary_outfile\n", RESET;
  my $compress_cmd = "gzip -f $summary_outfile";
  system($compress_cmd);


  #Finally go through the original read records file and update the read status values where appropriate
  print BLUE, "\n\nUpdate read status values in read record file: $read_records_file\n", RESET;
  print LOG "\n\nUpdate read status values in read record file: $read_records_file\n";

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
  print LOG "\n\n$message";

  untie(%read_pairs_db);

  my $rm_cmd2 = "rm -f $read_pairs_db_file";
  system($rm_cmd2);

  my $percent_u = sprintf("%.2f", (($repeat_u_count/$read_count)*100));
  my $percent_a = sprintf("%.2f", (($repeat_a_count/$read_count)*100));

  print BLUE, "\n\nAssigned $repeat_u_count reads ($percent_u%) to the status 'Repeat_U', and $repeat_a_count ($percent_a%) reads to the status 'Repeat_A'\n\n", RESET;
  print LOG "\n\nAssigned $repeat_u_count reads ($percent_u%) to the status 'Repeat_U', and $repeat_a_count ($percent_a%) reads to the status 'Repeat_A'\n\n";

  my $cmd_gzip = "gzip -f $new_read_record_file";
  print BLUE, "\n\nCompressing the new read records file:\n\t$cmd_gzip\n\n", RESET;
  print LOG "\n\nCompressing the new read records file:\n\t$cmd_gzip\n\n";
  system ($cmd_gzip);

  my $cmd_mv = "mv $new_read_record_file".".gz"." $read_records_file";
  print BLUE, "\n\nOverwriting existing read records file:\n\t$cmd_mv\n\n", RESET;
  print LOG "\n\nOverwriting existing read records file:\n\t$cmd_mv\n\n";
  system($cmd_mv);

  return();
}


