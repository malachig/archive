#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#The purpose of this script is to parse blast results files for Solexa Paired Reads that were BLASTed against all EnsEMBL transcripts

#1.) First get a hash from ALEXA containing all Transcript-to-Gene mappings

#2.)  Parse blast results files one at a time
#2a.) Go through the specified blast results files (in the specified directory)
#    - For each blast results file go through the hits for a single read
#    - Store the best hit at the GENE level (i.e. keep the best hit for the transcripts of a gene)
#    - If the best and second best hit have the same bit score, consider this a tie and make note of it as an ambiguous gene hit

#2b.) Print the best hit to the 'top_hits' file
#    - Read_ID, Hit_Type (i.e. Ambiguous or Top Hit), Gene_ID, Trans_ID (best hit),  AlignLength, Percent Identity, Subject Start, Subject End

#3.) At this point, also calculate the transcript position of the hit.
#    - Parse the top_hits file and append this info to each line
#    - Take the centre of the query hit coordinates
#    - Determine the relative position of this hit in the TRANSCRIPT that it mapped to (i.e. percent distance from 3' end)
#    - Only do this for 'Top Hits', this value will be NA for Ambiguous hits
#    - For convenience, also write out the transcript size for each hit
#    - The new columns will be called 'TranscriptSize' and 'RelativePosition'

#4.) At this point, also calculate the CHROMOSOME coordinates (gapped) of each top hit
#    - Once again this can be done while parsing the top_hits file.
#    - Note that it is possible for one or both reads to span an exon junction.
#    - So the transcript coordinates from mapping to EnsEMBL transcripts have to be mapped back to the genome by taking the introns into account
#    - This is similar to display protein feature from EnsEMBL on the UCSC genome browser (these are also reported relative to the transcript)
#    - Append two columns to the top_hits file: start_coords (array of chromosome start positions) and end_coord (array of chromosome end positions)

#5.) Join top hits results for each read of a pair to a single line of data and print as a summary of the blast results
#    - Count entries in the read records file.  Store all read IDs as an array.
#    - Note that not all read records will have a top hits record.  There may not have been any blast hits to either read of the pair
#    - Parse a BLOCK of full read records.  Create a hash and key on Read_ID
#    - Now parse the top hits file and build a second hash for all individual reads that correspond to those in the BLOCK
#    - Once these two hashes are built go through the read records hash and print out corresponding top hits to a new summary file
#      - If a particular read of a read pair in the read records hash does not have a top hits record, fill in values with 'NA'
#      - Each read record will have topGeneHit value which is either the Gene_ID, 'Ambiguous', or 'None'
#    - Initialize both hashes and start processing the next block

#    - Columns to be added to the EnsEMBL transcript blast results summary file:
#    - Read_ID, DistanceBetweenReads_Genomic, DistanceBetweenReads_Transcript, R1_ID, R1_HitType, R1_GeneID, R1_GeneName, R1_AlignmentLength, R1_PercentIdentity, R1_BitScore, R1_TranscriptSize, R1_RelativePosition, R1_Chromosome, R1_ChrStartCoords, R1_ChrEndCoords, R2_ID, R2_HitType, R2_GeneID, R2_GeneName, R2_AlignmentLength, R2_PercentIdentity, R2_BitScore, R2_TranscriptSize, R2_RelativePosition, R2_Chromosome, R2_ChrStartCoords, R2_ChrEndCoords

#    - Organize these columns so that all the Read1 data is together, followed by all the Read2 data

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
use utilities::ALEXA_DB qw(:all);

my $database = '';
my $server = '';
my $user = '';
my $password = '';
my $working_dir = '';
my $read_records_infile = '';
my $blast_results_dir = '';
my $reads_analyzed_per_block = '';
my $total_reads_analyzed = '';
my $summary_outfile = '';
my $filter_pseudogene_hits = '';
my $min_bit_score = '';    #Hits below the minimum Bit Score will not be considered reliable ...
my $upper_bit_score = '';  #The Bit Score used to assign a read to this source database (should reflect a perfect or near perfect alignment)
my $min_align_length = '';
my $logfile = '';

GetOptions ('database=s'=>\$database,'server=s'=>\$server, 'user=s'=>\$user, 'password=s'=>\$password,
	    'read_records_infile=s'=>\$read_records_infile, 'blast_results_dir=s'=>\$blast_results_dir,
	    'reads_analyzed_per_block=i'=>\$reads_analyzed_per_block, 'total_reads_analyzed=i'=>\$total_reads_analyzed,
	    'working_dir=s'=>\$working_dir, 'summary_outfile=s'=>\$summary_outfile,
	    'filter_pseudogene_hits=s'=>\$filter_pseudogene_hits, 'min_bit_score=f'=>\$min_bit_score, 'upper_bit_score=f'=>\$upper_bit_score, 'min_align_length=i'=>\$min_align_length,
	    'logfile=s'=>\$logfile);

#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tThis parses blast results for Solexa paired reads BLASTED against all ALEXA/EnsEMBL transcripts", RESET;
print GREEN, "\n\tThis script takes a read record file and blast results files as input, parses the results and creates a summary file", RESET;
print GREEN, "\n\tIt also updates the master read record file to assign perfect matching reads to this source (if they have not already been assigned)", RESET;
print GREEN, "\n\tSpecify the ALEXA database and server to query using: --database and --server", RESET;
print GREEN, "\n\tSpecify the ALEXA user and password for access using: --user and --password\n", RESET;
print GREEN, "\n\tSpecify a working directory for temp files using: --working_dir", RESET;
print GREEN, "\n\tSpecify a tab-delimited master read records file using: --read_records_infile", RESET;
print GREEN, "\n\tSpecify a directory containing tab-delimited blast results using: --blast_results_dir", RESET;
print GREEN, "\n\tSpecify the number of reads that were blasted to create each blast results file using: --reads_analyzed_per_block", RESET;
print GREEN, "\n\t\tThis number will be used to get an initial sense of the mapping success rate", RESET;
print GREEN, "\n\t\tCheck the corresponding fasta file directory to confirm the number of reads in each block", RESET;
print GREEN, "\n\tSimilarly specify the grand total number of reads blasted using: --total_reads_analyzed", RESET;
print GREEN, "\n\tSpecify the name of the resulting summary of blast hits for each read pair using: --summary_outfile", RESET;
print GREEN, "\n\tIf you would like to filter out hits to pseudogenes use: --filter_pseudogene_hits=yes", RESET;
print GREEN, "\n\tSpecify the minimum bit score for a blast hit to be considered reliable for resolving ambiguous hits with pairing info using: --min_bit_score", RESET;
print GREEN, "\n\t\tNOTE: that this min_bit_score is only used in this very limited sense - actual quality filtering will be done later!", RESET;
print GREEN, "\n\tSpecify the minimum bit score required for a read to be assigned to this source database using: --upper_bit_score", RESET;
print GREEN, "\n\t\tThis cutoff should represent a perfect or near perfect BLAST hit - if the cutoff is met the read status will be changed in the master read record file", RESET;
print GREEN, "\n\t\tYou can try different cutoffs with the same BLAST results data and the read record file will be updated", RESET;
print GREEN, "\n\tSpecify the minimum alignment length required for a read to be a assigned to this source database using: --min_align_length", RESET;
print GREEN, "\n\tSpecify a log file to use for script output using: --logfile", RESET;

print GREEN, "\n\nExample: parseEnsemblTranscriptBlastResults.pl  --database=ALEXA_hs_49_36k  --server=jango.bcgsc.ca  --user=viewer  --password=viewer  --working_dir=/projects/malachig/solexa/blast_results/HS04391/temp/   --read_records_infile=/projects/malachig/solexa/read_records/HS04391/20821AAXX_Lane1.txt  --blast_results_dir=/projects/malachig/solexa/blast_results/HS04391/20821AAXX_Lane1/ensembl_transcripts_v49/  --reads_analyzed_per_block=125000  --total_reads_analyzed=5738676  --summary_outfile=/projects/malachig/solexa/read_records/HS04391/ENST_v49/20821AAXX_Lane1_ENST_v49.txt  --filter_pseudogene_hits=no  --min_bit_score=40.0  --upper_bit_score=71.9  --min_align_length=32  --logfile=/projects/malachig/solexa/logs/HS04391/20821AAXX_Lane1/parseEnsemblTranscriptBlastResults_LOG.txt\n\n", RESET;


unless ($database && $server && $user && $password && $read_records_infile && $blast_results_dir && $reads_analyzed_per_block && $total_reads_analyzed && $working_dir && $summary_outfile && ($filter_pseudogene_hits =~ /^yes|^no/i) && $logfile && $min_bit_score && $upper_bit_score && $min_align_length){
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

#1.) First get a hash from ALEXA containing all Transcript-to-Gene mappings

#Establish connection with the Alternative Splicing Expression database
my $alexa_dbh = &connectDB('-database'=>$database, '-server'=>$server, '-user'=>$user, '-password'=>$password);

my $genes_ref;
my $gene_transcripts_ref;
my $gene_exon_content_ref;
my %transcript_gene_map;
&getBasicGeneInfo();

#Close database connection
$alexa_dbh->disconnect();


#2.)  Parse blast results files one at a time
my $unambig_readpair_reads;
my $ambig_readpair_reads;
my $unambig_individual_reads;
my $ambig_individual_reads;

#Grand total counters
my $gt_blast_results_processed = 0;
my $gt_filtered_pseudogene_hits = 0;
my $gt_redundant_hits = 0;
my $gt_blast_results_stored = 0;
my $gt_read_ids_processed = 0;
my $gt_unambig_readpair_reads = 0;
my $gt_ambig_readpair_reads = 0;
my $gt_unambig_individual_reads = 0;
my $gt_ambig_individual_reads = 0;

my $enst_u_count = 0;
my $enst_a_count = 0;

my $library_name;
if ($read_records_infile =~ /\/(\w+_Lane\d+)\.txt/){
  $library_name = $1;
}else{
  print RED, "\nCould not determine library name from read records infile: $read_records_infile\n\n", RESET;
  exit();
}

#Open global Berkley DB files for reading or writing
my %top_hits_db;
my $top_hits_db_file = "$working_dir"."$library_name"."_TopHits.btree";
system ("rm -f $top_hits_db_file");
tie(%top_hits_db, 'BerkeleyDB::Btree', -Cachesize => 128000000, -Filename=> $top_hits_db_file, -Flags => DB_CREATE) or die "can't open file $top_hits_db_file: $! $BerkeleyDB::Error\n";

&parseBlastFiles('-blast_results_dir'=>$blast_results_dir);
my $message = &memoryUsage();
print YELLOW, "\n\n$message", RESET;
print LOG "\n\n$message";


#3.) At this point, also calculate the transcript position of the hit.
#    - Parse the db and append this info to each line
#    - Take the centre of the query hit coordinates
#    - Determine the relative position of this hit in the TRANSCRIPT that it mapped to (i.e. percent distance from 3' end)
#    - Only do this for 'Top Hits', this value will be NA for Ambiguous hits
#    - For convenience, also write out the transcript size for each hit
#    - The new columns will be called 'TranscriptSize' and 'RelativePosition'
&determineTranscriptPositions();
$message = &memoryUsage();
print YELLOW, "\n\n$message", RESET;
print LOG "\n\n$message";


#4.) At this point, also calculate the CHROMOSOME coordinates (gapped) of each top hit
#    - Once again this can be done while parsing the top_hits file.
#    - Note that it is possible for one or both reads to span an exon junction.
#    - So the transcript coordinates from mapping to EnsEMBL transcripts have to be mapped back to the genome by taking the introns into account
#    - This is similar to display protein feature from EnsEMBL on the UCSC genome browser (these are also reported relative to the transcript)
#    - Append two columns to the top_hits file: start_coords (array of chromosome start positions) and end_coord (array of chromosome end positions)
&getHitChromosomeCoordinates();
$message = &memoryUsage();
print YELLOW, "\n\n$message", RESET;
print LOG "\n\n$message";


#5.) Create a summary of top hits and associated info with both reads of a pair on a single line
#    - Count entries in the read records file.  Store all read IDs as an array.
#    - Note that not all read records will have a top hits record.  There may not have been any blast hits to either read of the pair
#    - Parse a BLOCK of full read records.  Create a hash and key on Read_ID
#    - Now parse the top hits file and build a second hash for all individual reads that correspond to those in the BLOCK
#    - Once these two hashes are built go through the read records hash and print out corresponding top hits to a new summary file
#      - If a particular read of a read pair in the read records hash does not have a top hits record, fill in values with 'NA'
#      - Each read record will have topHitType value which is either the Gene_ID, 'Ambiguous', or 'None'
#    - Initialize both hashes and start processing the next block

#    - Columns that will be in the summary file:
#    - [Read_ID, R1_ID] then  [R1_HitType (TopHit, Ambiguous, None), R1_GeneID, R1_GeneName, R1_AlignmentLength, R1_PercentIdentity, R1_TranscriptSize, R1_RelativePosition, R1_Chromosome, R1_ChrStartCoords, R1_ChrEndCoords]
#    - Then repeat for R2
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
#1.) First get a hash from ALEXA containing all Transcript-to-Gene mappings
#############################################################################################################################
sub getBasicGeneInfo{

  my @gene_ids = @{&getAllGenes ('-dbh'=>$alexa_dbh, '-gene_type'=>'All')};

  #Get the gene info for all genes
  print BLUE, "\nGetting gene data", RESET;
  print LOG "\nGetting gene data";
  my $gene_storable = "$database"."_AllGenes_GeneInfo_NoSeq.storable";
  $genes_ref = &getGeneInfo ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids, '-sequence'=>"no", '-storable'=>$gene_storable);

  #Get the transcript info for all transcripts of these genes
  print BLUE, "\nGetting transcript data as well as exons for each transcript", RESET;
  print LOG "\nGetting transcript data as well as exons for each transcript";
  my $trans_storable = "$database"."_AllGenes_TranscriptInfo_NoSeq.storable";
  $gene_transcripts_ref = &getTranscripts('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids, '-sequence'=>"no", '-storable'=>$trans_storable);

  #Get exon content for all genes
  print BLUE, "\nGetting exon content data for each gene", RESET;
  print LOG "\nGetting exon content data for each gene";
  my $ec_storable = "$database"."_AllGenes_ExonContent.storable";
  $gene_exon_content_ref = &getExonContent ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids, '-storable'=>$ec_storable);

  foreach my $gene_id (sort keys %{$gene_transcripts_ref}){
    my $transcripts_ref = $gene_transcripts_ref->{$gene_id}->{transcripts};

    foreach my $trans_id (keys %{$transcripts_ref}){
      $transcript_gene_map{$trans_id}{gene_id} = $gene_id;

      #Calculate the size of each transcript by adding up the size of its exons
      my $size = 0;
      my $exons_ref = $transcripts_ref->{$trans_id}->{exons};

      foreach my $exon_id (keys %{$exons_ref}){
	my $e_size = ($exons_ref->{$exon_id}->{exon_end} - $exons_ref->{$exon_id}->{exon_start})+1;
	$size += $e_size;
      }
      $transcript_gene_map{$trans_id}{transcript_size} = $size;
    }
  }
  return();
}


#############################################################################################################################
#2.) Parse blast results files one at a time                                                                                #
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

    $file_count++;
    print BLUE, "\n\t$file_path was added to the list of files to be processed", RESET;
    print LOG "\n\t$file_path was added to the list of files to be processed";

    $blast_files{$file_count}{path} = $file_path;
  }

  #Parse blast files one at a time and write top hits to a berkley db file as each read is processed
  my %temp_blast_results;
  my %stored_blast_results;

  my $num_files = keys %blast_files;
  print BLUE, "\n\nBegin parsing $num_files blast results files", RESET;
  print LOG "\n\nBegin parsing $num_files blast results files";

  foreach my $file_count (sort {$blast_files{$a}{path} cmp $blast_files{$b}{path}} keys %blast_files){

    my $blast_results_processed = 0;
    my $filtered_pseudogene_hits = 0;
    my $redundant_hits = 0;
    my $blast_results_stored = 0;
    my $read_ids_processed = 0;

    $unambig_readpair_reads = 0;
    $ambig_readpair_reads = 0;
    $unambig_individual_reads = 0;
    $ambig_individual_reads = 0;

    my $results_file = $blast_files{$file_count}{path};
    my $working_read_pair_id;
    my $read_pair_id;

    print YELLOW, "\n\n\tParsing $results_file for blast results", RESET;
    print LOG "\n\n\tParsing $results_file for blast results";

    #2a.) Go through the hits for each single read 
    #    - Store the best hit at the GENE level (i.e. keep the best hit for the transcripts of a gene)
    #    - If the best and second best hit have the same alignment length and score, consider this a tie and make note of it as an ambiguous gene hit
    open (BLAST_RESULTS, "zcat $results_file |") || die "\nCould not open blast results file: $results_file\n\n";

    my $line_count = 0;
    while (<BLAST_RESULTS>){
      if ($_ =~ /^\s*$/){next();}
      chomp($_);
      my @line = split ("\t", $_);
      my $read_id = $line[0];

      if ($line_count == 0){
        if ($read_id =~ /([a-zA-Z0-9]+_\d+_\d+_\d+_\d+)_(R\d)/){
          $working_read_pair_id = $1;
        }else{
          print RED, "\nCould not extract first read pair ID from the individual read ID: $read_id\n\n", RESET;
          exit();
        }
        my $line_elements = scalar(@line);

        unless ($line_elements == 12){
          print RED, "\n\t$results_file does not appear to be a blast tabular output file!!\n\n", RESET;
          print LOG "\n\t$results_file does not appear to be a blast tabular output file!!\n\n";
          exit();
        }
      }
      
      $line_count++;
      my $trans_id = $line[1];
      my $percent_identity = $line[2];
      my $alignment_length = $line[3];
      my $subject_start = $line[8]; #Position on the mRNA/EST sequence at which the read alignment begins
      my $subject_end = $line[9];   #Position on the mRNA/EST sequence at which the read alignment ends
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
      my $subject_size = $transcript_gene_map{$trans_id}{transcript_size};
      if ($subject_start < 1 || $subject_end < 1 || $subject_start > $subject_size || $subject_end > $subject_size){
        next();
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


      #Get the gene ID for this transcript hit
      my $gene_id = $transcript_gene_map{$trans_id}{gene_id};

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

      #If the user specified to filter out pseudogene hits, do this now
      if ($filter_pseudogene_hits =~ /^yes/i){
	my $gene_type = $genes_ref->{$gene_id}->{gene_type};
	if ($gene_type =~ /pseudo/){
	  $filtered_pseudogene_hits++;
	  $gt_filtered_pseudogene_hits++;
	  next();
	}
      }

      #Note: using the following datastructure will not allow multiple hits by one read to the same gene to be recorded
      #Similarly if there are multiple hits to different part of a single transcript, only one will be stored!
      #Only the most significant hit (longest) will be stored when this this occurs
      #If there are multiple hits with the same length (the one with a higher bit score will be chosen)

      #See if there are any hits recorded for this read already (one set for read1 the other for read2)
      if ($temp_blast_results{$read_id}){
	my $hits_ref = $temp_blast_results{$read_id}{blast_hits};

	#See if there has already been a hit for this gene (ie. multiple hits of one read to different transcripts of the same gene - or the same transcript even)
	if ($hits_ref->{$gene_id}){

	  #Replace the old hit record for this gene if this hit has a better bit score than the one already recorded
	  if ($bit_score > $hits_ref->{$gene_id}->{bs}){
	    $hits_ref->{$gene_id}->{pi} = $percent_identity;
	    $hits_ref->{$gene_id}->{al} = $alignment_length;
	    $hits_ref->{$gene_id}->{ss} = $subject_start;
	    $hits_ref->{$gene_id}->{se} = $subject_end;
	    $hits_ref->{$gene_id}->{bs} = $bit_score;
	    $hits_ref->{$gene_id}->{ti} = $trans_id;
            $hits_ref->{$gene_id}->{hs} = $hit_strand;
	  }else{
	    $redundant_hits++;
	    $gt_redundant_hits++;
	    next(); #Skip less significant hits to the same sequence
	  }

	}else{
	  #This hit is to a new gene ID
	  $hits_ref->{$gene_id}->{pi} = $percent_identity;
	  $hits_ref->{$gene_id}->{al} = $alignment_length;
	  $hits_ref->{$gene_id}->{ss} = $subject_start;
	  $hits_ref->{$gene_id}->{se} = $subject_end;
	  $hits_ref->{$gene_id}->{bs} = $bit_score;
	  $hits_ref->{$gene_id}->{ti} = $trans_id;
          $hits_ref->{$gene_id}->{hs} = $hit_strand;
	  $blast_results_stored++;
	  $gt_blast_results_stored++;
	}

      }else{
	#create the first hit record for this read
	my %hits;
	$hits{$gene_id}{pi} = $percent_identity;
	$hits{$gene_id}{al} = $alignment_length;
	$hits{$gene_id}{ss} = $subject_start;
	$hits{$gene_id}{se} = $subject_end;
	$hits{$gene_id}{bs} = $bit_score;
	$hits{$gene_id}{ti} = $trans_id;
        $hits{$gene_id}{hs} = $hit_strand;
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
    print BLUE, "\n\tFiltered out $filtered_pseudogene_hits of these hits because they were hits to pseudogenes", RESET;
    print BLUE, "\n\tFound $redundant_hits of equal or lesser quality to single gene sequences (i.e. caused by multiple transcripts mostly)", RESET;
    print BLUE, "\n\tStored a total of $blast_results_stored blast results corresponding to $read_ids_processed reads (only the top hit will be printed)\n", RESET;

    print BLUE, "\n\tFound hits for ($percent_reads_with_hits%) of $reads_analyzed_per_block reads in this block of reads (block size specified by user)", RESET;
    print BLUE, "\n\tReads which could be assigned to a gene as a PAIR:", RESET;
    print BLUE, "\n\t\tFound a total of $unambig_readpair_reads ($percent_reads_unambig_as_pair%) reads belonging to unambiguous read-pair mappings", RESET;
    print BLUE, "\n\t\tFound a total of $ambig_readpair_reads ($percent_reads_ambig_as_pair%) reads belonging to ambiguous read-pair mappings", RESET;
    print BLUE, "\n\tReads which had to be assigned to a gene INDIVIDUALLY:", RESET;
    print BLUE, "\n\t\tFound a total of $unambig_individual_reads ($percent_unambig_individual_reads%) individual reads mapping unambiguously", RESET;
    print BLUE, "\n\t\tFound a total of $ambig_individual_reads mapping ($percent_ambig_individual_reads%) ambiguously", RESET;

    print LOG "\n\tProcessed $blast_results_processed blast results";
    print LOG "\n\tFiltered out $filtered_pseudogene_hits of these hits because they were hits to pseudogenes";
    print LOG "\n\tFound $redundant_hits of equal or lesser quality to single gene sequences (i.e. caused by multiple transcripts mostly)";
    print LOG "\n\tStored a total of $blast_results_stored blast results corresponding to $read_ids_processed reads (only the top hit will be printed)\n";

    print LOG "\n\tFound hits for ($percent_reads_with_hits%) of $reads_analyzed_per_block reads in this block of reads (block size specified by user)";
    print LOG "\n\tReads which could be assigned to a gene as a PAIR:";
    print LOG "\n\t\tFound a total of $unambig_readpair_reads ($percent_reads_unambig_as_pair%) reads belonging to unambiguous read-pair mappings";
    print LOG "\n\t\tFound a total of $ambig_readpair_reads ($percent_reads_ambig_as_pair%) reads belonging to ambiguous read-pair mappings";
    print LOG "\n\tReads which had to be assigned to a gene INDIVIDUALLY:";
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
  print BLUE, "\n\tFiltered out $gt_filtered_pseudogene_hits of these hits because they were hits to pseudogenes", RESET;
  print BLUE, "\n\tFound $gt_redundant_hits of equal or lesser quality to single gene sequences (i.e. caused by multiple transcripts mostly)", RESET;
  print BLUE, "\n\tStored a total of $gt_blast_results_stored blast results corresponding to $gt_read_ids_processed reads (only the top hit will be printed)\n", RESET;

  print BLUE, "\n\tFound hits for ($gt_percent_reads_with_hits%) of $total_reads_analyzed total reads analyzed (specified by user)", RESET;
  print BLUE, "\n\tReads which could be assigned to a gene as a PAIR:", RESET;
  print BLUE, "\n\t\tFound a total of $gt_unambig_readpair_reads ($gt_percent_reads_unambig_as_pair%) reads belonging to unambiguous read-pair mappings", RESET;
  print BLUE, "\n\t\tFound a total of $gt_ambig_readpair_reads ($gt_percent_reads_ambig_as_pair%) reads belonging to ambiguous read-pair mappings", RESET;
  print BLUE, "\n\tReads which had to be assigned to a gene INDIVIDUALLY:", RESET;
  print BLUE, "\n\t\tFound a total of $gt_unambig_individual_reads ($gt_percent_unambig_individual_reads%) individual reads mapping unambiguously", RESET;
  print BLUE, "\n\t\tFound a total of $gt_ambig_individual_reads ($gt_percent_ambig_individual_reads%) reads mapping ambiguously", RESET;
  print BLUE, "\n\tNOTE: Because of cases were either R1 or R2 has NO hits, these numbers will not correspond exactly to the number of lines in your file!", RESET;


  print LOG "\n\nGRAND TOTAL STATS";
  print LOG "\n\tProcessed $gt_blast_results_processed blast results";
  print LOG "\n\tFiltered out $gt_filtered_pseudogene_hits of these hits because they were hits to pseudogenes";
  print LOG "\n\tFound $gt_redundant_hits of equal or lesser quality to single gene sequences (i.e. caused by multiple transcripts mostly)";
  print LOG "\n\tStored a total of $gt_blast_results_stored blast results corresponding to $gt_read_ids_processed reads (only the top hit will be printed)\n";

  print LOG "\n\tFound hits for ($gt_percent_reads_with_hits%) of $total_reads_analyzed total reads analyzed (specified by user)";
  print LOG "\n\tReads which could be assigned to a gene as a PAIR:";
  print LOG "\n\t\tFound a total of $gt_unambig_readpair_reads ($gt_percent_reads_unambig_as_pair%) reads belonging to unambiguous read-pair mappings";
  print LOG "\n\t\tFound a total of $gt_ambig_readpair_reads ($gt_percent_reads_ambig_as_pair%) reads belonging to ambiguous read-pair mappings";
  print LOG "\n\tReads which had to be assigned to a gene INDIVIDUALLY:";
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

  #First attempt to map both reads to a single gene - only if this fails will they be treated as individual reads

  #Generate a complete list of genes hit by EITHER read
  #Calculate a combined score for cases where both reads hit a single gene
  my %gene_hit_list;
  foreach my $read (keys %read_pair_blast_results){
    my $hits_ref = $read_pair_blast_results{$read}{blast_hits};

    foreach my $gene_id (keys %{$hits_ref}){
      if ($gene_hit_list{$gene_id}){
	$gene_hit_list{$gene_id}{combined_score} += $hits_ref->{$gene_id}->{bs};
	$gene_hit_list{$gene_id}{paired_hit} = "yes";
	if ($hits_ref->{$gene_id}->{bs} < $min_bit_score){
	  $gene_hit_list{$gene_id}{both_reads_pass} = "no";
	}
      }else{
	$gene_hit_list{$gene_id}{combined_score} = $hits_ref->{$gene_id}->{bs};
	$gene_hit_list{$gene_id}{paired_hit} = "no";
	if ($hits_ref->{$gene_id}->{bs} >= $min_bit_score){
	  $gene_hit_list{$gene_id}{both_reads_pass} = "yes";
	}else{
	  $gene_hit_list{$gene_id}{both_reads_pass} = "no";
	}
      }
    }
  }

  #Determine the rank of gene pair hits according to combined bit score - where the same gene was hit by both reads
  my $rank = 0;
  my %ranked_gene_pair_hits;
  foreach my $gene_id (sort {$gene_hit_list{$b}{combined_score}<=>$gene_hit_list{$a}{combined_score}} keys %gene_hit_list){
    unless (($gene_hit_list{$gene_id}{paired_hit} eq "yes") &&  ($gene_hit_list{$gene_id}{both_reads_pass} eq "yes")){
      next();
    }
    $rank++;
    $ranked_gene_pair_hits{$rank}{gene_id} = $gene_id;
    $ranked_gene_pair_hits{$rank}{combined_score} = $gene_hit_list{$gene_id}{combined_score};
  }

  my $read_pair_gene_count = keys %ranked_gene_pair_hits;

  #Possible outcomes...
  #A.) There is a single gene where both reads map to the same gene, and both have bit scores greater than some cutoff ($min_bit_score)
  #    - Print out both reads of the pair using the hits to this single target gene and return().

  #B.) There are multiple genes where both reads map to a single gene, and both reads have bit scores larger than some cutoff ($min_bit_score)
  #   - If one of these genes has a clear 'best' combined bit score (simply the highest? or at least x larger than the next best?)
  #     - Print out both reads of the pair using the hits to this single target gene and return().
  #   - If multiple genes have the same combined score, then the hit is ambiguous.
  #     - Mark both reads as ambiguous, print them out and return().

  if ($read_pair_gene_count == 1){
    #Only one gene found with quality hits to both reads
    my $gene_id = $ranked_gene_pair_hits{1}{gene_id};

    my $read1_hits_ref = $read_pair_blast_results{$read1_id}{blast_hits};
    my $read2_hits_ref = $read_pair_blast_results{$read2_id}{blast_hits};

    #DEBUG
    unless ($read1_hits_ref->{$gene_id} && $read2_hits_ref->{$gene_id}){
      print RED, "\n$read_pair_id\tSingle Paired Read - Unambiguous.  Gene ID: $gene_id\tREAD1: $read1_id\tREAD2: $read2_id", RESET;
      #print Dumper $read1_hits_ref;
      #print Dumper $read2_hits_ref;
      print Dumper %read_pair_blast_results;
      print Dumper %gene_hit_list;
      print Dumper %ranked_gene_pair_hits;
      exit();
    }

    #    - Read_ID, Hit_Type (i.e. Ambiguous or Top Hit), Gene_ID, Trans_ID (best hit),  AlignLength, Percent Identity, BitScore, Subject Start, Subject End
    $top_hits_db{$read1_id} = "Top_Hit\t$gene_id\t$read1_hits_ref->{$gene_id}->{ti}\t$read1_hits_ref->{$gene_id}->{al}\t$read1_hits_ref->{$gene_id}->{pi}\t$read1_hits_ref->{$gene_id}->{bs}\t$read1_hits_ref->{$gene_id}->{ss}\t$read1_hits_ref->{$gene_id}->{se}\t$read1_hits_ref->{$gene_id}->{hs}";

    $top_hits_db{$read2_id} = "Top_Hit\t$gene_id\t$read2_hits_ref->{$gene_id}->{ti}\t$read2_hits_ref->{$gene_id}->{al}\t$read2_hits_ref->{$gene_id}->{pi}\t$read2_hits_ref->{$gene_id}->{bs}\t$read2_hits_ref->{$gene_id}->{ss}\t$read2_hits_ref->{$gene_id}->{se}\t$read2_hits_ref->{$gene_id}->{hs}";

    $unambig_readpair_reads += 2;
    $gt_unambig_readpair_reads += 2;

    #print YELLOW, "\n\tSingle gene - two quality hits", RESET;
    return();

  }elsif($read_pair_gene_count > 1){

    #Multiple genes found that have quality hits to both reads.  See if there is a clear winner according to combined score

    my $read1_hits_ref = $read_pair_blast_results{$read1_id}{blast_hits};
    my $read2_hits_ref = $read_pair_blast_results{$read2_id}{blast_hits};

    if ($ranked_gene_pair_hits{1}{combined_score} > $ranked_gene_pair_hits{2}{combined_score}){

      my $gene_id = $ranked_gene_pair_hits{1}{gene_id};

      #DEBUG
      unless ($read1_hits_ref->{$gene_id} && $read2_hits_ref->{$gene_id}){
	print RED, "\nMultiple Paired Reads - Unambiguous.  Gene ID: $gene_id\tREAD1: $read1_id\tREAD2: $read2_id", RESET;
	print Dumper $read1_hits_ref;
	print Dumper $read2_hits_ref;
	exit();
      }

      $top_hits_db{$read1_id} = "Top_Hit\t$gene_id\t$read1_hits_ref->{$gene_id}->{ti}\t$read1_hits_ref->{$gene_id}->{al}\t$read1_hits_ref->{$gene_id}->{pi}\t$read1_hits_ref->{$gene_id}->{bs}\t$read1_hits_ref->{$gene_id}->{ss}\t$read1_hits_ref->{$gene_id}->{se}\t$read1_hits_ref->{$gene_id}->{hs}";
      $top_hits_db{$read2_id} = "Top_Hit\t$gene_id\t$read2_hits_ref->{$gene_id}->{ti}\t$read2_hits_ref->{$gene_id}->{al}\t$read2_hits_ref->{$gene_id}->{pi}\t$read2_hits_ref->{$gene_id}->{bs}\t$read2_hits_ref->{$gene_id}->{ss}\t$read2_hits_ref->{$gene_id}->{se}\t$read2_hits_ref->{$gene_id}->{hs}";

      $unambig_readpair_reads += 2;
      $gt_unambig_readpair_reads += 2;

      #print YELLOW, "\n\tMultiple genes - two quality hits - but a 'best' one was found", RESET;
      return();
    }else{

      my $gene_id = $ranked_gene_pair_hits{1}{gene_id};

      #DEBUG
      unless ($read1_hits_ref->{$gene_id} && $read2_hits_ref->{$gene_id}){
	print RED, "\nPaired Read - Ambiguous.  Gene ID: $gene_id\tREAD1: $read1_id\tREAD2: $read2_id", RESET;
	print Dumper $read1_hits_ref;
	print Dumper $read2_hits_ref;
	exit();
      }

      $top_hits_db{$read1_id} = "Ambiguous\tNA\tNA\t$read1_hits_ref->{$gene_id}->{al}\t$read1_hits_ref->{$gene_id}->{pi}\t$read1_hits_ref->{$gene_id}->{bs}\tNA\tNA\tNA";
      $top_hits_db{$read2_id} = "Ambiguous\tNA\tNA\t$read2_hits_ref->{$gene_id}->{al}\t$read2_hits_ref->{$gene_id}->{pi}\t$read2_hits_ref->{$gene_id}->{bs}\tNA\tNA\tNA";

      $ambig_readpair_reads += 2;
      $gt_ambig_readpair_reads += 2;

      #print YELLOW, "\n\tMultiple genes - two quality hits - ambiguous", RESET;
      return();
    }
  }

  #C.) There are no genes where both reads map to a single gene.  OR the BIT SCORE of one or both of these reads is below the cutoff.
  #   - Treat these as individual reads.  Use the bit score to identify the 'top hit'.  If there are multiple equal top hits, mark the read as Ambiguous
  #   - Once both reads have been processed, return().

  #print YELLOW, "\n\tNo Single Gene for Read Pair - Treating reads individually", RESET;

  foreach my $working_read_id (sort {$read_pair_blast_results{$a}{read_num} <=> $read_pair_blast_results{$b}{read_num}} keys %read_pair_blast_results){

    #DEBUG
    #unless ($read_pair_blast_results{$working_read_id}{read_num} =~ /\d+/){
    #  print Dumper %read_pair_blast_results;
    #  exit();
    #}

    my $hits_ref = $read_pair_blast_results{$working_read_id}{blast_hits};
    my $hit_count = keys %{$hits_ref};

    #Create a sorted list to help determine the ranking of hits to each gene according to Bit Scores
    my %sorted_list;
    foreach my $gene_id (sort keys %{$hits_ref}){
      my $bs = $hits_ref->{$gene_id}->{bs};

      if ($sorted_list{$bs}){

	my $list_ref = $sorted_list{$bs}{gene_list};
	push (@{$list_ref}, $gene_id);
      }else{
	my @gene_list;
	push (@gene_list, $gene_id);
	$sorted_list{$bs}{gene_list} = \@gene_list;
      }
    }

  BS:foreach my $bs (sort {$b <=> $a} keys %sorted_list){
      my @gene_list = @{$sorted_list{$bs}{gene_list}};

      #D.) Print the best hits for individual reads to the 'top_hits' file

      my $top_gene_count = scalar(@gene_list);
      my $gene_id = $gene_list[0];

      if ($top_gene_count == 1){

	#DEBUG
	unless ($hits_ref->{$gene_id}){
	  print RED, "\nSingle Read - Unambiguous.  Gene ID: $gene_id\tREAD: $working_read_id", RESET;
	  print Dumper $hits_ref;
	  exit();
	}


	$top_hits_db{$working_read_id} = "Top_Hit\t$gene_id\t$hits_ref->{$gene_id}->{ti}\t$hits_ref->{$gene_id}->{al}\t$hits_ref->{$gene_id}->{pi}\t$hits_ref->{$gene_id}->{bs}\t$hits_ref->{$gene_id}->{ss}\t$hits_ref->{$gene_id}->{se}\t$hits_ref->{$gene_id}->{hs}";
	$unambig_individual_reads++;
	$gt_unambig_individual_reads++;
	last (BS);

      }else{

	#DEBUG
	unless ($hits_ref->{$gene_id}){
	  print RED, "\nSingle Read - Ambiguous.  Gene ID: $gene_id\tREAD: $working_read_id", RESET;
	  print Dumper $hits_ref;
	  exit();
	}

	$top_hits_db{$working_read_id} = "Ambiguous\tNA\tNA\t$hits_ref->{$gene_id}->{al}\t$hits_ref->{$gene_id}->{pi}\t$hits_ref->{$gene_id}->{bs}\tNA\tNA\tNA";
	$ambig_individual_reads++;
	$gt_ambig_individual_reads++;
	last (BS);

      }
    }
  }
  return();
}


#############################################################################################################################
#3.) Calculate the transcript position of the hit
#############################################################################################################################
sub determineTranscriptPositions{
  my %args = @_;

  #    - Parse the top_hits and append this info to each line
  print BLUE, "\n\nParsing top hits berkley DB - calculating relative position of top hits and adding to the top hits db\n\n", RESET;
  print LOG "\n\nParsing top hits berkley DB - calculating relative position of top hits and adding to the top hits db\n\n";

  #    - Take the centre of the query hit coordinates
  #    - Determine the relative position of this hit in the TRANSCRIPT that it mapped to (i.e. percent distance from 3' end)
  #    - Only do this for 'Top Hits', this value will be NA for Ambiguous hits
  #    - For convenience, also write out the transcript size for each hit
  #    - The new columns will be called 'TranscriptSize' and 'RelativePosition'


  #Format of top hits db records: 
  # "Top_Hit\t$gene_id\t$hits_ref->{$gene_id}->{ti}\t$hits_ref->{$gene_id}->{al}\t$hits_ref->{$gene_id}->{pi}\t$hits_ref->{$gene_id}->{bs}\t$hits_ref->{$gene_id}->{ss}\t$hits_ref->{$gene_id}->{se}\t$hits_ref->{$gene_id}->{hs}";

  my $counter = 0;

  while (my ($read_id, $string) = each %top_hits_db){
    my $counter++;

    if ($counter == 100000){
      $counter = 0;
      $| = 1; print BLUE, ".", RESET; $| = 0;
    }

    chomp($string);
    my @line = split("\t", $string);
    my $hit_type = $line[0];

    if ($hit_type eq "Ambiguous"){
      $string = "$string"."\tNA\tNA";
      $top_hits_db{$read_id} = "$string";
      next();
    }

    my $trans_id = $line[2];
    my $subject_start = $line[6];
    my $subject_end = $line[7];

    my $transcript_size = $transcript_gene_map{$trans_id}{transcript_size};
    my $centre = $subject_end - (($subject_end-$subject_start)/2);
    my $relative_position = ($centre/$transcript_size)*100;
    my $relative_position_f = sprintf("%.1f", $relative_position);

    $string = "$string"."\t$transcript_size\t$relative_position_f";
    $top_hits_db{$read_id} = "$string";
  }

  return();
}


#############################################################################################################################
#4.) At this point, also calculate the CHROMOSOME coordinates (gapped) of each top hit
#############################################################################################################################
sub getHitChromosomeCoordinates{
  my %args = @_;

  print BLUE, "\n\nCalculating chromosome coordinates for all top hits", RESET;
  print LOG "\n\nCalculating chromosome coordinates for all top hits";

  #- Once again this can be done while parsing the top_hits file.
  #- Note that it is possible for one or both reads to span an exon junction.
  #- So the transcript coordinates from mapping to EnsEMBL transcripts have to be mapped back to the genome by taking the introns into account
  #- This is similar to display protein feature from EnsEMBL on the UCSC genome browser (these are also reported relative to the transcript)
  #- Append two columns to the top_hits file: start_coords (array of chromosome start positions) and end_coord (array of chromosome end positions)

  #A.) First step is to get chromosome coordinates for all exons of all genes
  print BLUE, "\nFirst calculating chromosome coordinates for all exons of all genes\n", RESET;
  print LOG "\nFirst calculating chromosome coordinates for all exons of all genes\n";

  my $counter = 0;
  foreach my $gene_id (sort keys %{$gene_transcripts_ref}){
    my $transcripts_ref = $gene_transcripts_ref->{$gene_id}->{transcripts};

    my $chromosome = $genes_ref->{$gene_id}->{chromosome};
    my $chr_strand = $genes_ref->{$gene_id}->{chr_strand};
    my $chr_start = $genes_ref->{$gene_id}->{chr_start};
    my $chr_end = $genes_ref->{$gene_id}->{chr_end};
    my $gene_start = $genes_ref->{$gene_id}->{gene_start};
    my $gene_end = $genes_ref->{$gene_id}->{gene_end};

    foreach my $trans_id (keys %{$transcripts_ref}){
      $transcript_gene_map{$trans_id}{gene_id} = $gene_id;

      #Calculate the size of each transcript by adding up the size of its exons
      my $size = 0;
      my $exons_ref = $transcripts_ref->{$trans_id}->{exons};

      foreach my $exon_id (keys %{$exons_ref}){
        $counter++;
        if ($counter == 1000){
          $counter = 0;
          $| = 1; print BLUE, ".", RESET; $| = 0;
        }

	my $start = $exons_ref->{$exon_id}->{exon_start};
	my $end = $exons_ref->{$exon_id}->{exon_end};

	#Make sure the supplied coordinates are actually within the specified gene
	unless ($start >= $gene_start-1 && $start <= $gene_end+1){
	  print RED, "\nStart coordinate ($start) provided to convertGeneCoordinates() does not appear valid for gene_id $gene_id\n\n", RESET;
	  exit();
	}
	unless ($end >= $gene_start-1 && $end <= $gene_end+1){
	  print RED, "\nEnd coordinate ($end) provided to convertGeneCoordinates() does not appear valid for gene_id $gene_id\n\n", RESET;
	  exit();
	}

	#Convert provided gene coordinates to coordinates relative to the chromosome
	if ($chr_strand == 1){
	  my $query_chr_start = $chr_start + $start - 1;
	  my $query_chr_end = $chr_start + $end - 1;

	  #Make sure the start and end are reported such that start is always smaller than end
	  my $temp;
	  if ($query_chr_start > $query_chr_end){
	    $temp = $query_chr_start;
	    $query_chr_start = $query_chr_end;
	    $query_chr_end = $temp;
	  }

	  #print YELLOW, "\n$gene_id\t$chromosome\texon: $exon_id\tstart: $query_chr_start\tend: $query_chr_end\tstrand: +", RESET;

	  $exons_ref->{$exon_id}->{chr_start} = $query_chr_start;
	  $exons_ref->{$exon_id}->{chr_end} = $query_chr_end;
	  $exons_ref->{$exon_id}->{strand} = "+";
	  $exons_ref->{$exon_id}->{size} = ($query_chr_end - $query_chr_start)+1;

	}elsif ($chr_strand == -1){

	  my $query_chr_start = $chr_end - $end + 1;
	  my $query_chr_end = $chr_end - $start + 1;

	  #Make sure the start and end are reported such that start is always smaller than end
	  my $temp;
	  if ($query_chr_start > $query_chr_end){
	    $temp = $query_chr_start;
	    $query_chr_start = $query_chr_end;
	    $query_chr_end = $temp;
	  }

	  #print YELLOW, "\n$gene_id\t$chromosome\texon: $exon_id\tstart: $query_chr_start\tend: $query_chr_end\tstrand: -", RESET;

	  $exons_ref->{$exon_id}->{chr_start} = $query_chr_start;
	  $exons_ref->{$exon_id}->{chr_end} = $query_chr_end;
	  $exons_ref->{$exon_id}->{strand} = "-";
	  $exons_ref->{$exon_id}->{size} = ($query_chr_end - $query_chr_start)+1;

	}else{
	  print RED, "\nStrand format: $chr_strand not understood !\n\n", RESET;
	  exit();
	}
      }
    }
  }


  #B.) At this time, also get the chromosome coordinates for all exon-content coordinates
  print BLUE, "\nNow calculating chromosome coordinates for the EXON CONTENT of each gene\n", RESET;
  print LOG "\nNow calculating chromosome coordinates for the EXON CONTENT of each gene\n";
  $counter = 0;
  foreach my $gene_id (sort keys %{$gene_exon_content_ref}){

    my $chromosome = $genes_ref->{$gene_id}->{chromosome};
    my $chr_strand = $genes_ref->{$gene_id}->{chr_strand};
    my $chr_start = $genes_ref->{$gene_id}->{chr_start};
    my $chr_end = $genes_ref->{$gene_id}->{chr_end};
    my $gene_start = $genes_ref->{$gene_id}->{gene_start};
    my $gene_end = $genes_ref->{$gene_id}->{gene_end};

    #Calculate the size of each transcript by adding up the size of its exons
    my $size = 0;
    my $exon_content_ref = $gene_exon_content_ref->{$gene_id}->{exon_content};

    foreach my $exon_id (keys %{$exon_content_ref}){

      $counter++;
      if ($counter == 1000){
        $counter = 0;
        $| = 1; print BLUE, ".", RESET; $| = 0;
      }

      my $start = $exon_content_ref->{$exon_id}->{start};
      my $end = $exon_content_ref->{$exon_id}->{end};

      #Make sure the supplied coordinates are actually within the specified gene
      unless ($start >= $gene_start-1 && $start <= $gene_end+1){
	print RED, "\nStart coordinate ($start) provided to convertGeneCoordinates() does not appear valid for gene_id $gene_id\n\n", RESET;
	exit();
      }
      unless ($end >= $gene_start-1 && $end <= $gene_end+1){
	print RED, "\nEnd coordinate ($end) provided to convertGeneCoordinates() does not appear valid for gene_id $gene_id\n\n", RESET;
	exit();
      }

      #Convert provided gene coordinates to coordinates relative to the chromosome
      if ($chr_strand == 1){
	my $query_chr_start = $chr_start + $start - 1;
	my $query_chr_end = $chr_start + $end - 1;

	#Make sure the start and end are reported such that start is always smaller than end
	my $temp;
	if ($query_chr_start > $query_chr_end){
	  $temp = $query_chr_start;
	  $query_chr_start = $query_chr_end;
	  $query_chr_end = $temp;
	}

	#print YELLOW, "\n$gene_id\t$chromosome\texon: $exon_id\tstart: $query_chr_start\tend: $query_chr_end\tstrand: +", RESET;

	$exon_content_ref->{$exon_id}->{chr_start} = $query_chr_start;
	$exon_content_ref->{$exon_id}->{chr_end} = $query_chr_end;
	$exon_content_ref->{$exon_id}->{strand} = "+";
	$exon_content_ref->{$exon_id}->{size} = ($query_chr_end - $query_chr_start)+1;

      }elsif ($chr_strand == -1){

	my $query_chr_start = $chr_end - $end + 1;
	my $query_chr_end = $chr_end - $start + 1;

	#Make sure the start and end are reported such that start is always smaller than end
	my $temp;
	if ($query_chr_start > $query_chr_end){
	  $temp = $query_chr_start;
	  $query_chr_start = $query_chr_end;
	  $query_chr_end = $temp;
	}

	#print YELLOW, "\n$gene_id\t$chromosome\texon: $exon_id\tstart: $query_chr_start\tend: $query_chr_end\tstrand: -", RESET;

	$exon_content_ref->{$exon_id}->{chr_start} = $query_chr_start;
	$exon_content_ref->{$exon_id}->{chr_end} = $query_chr_end;
	$exon_content_ref->{$exon_id}->{strand} = "-";
	$exon_content_ref->{$exon_id}->{size} = ($query_chr_end - $query_chr_start)+1;

      }else{
	print RED, "\nStrand format: $chr_strand not understood !\n\n", RESET;
	exit();
      }
    }
  }


  #C.) Parse top hits db and for each hit with a valid top hit calculate the chromosome coordinates of the hit
  #    - For Ambiguous hits, print NA values
  print BLUE, "\nNow calculating chromosome coordinates for all read hits\n", RESET;
  print LOG "\nNow calculating chromosome coordinates for all read hits\n";

  #First we need to take the position of the hit within its top transcript (subject start & end) and determine the corresponding chromosome coordinates

  #Format of top hits db records: 
  # "Top_Hit\t$gene_id\t$hits_ref->{$gene_id}->{ti}\t$hits_ref->{$gene_id}->{al}\t$hits_ref->{$gene_id}->{pi}\t$hits_ref->{$gene_id}->{bs}\t$hits_ref->{$gene_id}->{ss}\t$hits_ref->{$gene_id}->{se}\t$hits_ref->{$gene_id}->{hs}\t$transcript_size\t$relative_position_f"";

  $counter = 0;
  while (my ($read_id, $string) = each %top_hits_db){

    $counter++;
    if ($counter == 100000){
      $counter = 0;
      $| = 1; print BLUE, ".", RESET; $| = 0;
    }

    chomp($string);
    my @line = split("\t", $string);

    #Print NA values for Ambiguous hits
    my $hit_type = $line[0];
    if ($hit_type eq "Ambiguous"){
      $string = "$string"."\tNA\tNA\tNA";
      $top_hits_db{$read_id} = $string;
      next();
    }

    my $gene_id = $line[1];
    my $trans_id = $line[2];
    my $align_length = $line[3];
    my $subject_start = $line[6];
    my $subject_end = $line[7];

    #Get the corresponding, Gene/Transcript/Exons records for the gene hit by this read
    my $transcripts_ref = $gene_transcripts_ref->{$gene_id}->{transcripts};
    my $exons_ref = $transcripts_ref->{$trans_id}->{exons};
    my $chromosome = $genes_ref->{$gene_id}->{chromosome};
    my $chr_strand = $genes_ref->{$gene_id}->{chr_strand};

    my $ss_chromosome_coord;
    my $se_chromosome_coord;

    if ($chr_strand eq "1"){
      #Reads hitting a transcript from a +ve strand gene

      #Go through the exons (in order) and find the coordinate of the start and end position of the read hit within the transcript
      my $ss_working = $subject_start;
      foreach my $e_id (sort {$exons_ref->{$a}->{chr_start} <=> $exons_ref->{$b}->{chr_start}} keys %{$exons_ref}){
	$ss_working -= $exons_ref->{$e_id}->{size};
	if ($ss_working <= 0){
	  $ss_chromosome_coord = $exons_ref->{$e_id}->{chr_end} + $ss_working;
	  last();
	}
      }

      my $se_working = $subject_end;
      foreach my $e_id (sort {$exons_ref->{$a}->{chr_start} <=> $exons_ref->{$b}->{chr_start}} keys %{$exons_ref}){
	$se_working -= $exons_ref->{$e_id}->{size};
	if ($se_working <= 0){
	  $se_chromosome_coord = $exons_ref->{$e_id}->{chr_end} + $se_working;
	  last();
	}
      }

    }else{
      #Reads hitting a transcript from a -ve strand gene

      #Go through the exons (in order) and find the coordinate of the start and end position of the read hit within the transcript
      my $ss_working = $subject_start;
      foreach my $e_id (sort {$exons_ref->{$b}->{chr_start} <=> $exons_ref->{$a}->{chr_start}} keys %{$exons_ref}){
	$ss_working -= $exons_ref->{$e_id}->{size};
	if ($ss_working <= 0){
	  $ss_chromosome_coord = $exons_ref->{$e_id}->{chr_start} - $ss_working;
	  last();
	}
      }

      my $se_working = $subject_end;
      foreach my $e_id (sort {$exons_ref->{$b}->{chr_start} <=> $exons_ref->{$a}->{chr_start}} keys %{$exons_ref}){
	$se_working -= $exons_ref->{$e_id}->{size};
	if ($se_working <= 0){
	  $se_chromosome_coord = $exons_ref->{$e_id}->{chr_start} - $se_working;
	  last();
	}
      }
      #Once again make sure the smaller coordinate is reported as the start position
      if ($ss_chromosome_coord > $se_chromosome_coord){
	my $temp = $ss_chromosome_coord;
	$ss_chromosome_coord = $se_chromosome_coord;
	$se_chromosome_coord = $temp;
      }
    }

    #print YELLOW, "\nDEBUG: $read_id has length $align_length and maps to: chr$chromosome:$ss_chromosome_coord-$se_chromosome_coord ($chr_strand)", RESET;

    #Finally, we need to resolve cases where the hit is crossing one or more exon junctions...
    #Once this is done save the read 'starts' and 'ends' relative to the chromosome in two hashes
    my @read_starts;
    my @read_ends;

    #Figure out which exon(s) the hit start and end occur in
    my $first_exon_id;
    my $last_exon_id;
    foreach my $e_id (sort {$exons_ref->{$a}->{chr_start} <=> $exons_ref->{$b}->{chr_start}} keys %{$exons_ref}){
      if (($ss_chromosome_coord >= $exons_ref->{$e_id}->{chr_start}) && ($ss_chromosome_coord <= $exons_ref->{$e_id}->{chr_end})){
	$first_exon_id = $e_id;
      }
      if (($se_chromosome_coord >= $exons_ref->{$e_id}->{chr_start}) && ($se_chromosome_coord <= $exons_ref->{$e_id}->{chr_end})){
	$last_exon_id = $e_id;
      }
    }

    if ($first_exon_id == $last_exon_id){
      #If the hit start and end are both within a single exon then we already know the read_starts and read_ends
      push(@read_starts, $ss_chromosome_coord);
      push (@read_ends, $se_chromosome_coord);
    }else{

      #Build a list of the exons spanned
      my %exons_spanned;
      my $first_exon_found = 0;
      foreach my $e_id (sort {$exons_ref->{$a}->{chr_start} <=> $exons_ref->{$b}->{chr_start}} keys %{$exons_ref}){
	if ($e_id == $first_exon_id){
	  $exons_spanned{$e_id}{chr_start} = $exons_ref->{$e_id}->{chr_start};
	  $exons_spanned{$e_id}{chr_end} = $exons_ref->{$e_id}->{chr_end};
	  $first_exon_found = 1;
	  next();
	}
	if ($e_id == $last_exon_id){
	  $exons_spanned{$e_id}{chr_start} = $exons_ref->{$e_id}->{chr_start};
	  $exons_spanned{$e_id}{chr_end} = $exons_ref->{$e_id}->{chr_end};
	  last();
	}
	if ($first_exon_found == 1){
	  $exons_spanned{$e_id}{chr_start} = $exons_ref->{$e_id}->{chr_start};
	  $exons_spanned{$e_id}{chr_end} = $exons_ref->{$e_id}->{chr_end};
	}
      }

      #Now go through the list of exon_spanned and record the read_starts and ends
      push(@read_starts, $ss_chromosome_coord);
      foreach my $e_id (sort {$exons_spanned{$a}->{chr_start} <=> $exons_spanned{$b}->{chr_start}} keys %exons_spanned){
	if ($e_id == $first_exon_id){
	  push(@read_ends, $exons_spanned{$e_id}{chr_end});
	}elsif ($e_id == $last_exon_id){
	  push(@read_starts, $exons_spanned{$e_id}{chr_start})
	}else{
	  push(@read_starts, $exons_spanned{$e_id}{chr_start});
	  push(@read_ends, $exons_spanned{$e_id}{chr_end});
	}
      }
      push (@read_ends, $se_chromosome_coord);
    }
    #print YELLOW, "\n\tStarts = @read_starts\tEnds = @read_ends", RESET;
    $string = "$string"."\tchr$chromosome\t@read_starts\t@read_ends";
    $top_hits_db{$read_id} = $string;

  }
  return();
}


#############################################################################################################################
#5.) Create a summary of top hits and associated info with both reads of a pair on a single line
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
  tie(%read_pairs_db, 'BerkeleyDB::Btree', -Cachesize => 256000000, -Filename=> $read_pairs_db_file, -Flags => DB_CREATE) or die "can't open file $read_pairs_db_file: $! $BerkeleyDB::Error\n";

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
      $| = 1; print BLUE, ".", RESET; $| = 0;
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

  print SUMMARY "Read_ID\tDistanceBetweenReads_Genomic\tDistanceBetweenReads_Transcript\tR1_ID\tR1_HitType\tR1_GeneID\tR1_GeneName\tR1_AlignmentLength\tR1_PercentIdentity\tR1_BitScore\tR1_TranscriptSize\tR1_RelativePosition\tR1_Chromosome\tR1_Strand\tR1_ChrStartCoords\tR1_ChrEndCoords\tR2_ID\tR2_HitType\tR2_GeneID\tR2_GeneName\tR2_AlignmentLength\tR2_PercentIdentity\tR2_BitScore\tR2_TranscriptSize\tR2_RelativePosition\tR2_Chromosome\tR2_Strand\tR2_ChrStartCoords\tR2_ChrEndCoords\n";


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

    #Format of top hits db records: 
    # "Top_Hit\t$gene_id\t$hits_ref->{$gene_id}->{ti}\t$hits_ref->{$gene_id}->{al}\t$hits_ref->{$gene_id}->{pi}\t$hits_ref->{$gene_id}->{bs}\t$hits_ref->{$gene_id}->{ss}\t$hits_ref->{$gene_id}->{se}\t$hits_ref->{$gene_id}->{hs}\t$transcript_size\t$relative_position_f\tchr$chromosome\t@read_starts\t@read_ends";


    #Determine the fragment size based on the alignment of the two reads.  Only consider cases where both reads have 'Top_Hits'
    #Use the outer most coordinates to calculate the distance between reads in the genome and transcriptome (in this case should be the fragment size)
    my $distance_genomic = "NA";
    my $distance_transcript = "NA";

    if ($top_hits_db{$read1_id} && $top_hits_db{$read2_id}){
      my $string1 = $top_hits_db{$read1_id};
      my @string1 = split("\t", $string1);
      my $string2 = $top_hits_db{$read2_id};
      my @string2 = split("\t", $string2);

      my $type1 = $string1[0];
      my $type2 = $string2[0];
      my $gene1 = $string1[1];
      my $gene2 = $string2[1];
      my $chr1 = $string1[11];
      my $chr2 = $string2[11];
      my @starts1 = split(" ", $string1[12]);
      my @starts2 = split(" ", $string2[12]);
      my @ends1 = split(" ", $string1[13]);
      my @ends2 = split(" ", $string2[13]);

      #print RED, "DEBUG: \ns1: starts1\ns2: @starts2\ne1: @ends1\ne2: @ends2\nct: @coords_temp\ncs: @coords_sort\ngs: $grand_start\nge: $grand_end\n", RESET;
      #exit();

      if ($type1 eq "Top_Hit" && $type2 eq "Top_Hit"){

        #Calculate the distance between Read1 and Read2 in the GENOME and TRANSCRIPTOME
        my @coords_temp = (@starts1, @starts2, @ends1, @ends2);
        my @coords_sort = sort {$a <=> $b}(@coords_temp);
        my $grand_start = $coords_sort[0];
        my $grand_end = $coords_sort[(scalar(@coords_sort)-1)];

        #Calculate the distance between Read1 and Read2 in the GENOME
	if ($chr1 eq $chr2){
	  $distance_genomic = ($grand_end - $grand_start)+1;
	}else{
	  $distance_genomic = "Chromosome_Mismatch";
	}

        #Calculate the distance between Read1 and Read2 in the GENOME
	if ($gene1 == $gene2){
	  $distance_transcript = &calculateTranscriptDistance('-gene_id'=>$gene1, '-read1_start'=>$grand_start, '-read2_end'=>$grand_end);
	}else{
	  $distance_transcript = "Gene_Mismatch";
	}
      }
    }

    #Print out the read_id
    print SUMMARY "$read_pair_id\t$distance_genomic\t$distance_transcript";


    #READ 1#
    if ($top_hits_db{$read1_id}){
      my $string = $top_hits_db{$read1_id};
      my @string = split("\t", $string);

      my $hit_type = $string[0];
      my $subject_id = $string[1];
      my $align_length = $string[3];
      my $percent_id = $string[4];
      my $bit_score = $string[5];
      my $hit_strand = $string[8];
      my $trans_size = $string[9];
      my $position = $string[10];
      my $chr = $string[11];
      my $starts = $string[12];
      my $ends = $string[13];

      my $gene_name = "NA";
      unless ($subject_id eq "NA"){
        $gene_name = $genes_ref->{$subject_id}->{gene_name};
      }

      #print SUMMARY "\tR1_ID\tR1_HitType\tR1_GeneID\tR1_GeneName\tR1_AlignmentLength\tR1_PercentIdentity\tR1_BitScore\tR1_TranscriptSize\tR1_RelativePosition\tR1_Chromosome\tR1_Strand\tR1_ChrStartCoords\tR1_ChrEndCoords";


      print SUMMARY "\t$read1_id\t$hit_type\t$subject_id\t$gene_name\t$align_length\t$percent_id\t$bit_score\t$trans_size\t$position\t$chr\t$hit_strand\t$starts\t$ends";

      #If the upper bit score criteria is met and the status of this read is still 'Unassigned' in the read record file, change the status
      #Possible status is 'ENST_U' or 'ENST_A'

      #If a read has already been assigned to 'Repeat_U' or 'Repeat_A', start by resetting these to 'Unassigned'
      if ($read1_status eq "ENST_U" || $read1_status eq "ENST_A"){
	$read1_status = "Unassigned";
      }

      #Make sure that the status of this reads is still 'Unassigned'
      if ($read1_status eq "Unassigned"){

	#Now, if the top hit for this read meets the user specified threshold, change its status
	#- The resulting status will depend on whether the hit was ambiguous or not ...
	if (($hit_type eq "Top_Hit") && ($bit_score >= $upper_bit_score) && ($align_length >= $min_align_length)){
	  $enst_u_count++;
	  $read1_status = "ENST_U";
	}elsif(($hit_type eq "Ambiguous") && ($bit_score >= $upper_bit_score) && ($align_length >= $min_align_length)){
	  $enst_a_count++;
	  $read1_status = "ENST_A";
	}
      }
    }else{
      #Read 1 is not defined but READ 2 was
      print SUMMARY "\t$read1_id\tNone\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA";
    }

    #READ 2#
    if ($top_hits_db{$read2_id}){
      my $string = $top_hits_db{$read2_id};
      my @string = split("\t", $string);

      my $hit_type = $string[0];
      my $subject_id = $string[1];
      my $align_length = $string[3];
      my $percent_id = $string[4];
      my $bit_score = $string[5];
      my $hit_strand = $string[8];
      my $trans_size = $string[9];
      my $position = $string[10];
      my $chr = $string[11];
      my $starts = $string[12];
      my $ends = $string[13];

      my $gene_name = "NA";
      unless ($subject_id eq "NA"){
        $gene_name = $genes_ref->{$subject_id}->{gene_name};
      }


      #print SUMMARY "\tR2_ID\tR2_HitType\tR2_GeneID\tR2_GeneName\tR2_AlignmentLength\tR2_PercentIdentity\tR2_BitScore\tR2_TranscriptSize\tR2_RelativePosition\tR2_Chromosome\tR2_Chromosome\tR2_ChrStartCoords\tR2_ChrEndCoords\n";

      print SUMMARY "\t$read2_id\t$hit_type\t$subject_id\t$gene_name\t$align_length\t$percent_id\t$bit_score\t$trans_size\t$position\t$chr\t$hit_strand\t$starts\t$ends\n";

      #If the upper bit score criteria is met and the status of this read is still 'Unassigned' in the read record file, change the status
      #Possible status is 'Repeat_U' or 'Repeat_A'

      #If a read has already been assigned to 'Repeat_U' or 'Repeat_A', start by resetting these to 'Unassigned'
      if ($read2_status eq "ENST_U" || $read2_status eq "ENST_A"){
	$read2_status = "Unassigned";
      }

      #Make sure that the status of this reads is still 'Unassigned'
      if ($read2_status eq "Unassigned"){

	#Now, if the top hit for this read meets the user specified threshold, change its status
	#- The resulting status will depend on whether the hit was ambiguous or not ...
	if (($hit_type eq "Top_Hit") && ($bit_score >= $upper_bit_score) && ($align_length >= $min_align_length)){
	  $enst_u_count++;
	  $read2_status = "ENST_U";
	}elsif(($hit_type eq "Ambiguous") && ($bit_score >= $upper_bit_score) && ($align_length >= $min_align_length)){
	  $enst_a_count++;
	  $read2_status = "ENST_A";
	}
      }
    }else{
      #Read 2 is not defined but Read 1 was
      print SUMMARY "\t$read2_id\tNone\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n";

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
  print BLUE, "\n\nUpdate read status values in read record file: $read_records_infile\n", RESET;
  print LOG "\n\nUpdate read status values in read record file: $read_records_infile\n";

  my $new_read_record_file = "$base"."2";
  $header = 0;
  open (READS, "zcat $read_records_infile |") || die "\nCould not open read records file: $read_records_infile\n\n";
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

  my $percent_u = sprintf("%.2f", (($enst_u_count/$read_count)*100));
  my $percent_a = sprintf("%.2f", (($enst_a_count/$read_count)*100));

  print BLUE, "\n\nAssigned $enst_u_count reads ($percent_u%) to the status 'ENST_U', and $enst_a_count ($percent_a%) reads to the status 'ENST_A'\n\n", RESET;
  print LOG "\n\nAssigned $enst_u_count reads ($percent_u%) to the status 'ENST_U', and $enst_a_count ($percent_a%) reads to the status 'ENST_A'\n\n";

  print BLUE, "\n\nCompress updated read record file: $new_read_record_file\n", RESET;
  print LOG "\n\nCompress updated read record file: $new_read_record_file\n";

  my $cmd_gzip = "gzip -f $new_read_record_file";
  system ($cmd_gzip);

  my $cmd_mv = "mv $new_read_record_file".".gz"." $read_records_infile";
  print BLUE, "\n\nOverwriting existing read records file:\n\t$cmd_mv\n\n", RESET;
  print LOG "\n\nOverwriting existing read records file:\n\t$cmd_mv\n\n";
  system($cmd_mv);

  return();
}



#######################################################################################################################################
#Subroutine to help calculate the distance between reads within the transcript
#######################################################################################################################################
sub calculateTranscriptDistance{
  my %args = @_;
  my $gene_id = $args{'-gene_id'};
  my $read1_start = $args{'-read1_start'};
  my $read2_end = $args{'-read2_end'};

  my $transcript_distance = 0;

  #Figure out which exon content region(s) the read1 and read2 starts occur in
  my $read1_exon_id;
  my $read2_exon_id;

  my $exon_content_ref = $gene_exon_content_ref->{$gene_id}->{exon_content};

  foreach my $ec_id (sort {$exon_content_ref->{$a}->{chr_start} <=> $exon_content_ref->{$b}->{chr_start}} keys %{$exon_content_ref}){
    if (($read1_start >= $exon_content_ref->{$ec_id}->{chr_start}) && ($read1_start <= $exon_content_ref->{$ec_id}->{chr_end})){
      $read1_exon_id = $ec_id;
    }
    if (($read2_end >= $exon_content_ref->{$ec_id}->{chr_start}) && ($read2_end <= $exon_content_ref->{$ec_id}->{chr_end})){
      $read2_exon_id = $ec_id;
    }
  }

  #Now that the exon containing the start of read1 and read2 is known, calculate the amount of exon content between them
  if ($read1_exon_id == $read2_exon_id){

    #If the read1_start and read2_end are both within the same exon, then the distance between them is straight forward
    $transcript_distance = ($read2_end - $read1_start)+1;

  }else{

    #Build a list of the exons spanned
    my %exons_spanned;
    my $first_exon_found = 0;
    foreach my $ec_id (sort {$exon_content_ref->{$a}->{chr_start} <=> $exon_content_ref->{$b}->{chr_start}} keys %{$exon_content_ref}){
      if ($ec_id == $read1_exon_id){
	$exons_spanned{$ec_id}{chr_start} = $exon_content_ref->{$ec_id}->{chr_start};
	$exons_spanned{$ec_id}{chr_end} = $exon_content_ref->{$ec_id}->{chr_end};
	$first_exon_found = 1;
	next();
      }
      if ($ec_id == $read2_exon_id){
	$exons_spanned{$ec_id}{chr_start} = $exon_content_ref->{$ec_id}->{chr_start};
	$exons_spanned{$ec_id}{chr_end} = $exon_content_ref->{$ec_id}->{chr_end};
	last();
      }
      if ($first_exon_found == 1){
	$exons_spanned{$ec_id}{chr_start} = $exon_content_ref->{$ec_id}->{chr_start};
	$exons_spanned{$ec_id}{chr_end} = $exon_content_ref->{$ec_id}->{chr_end};
      }
    }

    #Now go through the list of exon_spanned and calculate the total transcript distance
    foreach my $e_id (sort {$exons_spanned{$a}->{chr_start} <=> $exons_spanned{$b}->{chr_start}} keys %exons_spanned){
      if ($e_id == $read1_exon_id){
	#Add distance from read1_start position to end of the exon it resides in
	my $exon_distance = ($exons_spanned{$e_id}{chr_end} - $read1_start)+1;
	$transcript_distance += $exon_distance;

      }elsif ($e_id == $read2_exon_id){
	#Add distance from read2_end position to beginning of the exon it resides in
	my $exon_distance = ($read2_end - $exons_spanned{$e_id}{chr_start})+1;
	$transcript_distance += $exon_distance;

      }else{
	#Add size of exon
	my $exon_distance = ($exons_spanned{$e_id}{chr_end} - $exons_spanned{$e_id}{chr_start})+1;
	$transcript_distance += $exon_distance;
      }
    }
  }
  return($transcript_distance);
}



