#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#The purpose of this script is to parse blast results files for Solexa Paired Reads that were BLASTed against a database of human Intergenic regions

#NOTE that these are the result of mapping individual reads of a pair to a database of intergenic regions (masked, non-genic regions defined by EnsEMBL annotation of genes)
#- Reads which hit a repeat, exon, junction, boundary, or intron should have already been captured by mapping directly to these elements
#- This analysis should capture reads that hit entirely within intergenic regions between known EnsEMBL genes (or on chromosome arms)
#- Since intergenic regions are often quite large, it may be possible to map reads as a pair in some cases - but this is not required

#1.)  Parse blast results files one at a time
#1a.) Go through the specified blast results files (in the specified directory)
#    - For each blast results file go through the hits for a single read
#    - Store the best intergenic hits
#    - If the best and second best hit have the same bit score, consider this a tie and make note of it as an ambiguous hit


use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;
use IO::File; 
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
use utilities::ALEXA_DB qw(:all);

my $intergenic_db = '';
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

GetOptions ('intergenic_db=s'=>\$intergenic_db, 'read_records_infile=s'=>\$read_records_infile, 'blast_results_dir=s'=>\$blast_results_dir,
	    'reads_analyzed_per_block=i'=>\$reads_analyzed_per_block, 'total_reads_analyzed=i'=>\$total_reads_analyzed,
	    'working_dir=s'=>\$working_dir, 'summary_outfile=s'=>\$summary_outfile,
	    'min_bit_score=f'=>\$min_bit_score, 'upper_bit_score=f'=>\$upper_bit_score, 'min_align_length=i'=>\$min_align_length,
	    'logfile=s'=>\$logfile);

#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tThis parses blast results for Solexa paired reads BLASTED against a database of Introns", RESET;
print GREEN, "\n\tThis script takes a read record file and blast results files as input, parses the results and creates a summary file", RESET;
print GREEN, "\n\tIt also updates the master read record file to assign perfect matching reads to this source (if they have not already been assigned)", RESET;
print GREEN, "\n\tSpecify an annotated intergenic database file using: --intergenic_db", RESET;
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

print GREEN, "\n\nExample: parseBlastResults_Intergenic.pl  --intergenic_db=/projects/malachig/sequence_databases/hs_49_36k/ensembl_intergenics_hs_49_36k/intergenics_annotated.txt  --read_records_infile=/projects/malachig/solexa/read_records/HS04391/20821AAXX_Lane1.txt.gz  --blast_results_dir=/projects/malachig/solexa/blast_results/HS04391/20821AAXX_Lane1/ensembl_intergenics_v49/  --reads_analyzed_per_block=125000  --total_reads_analyzed=5738676  --working_dir=/projects/malachig/solexa/blast_results/HS04391/20821AAXX_Lane1/temp/  --summary_outfile=/projects/malachig/solexa/read_records/HS04391/Intergenics_v49/20821AAXX_Lane1_Intergenics_v49.txt  --min_bit_score=48.1  --min_align_length=32  --upper_bit_score=60.0  --logfile=/projects/malachig/solexa/logs/HS04391/20821AAXX_Lane1/parseBlastResults_Intergenic_LOG.txt\n\n", RESET;

unless ($intergenic_db && $read_records_infile && $blast_results_dir && $reads_analyzed_per_block && $total_reads_analyzed && $working_dir && $summary_outfile && $min_bit_score && $upper_bit_score && $min_align_length && $logfile){
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  exit();
}

#Check dirs before getting started
$working_dir = &checkDir('-dir'=>$working_dir, '-clear'=>"no");
$blast_results_dir = &checkDir('-dir'=>$blast_results_dir, '-clear'=>"no");

unless (-e $intergenic_db && ($intergenic_db =~ /\.gz/)){
  print RED, "\nAnnotation database file: $intergenic_db does not exist or is not compressed!!\n\n", RESET;
  exit();
}
unless ($read_records_infile =~ /(.*)\.gz$/){
  print RED, "\nFormat of read records infile name not understood: $read_records_infile - not compressed?\n\n", RESET;
  exit();
}

my $log_fh = IO::File->new(">$logfile") || die "\nCould not open log file: $logfile\n\n";

#1.) Create a map of intergenic_IDs to chromosome coordinates by parsing the intergenic annotation database file.
#    - This hash will be needed to calculate the chromosome coordinates for hits to the intergenic sequences
my %intergenics;
open(INTERGENIC_DB, "zcat $intergenic_db |") || die "\nCould not open intergenic db file: $intergenic_db\n\n";
my $intergenic_count = 0;
my $header = 1;
my %columns;

while(<INTERGENIC_DB>){

  chomp($_);
  my $line = $_;
  my @line = split("\t", $line);

  if ($header == 1){
    my $column_count = 0;
    foreach my $column (@line){
      $columns{$column}{column_pos} = $column_count;
      $column_count++;
    }
    $header = 0;
    next();
  }
  $intergenic_count++;
  my $chr = $line[$columns{'Chromosome'}{column_pos}];
  my $intergenic_id = $line[$columns{'Intergenic_ID'}{column_pos}];
  my $strand = $line[$columns{'Strand'}{column_pos}];
  my $start = 1;
  my $end = $line[$columns{'Base_Count'}{column_pos}];
  my $start_chr = $line[$columns{'Unit1_start_chr'}{column_pos}];
  my $end_chr = $line[$columns{'Unit1_end_chr'}{column_pos}];

  $intergenics{$intergenic_id}{count} = $intergenic_count;
  $intergenics{$intergenic_id}{chromosome} = $chr;
  $intergenics{$intergenic_id}{start} = $start;
  $intergenics{$intergenic_id}{end} = $end;
  $intergenics{$intergenic_id}{start_chr} = $start_chr;
  $intergenics{$intergenic_id}{end_chr} = $end_chr;
  $intergenics{$intergenic_id}{size} = ($end-$start)+1;
  $intergenics{$intergenic_id}{strand} = $strand;
}
close(INTERGENIC_DB);
my $message = &memoryUsage();
print YELLOW, "\n\n$message", RESET;
print $log_fh "\n\n$message";


#2.)  Parse blast results files one at a time
#3.)  Print the best hit to the 'top_hits' file

my $library_name;
if ($read_records_infile =~ /\/(\w+_Lane\d+)\.txt/){
  $library_name = $1;
}else{
  print RED, "\nCould not determine library name from read records infile: $read_records_infile\n\n", RESET;
  exit();
}

my $seq_type = "Intergenic";

#Open global Berkley DB files for reading or writing
my %top_hits_db;
my $top_hits_db_file = "$working_dir"."$library_name"."_TopHits.btree";
system ("rm -f $top_hits_db_file");
tie(%top_hits_db, 'BerkeleyDB::Btree', -Cachesize => 256000000, -Filename=> $top_hits_db_file, -Flags => DB_CREATE) or die "can't open file $top_hits_db_file: $! $BerkeleyDB::Error\n";

my $top_hits_stored = &parseBlastFiles('-seq_type'=>$seq_type, '-blast_results_dir'=>$blast_results_dir, '-top_hits_db'=>\%top_hits_db, '-min_bit_score'=>$min_bit_score, 
                                       '-reads_analyzed_per_block'=>$reads_analyzed_per_block, '-total_reads_analyzed'=>$total_reads_analyzed,
                                       '-log_file_handle'=>$log_fh,  '-seq_db'=>\%intergenics);

$message = &memoryUsage();
print YELLOW, "\n\n$message", RESET;
print $log_fh "\n\n$message";

#Cleanly close the berkeley DB files
untie(%top_hits_db);


#Retie top hits file as readonly
tie(%top_hits_db, 'BerkeleyDB::Btree', -Cachesize => 256000000, -Filename=> $top_hits_db_file, -Flags => DB_RDONLY ) or die "can't open file $top_hits_db_file: $! $BerkeleyDB::Error\n";


#4.) Create a summary of top hits and associated info with both reads of a pair on a single line
#    - Count entries in the read records file.  Store all read IDs as an array.
#    - Note that not all read records will have a top hits record.  There may not have been any blast hits to either read of the pair
#    - Parse a BLOCK of full read records.  Create a hash and key on Read_ID
#    - Using this hash and the berkeley DB of top_hits, created a joined read records summary file
#    - Also update the mapping status of the master read record file for the lane

#    - Columns that will be in the summary file:
#    - Read_ID, DistanceBetweenReads, R1_ID, R1_HitType, R1_Intron_ID, R1_AlignmentLength, R1_PercentIdentity, R1_BitScore, R1_Chr_start, R1_Chr_end, R2_ID, R2_HitType, R2_Intron_ID, R2_AlignmentLength, R2_PercentIdentity, R2_BitScore, R2_Chr_start, R2_Chr_end
my $hit_u_count = 0;
my $hit_a_count = 0;

&joinRecords('-subjects'=>\%intergenics,
             '-seq_type'=>$seq_type, '-upper_bit_score'=>$upper_bit_score, '-min_align_length'=>$min_align_length, '-top_hits_db'=>\%top_hits_db, '-hit_u_count'=>\$hit_u_count, '-hit_a_count'=>\$hit_a_count,
             '-read_records_file'=>$read_records_infile, '-summary_outfile'=>$summary_outfile, '-library_name'=>$library_name, '-working_dir'=>$working_dir,
             '-log_file_handle'=>$log_fh);

$message = &memoryUsage();
print YELLOW, "\n\n$message", RESET;
print $log_fh "\n\n$message";

print "\n\n";
print $log_fh "\n\n";

untie(%top_hits_db);
print BLUE, "\nDeleting top hits db file: $top_hits_db_file\n\n", RESET;
print $log_fh "\nDeleting top hits db file: $top_hits_db_file\n\n";

my $cmd = "rm -f $top_hits_db_file";
system($cmd);
$message = &memoryUsage();
print YELLOW, "\n\n$message", RESET;
print $log_fh "\n\n$message";

print $log_fh "\n\nSCRIPT COMPLETE\n\n";
close($log_fh);
print "\n\nSCRIPT COMPLETE\n\n";

exit();

