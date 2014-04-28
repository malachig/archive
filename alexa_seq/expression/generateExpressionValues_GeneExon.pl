#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#The purpose of this script is to generate gene and exon level expression values for all genes of a chromosome
#Exons are summarized at the level of 'Exon Regions' (i.e. in the case of overlaping exons with varying boundaries, each section is summarize seperately)
#The input is a tabular file summarzing the top hits of WTSS read pairs to a database of all EnsEMBL transcripts
#The user also inputs a previously generated Exon Region database file, consisting of annotated exon regions
#The outputs are gene and exon region expression value files (one line per gene or exon region) and UCSC wiggle tracks summarizing base level mapability and expression

#Each GENE will be assigned values as follows:

#READ COUNTS
#Raw_Read_Count - # number of unambiguous reads meeting a bit score cutoff
#Norm_Read_Count - Read count normalized to a user specified library size (use the # of quality, unambiguously mapped reads for the smaller library)

#Each GENE/EXON will be assigned values as follows:

#AVERAGE BASE COVERAGE
#Average_Base_Coverage - Total base coverage divided by number of bases in a gene/exon
#Norm1_Average_Base_Coverge - Normalized to user specified library size (use the # of quality, unambiguously mapped reads for the smaller library)
#Norm2_Average_Base_Coverage - Normalized to account for mappability of bases in a gene/exon


#STEPS
#1.) For the user specified chromosome, get all reads out of the input mapped reads file
#    - Store a Berkley DB file in a working dir (or as a flat file)

#2.) Get basic gene info from ALEXA (genes, exon content coords, etc.)

#3.) Import mappability scores
#    - From a Berkley DB

#4.) Build a base-by-base coverage object for each GENE (store in memory or as a Berkley DB)
#    - At the same time, Build a base-by-base coverage object for the entire CHROMOSOME

#    - GENE coverage object
#    - Use a gene_position key as follows: key{'GeneID_ChromosomePosition'} value{read_base_coverage}
#    - Initialize to 0 for all valid theoretical exonic bases

#    - CHROMOSOME coverage object
#    - Use a position key as follows: key{'ChromosomePosition'} value{read_base_coverage}
#    - Do NOT initialize - only bases actually covered by reads will be represented

#5.) Correct base count values at the GENE and CHROMOSOME level
#    - correct base counts to account for library size (normalize to a user specified library size - use the number of mapped reads in the smaller library)
#    - correct base counts to account for mapability

#6.) Import the exon regions database file

#7.) Summarize gene and exon level expression values.  
#    - Also create files that report percent coverage for each gene (percent of all sequenced bases) versus percent position in the gene (5'end -> 1%, 2%, ... 100% -> 3'end)

#8.) Generate a UCSC wig file for the entire chromosome
#    - Use the library size and mapability corrected base coverage values


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
use utilities::ALEXA_DB qw(:all);
use utilities::utility qw(:all);
use utilities::mapping qw(:all);

my $database = '';
my $server = '';
my $user = '';
my $password = '';
my $library = '';
my $library_name = '';
my $chr_filter = '';
my $read_record_dir = '';
my $mapped_reads_dir = '';
my $annotation_dir = '';
my $min_bit_score = '';
my $min_seq_coverage = '';
my $working_dir = '';
my $results_dir = '';
my $ucsc_dir = '';
my $ucsc_build = '';
my $web_path = '';
my $color_set = '';
my $cutoffs_file = '';
my $log_file = '';
my $berk = '';

GetOptions ('database=s'=>\$database,'server=s'=>\$server, 'user=s'=>\$user, 'password=s'=>\$password,
	    'library=s'=>\$library, 'library_name=s'=>\$library_name, 'ucsc_build=s'=>\$ucsc_build,
	    'chr_filter=s'=>\$chr_filter, 'read_record_dir=s'=>\$read_record_dir, 'mapped_reads_dir=s'=>\$mapped_reads_dir, 'annotation_dir=s'=>\$annotation_dir,
            'min_bit_score=f'=>\$min_bit_score, 'min_seq_coverage=f'=>\$min_seq_coverage, 'working_dir=s'=>\$working_dir, 'results_dir=s'=>\$results_dir, 
            'ucsc_dir=s'=>\$ucsc_dir, 'web_path=s'=>\$web_path, 'color_set=i'=>\$color_set, 'cutoffs_file=s'=>\$cutoffs_file, 'log_file=s'=>\$log_file, 'berk=i'=>\$berk);

#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the ALEXA database and server to query using: --database and --server", RESET;
print GREEN, "\n\tSpecify the ALEXA user and password for access using: --user and --password", RESET;
print GREEN, "\n\tSpecify the library ID for read data using: --library", RESET;
print GREEN, "\n\tSpecify the library name for read data using: --library_name", RESET;
print GREEN, "\n\tSpecify the chromosome to be processed using: --chr_filter", RESET;
print GREEN, "\n\tSpecify a directory containing files of reads mapped to EnsEMBL genes using: --mapped_reads_dir", RESET;
print GREEN, "\n\tSpecify a directory containing read record files for this library using:  --read_record_dir", RESET;
print GREEN, "\n\t\tThis directory should contain all of the lanes for a single library only - and the files must be compressed", RESET;
print GREEN, "\n\tSpecify a directory containing gene and exon region annotations using: --annotation_dir", RESET;
print GREEN, "\n\tSpecify the minimum bit score to a transcript for each read to be considered for the summary using: --min_bit_score", RESET;
print GREEN, "\n\tSpecify the minimum percent sequence coverage for an observed sequence to be considered expressed and written to the UCSC custom track using: --min_seq_coverage", RESET;
print GREEN, "\n\t\tThis is the percentage of bases of an Exon Region covered to 1X depth or greater", RESET;
print GREEN, "\n\tSpecify the working directory to write binary tree files to using: --working_dir", RESET;
print GREEN, "\n\tSpecify the a directory to write final GENE/EXON results files to using: --results_dir", RESET;
print GREEN, "\n\tSpecify the target UCSC directory for custom UCSC track files using: --ucsc_dir", RESET;
print GREEN, "\n\tSpecify the UCSC build to be used for links to UCSC tracks using: --ucsc_build (e.g. hg18)", RESET;
print GREEN, "\n\tSpecify the html web path to this directory using: --web_path", RESET;
print GREEN, "\n\tSpecify which hard coded color set to use using: --color_set (1 for LibraryA, 2 for LibraryB)", RESET;
print GREEN, "\n\tSpecify the path to a file containing expression cutoffs values using:  --cutoffs_file", RESET;
print GREEN, "\n\t\tIf these have not been calculated yet, use: --cutoffs=0", RESET;
print GREEN, "\n\tSpecify a log file using: --log_file", RESET;
print GREEN, "\n\tFor large memory jobs you can try: --berk=1 (will use a Berkeley DB for large memory jobs)", RESET;

print GREEN, "\n\nExample: generateExpressionValues_GeneExon.pl  --database=ALEXA_hs_53_36o  --server=jango.bcgsc.ca  --user=viewer  --password=viewer  --library=HS04391  --library_name=MIP101  --chr_filter='3:10:118167002-127874001'  --working_dir=/projects/malachig/solexa/read_records/HS04391/ENST_v49/Summary/temp/  --min_bit_score=60.0  --min_seq_coverage=75.0  --mapped_reads_dir=/projects/malachig/solexa/read_records/HS04391/ENST_v53/  --read_record_dir=/projects/malachig/solexa/read_records/HS04391/  --annotation_dir=/projects/malachig/sequence_databases/hs_53_36o/   --results_dir=/projects/malachig/solexa/read_records/HS04391/ENST_v53/Summary/gene_exon_results/  --ucsc_dir=/home/malachig/www/public/htdocs/solexa/HS04391/  --web_path=http://www.bcgsc.ca/people/malachig/htdocs/solexa/HS04391/  --color_set=1  --cutoffs_file=/projects/malachig/solexa/figures_and_stats/HS04391/Expression_v53/HS04391_NORM1_average_coverage_cutoffs.txt  --log_file=/projects/malachig/solexa/logs/HS04391/generateExpressionValues/GeneExon/generateExpressionValues_GeneExon_chr3_10.txt\n\n", RESET;

unless ($database && $server && $user && $password && $library && $library_name && $chr_filter && $annotation_dir && $mapped_reads_dir && $read_record_dir && $min_bit_score && $min_seq_coverage && $working_dir && $results_dir && $ucsc_dir && $ucsc_build && $web_path && $color_set && ($cutoffs_file || $cutoffs_file eq '0') && $log_file){
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  exit();
}
if ($berk){
  $berk = 1;
}else{
  $berk = 0;
}

#Scale all expression values to a constant value (values will be scaled to what they would be expected to be if the library contained 10 billion mapped bases)
#If a library has less than this, the normalized expression values will be increased and vice versa
#The actual number of mapped bases is determined by parsing the 'read record' files and counting ALL bases that have been unambiguously assigned to one of the following classes:
#'ENST_U', 'INTRON_U', 'INTERGENIC_U', 'NOVEL_JUNCTION_U', 'NOVEL_BOUNDARY_U' 
my $library_normalization_value = 10000000000;
my $library_size;
 
#Set priority values to determine ordering of custom tracks
my $priority_annotation = -1;
my $priority_expression = 30;
my $priority_wiggle = 60;

#Define some color profiles
my %colors;
$colors{'1'}{0} = "153,0,0";      #Brick Red
$colors{'1'}{1} = "153,51,0";     #Red-Brown
$colors{'1'}{2} = "0,51,102";     #Dark Blue
$colors{'2'}{0} = "153,0,0";      #Brick Red
$colors{'2'}{1} = "153,51,0";     #Red-Brown
$colors{'2'}{2} = "0,51,204";     #Medium Blue

unless ($colors{$color_set}){
  print RED, "\nColor set specified by user is not understood", RESET;
  print Dumper %colors;
  exit();
}

my $file_name_prefix;
my $region_number;
my $start_filter;
my $end_filter;
if ($chr_filter =~ /(.*)\:(\d+)\:(\d+)\-(\d+)/){
  $chr_filter = $1;
  $region_number = $2;
  $start_filter = $3;
  $end_filter = $4;
  $file_name_prefix = "$chr_filter"."_"."$region_number";
  unless ($end_filter > $start_filter){
    print RED, "\nStart of range must be smaller than end ($chr_filter)\n\n", RESET;
    exit();
  }
}else{
  print RED, "\nFormat of chr_filter not understood: $chr_filter (should be of the form:  Y:1:1-9999001)\n\n", RESET;
  exit();
}

if ($chr_filter eq "MT"){$chr_filter = "M"};

#Check working dirs before proceeding
#Multiple scripts may be writing to the directory, do not clear
#All files will be names chrN______
$mapped_reads_dir = &checkDir('-dir'=>$mapped_reads_dir, '-clear'=>"no");
$read_record_dir = &checkDir('-dir'=>$read_record_dir, '-clear'=>"no");
$working_dir = &checkDir('-dir'=>$working_dir, '-clear'=>"no");
$ucsc_dir = &checkDir('-dir'=>$ucsc_dir, '-clear'=>"no");
$results_dir = &checkDir('-dir'=>$results_dir, '-clear'=>"no");
$annotation_dir = &checkDir('-dir'=>$annotation_dir, '-clear'=>"no");

open (LOG, ">$log_file") || die "\nCould not open log file: $log_file\n\n";

#Global variable storing letters for labeling of Acceptors and Donors
#my @letters = ('aa' .. 'zz');
my @letters = ('a' .. 'z');


#Open global Berkley DB files for reading or writing
my %chr_reads;
my %ensembl_transcript_tracks;
my %gene_coverage;
my %gene_coverage_norm1;
my %chr_coverage;
my %chr_coverage_norm1;
my %ensembl_exon_annotation_tracks;
my %ensembl_exon_region_tracks;
my %expressed_exon_region_tracks;

my $rm_cmd = "rm -f $working_dir"."*chr"."$file_name_prefix"."_"."*";
if ($berk == 1){
  #Delete pre-existing berkley DB files for the current chromosome
  print YELLOW, "\nDeleting Berkley DB files and creating new\n\t($rm_cmd)\n\n", RESET;
  print LOG "\nDeleting Berkley DB files and creating new\n\t($rm_cmd)\n\n";
  system($rm_cmd);
 
  my $chr_reads_file = "$working_dir"."chr"."$file_name_prefix"."_MappedReads.btree";
  tie(%chr_reads, 'BerkeleyDB::Btree', -Cachesize => 256000000, -Filename=> $chr_reads_file , -Flags => DB_CREATE) or die "can't open file $chr_reads_file: $! $BerkeleyDB::Error\n";

  my $ensembl_transcript_tracks_file = "$working_dir"."chr"."$file_name_prefix"."_EnsEMBL_Transcripts.btree";
  tie(%ensembl_transcript_tracks, 'BerkeleyDB::Btree', -Cachesize => 256000000, -Filename=> $ensembl_transcript_tracks_file, -Flags => DB_CREATE) or die "can't open file $ensembl_transcript_tracks_file: $! $BerkeleyDB::Error\n";

  my $gene_coverage_file = "$working_dir"."chr"."$file_name_prefix"."_GeneCoverage.btree";
  tie(%gene_coverage, 'BerkeleyDB::Btree', -Cachesize => 256000000, -Filename=> $gene_coverage_file, -Flags => DB_CREATE) or die "can't open file $gene_coverage_file: $! $BerkeleyDB::Error\n";

  my $gene_coverage_norm1_file = "$working_dir"."chr"."$file_name_prefix"."_GeneCoverage_norm1.btree";
  tie(%gene_coverage_norm1, 'BerkeleyDB::Btree', -Cachesize => 256000000, -Filename=> $gene_coverage_norm1_file, -Flags => DB_CREATE) or die "can't open file $gene_coverage_norm1_file: $! $BerkeleyDB::Error\n";

  my $chr_coverage_file = "$working_dir"."chr"."$file_name_prefix"."_ChrCoverage.btree";
  tie(%chr_coverage, 'BerkeleyDB::Btree', -Cachesize => 256000000, -Filename=> $chr_coverage_file, -Flags => DB_CREATE) or die "can't open file $chr_coverage_file: $! $BerkeleyDB::Error\n";

  my $chr_coverage_norm1_file = "$working_dir"."chr"."$file_name_prefix"."_ChrCoverage_norm1.btree";
  tie(%chr_coverage_norm1, 'BerkeleyDB::Btree', -Cachesize => 256000000, -Filename=> $chr_coverage_norm1_file, -Flags => DB_CREATE) or die "can't open file $chr_coverage_norm1_file: $! $BerkeleyDB::Error\n";

  my $chr_exon_annotations_file = "$working_dir"."chr"."$file_name_prefix"."_ExonAnnotations.btree";
  tie(%ensembl_exon_annotation_tracks, 'BerkeleyDB::Btree', -Cachesize => 256000000, -Filename=> $chr_exon_annotations_file, -Flags => DB_CREATE) or die "can't open file $chr_exon_annotations_file: $! $BerkeleyDB::Error\n";

  my $chr_exon_region_file = "$working_dir"."chr"."$file_name_prefix"."_ExonRegions.btree";
  tie(%ensembl_exon_region_tracks, 'BerkeleyDB::Btree', -Cachesize => 256000000, -Filename=> $chr_exon_region_file, -Flags => DB_CREATE) or die "can't open file $chr_exon_annotations_file: $! $BerkeleyDB::Error\n";

  my $expressed_exon_region_file = "$working_dir"."chr"."$file_name_prefix"."_ExpressedExonRegions.btree";
  tie(%expressed_exon_region_tracks, 'BerkeleyDB::Btree', -Cachesize => 256000000, -Filename=> $expressed_exon_region_file, -Flags => DB_CREATE) or die "can't open file $expressed_exon_region_file: $! $BerkeleyDB::Error\n";
}


#0.) If specified, get the gene-by-gene expression cutoffs.
#    - These will be used to decide whether a particular sequence is expressed or not
my $gene_cutoffs_ref;
if ($cutoffs_file && -e $cutoffs_file){
  $gene_cutoffs_ref = &importExpressionCutoffs ('-cutoffs_file'=>$cutoffs_file);
}else{
  $cutoffs_file = 0;
  print YELLOW, "\nCutoffs file not specified - or not found, expression will be evaluated by percent coverage only\n\n", RESET;
}


#1-A.) For the user specified chromosome, get all reads out of the input mapped reads file
#    - Store a Berkley DB file in a working dir (or as a flat file)
my $grand_total_mapped_reads = 0;  #Count all the quality unambiguous mapped reads from ANY chromosome
&parseMappedReads('-indir'=>$mapped_reads_dir, '-chromosome'=>$chr_filter, '-range'=>"$start_filter-$end_filter", '-min_bit_score'=>$min_bit_score, '-working_dir'=>$working_dir);

#1-B.) Go through the read records files and remove reads from those stored if they have not been assigned to the correct class
my $read_class = "ENST_U";
$library_size = &importLibrarySize ('-read_records_dir'=>$read_record_dir);


#2.) Get basic gene info from ALEXA (genes, exon content coords, etc.)
my $genes_ref;
my $gene_transcripts_ref;
my $gene_exon_content_ref;
&getBasicGeneInfo('-chromosome'=>$chr_filter, '-range'=>"$start_filter-$end_filter");

#Import the exon regions database file 
my $exon_region_file = "$annotation_dir"."exonRegions/exonRegions_annotated.txt.gz";
my $gene_file = "$annotation_dir"."genes/genes_annotated.txt.gz";
unless (-e $exon_region_file && -e $gene_file){
  print RED, "\nCould not find gene or exon region annotation files:\n\t$gene_file\n\t$exon_region_file\n\n", RESET;
  exit();
}

#Exon regions
my $er_header;
my $exon_regions_ref = &importExonRegions('-infile'=>$exon_region_file, '-chromosome'=>$chr_filter, '-range'=>"$start_filter-$end_filter");

#For genes, simply add data to the existing genes_ref
my $gene_header;
&importGenes('-infile'=>$gene_file, '-chromosome'=>$chr_filter, '-range'=>"$start_filter-$end_filter");


#3.) Import mappability scores

#4.) Build a base-by-base coverage object for each GENE (store in memory or as a Berkley DB)
#    - At the same time, Build a base-by-base coverage object for the entire CHROMOSOME

#    - GENE coverage object
#    - Use a gene_position key as follows: key{'GeneID_ChromosomePosition'} value{read_base_coverage}
#    - Initialize to 0 for all valid theoretical exonic bases

#    - CHROMOSOME coverage object
#    - Use a position key as follows: key{'ChromosomePosition'} value{read_base_coverage}
#    - Do NOT initialize - only bases actually covered by reads will be represented
&buildCoverageObjects('-chromosome'=>$chr_filter, '-range'=>"$start_filter-$end_filter", '-working_dir'=>$working_dir);


#5.) Correct base count values at the GENE and CHROMOSOME level
#    - correct base counts to account for library size (normalize to a user specified library size - use the number of mapped reads in the smaller library)
#    - correct base counts to account for mappability
&correctBaseCounts('-normalization_value'=>$library_normalization_value, '-library_size'=>$library_size);


#7.) Summarize gene and exon level expression values
&printGeneExonSummary('-results_dir'=>$results_dir);


#8.) Generate a UCSC wig file for the entire chromosome
#    - Use the library size corrected base coverage values
&printUCSCTrack('-ucsc_dir'=>$ucsc_dir, '-chr'=>$chr_filter, '-library'=>$library);

#If neccessary, clean berkeley DB files
if ($berk == 1){
  #Cleanly untie the berkeley dbs
  untie (%ensembl_transcript_tracks);
  untie (%chr_coverage);
  untie (%chr_coverage_norm1);
  untie (%gene_coverage);
  untie (%gene_coverage_norm1);
  untie (%ensembl_exon_annotation_tracks);
  untie (%ensembl_exon_region_tracks);
  untie (%expressed_exon_region_tracks);
  untie (%chr_reads);
  #Delete temp berkley DB files created by this script in the working dir (for the target chromosome)
  print YELLOW, "\n\nScript complete, deleting Berkley DB files\n\t($rm_cmd)\n\n", RESET;
  print LOG "\n\nScript complete, deleting Berkley DB files\n\t($rm_cmd)\n\n";
  system($rm_cmd);
}

#Summarize the total memory usage at close (since Perl doesnt usually release memory ... this should be the max used by the script):
my $message = &memoryUsage();
print YELLOW, "\n\n$message", RESET; 
print LOG "\n\n$message"; 

print BLUE, "\n\nSCRIPT COMPLETE\n\n", RESET;
print LOG "\n\nSCRIPT COMPLETE\n\n";

close(LOG);

exit();


################################################################################################################################################
#For the user specified chromosome, get all reads out of the input mapped reads file (above a certain quality) and store seperately
################################################################################################################################################
sub parseMappedReads{
  my %args = @_;
  my $indir = $args{'-indir'};
  my $target_chr = $args{'-chromosome'};
  my $range = $args{'-range'};
  my $min_bit_score = $args{'-min_bit_score'};
  my $working_dir = $args{'-working_dir'};

  $target_chr = "chr"."$target_chr";

  my $start_filter;
  my $end_filter;
  if ($range =~ /(\d+)\-(\d+)/){
    $start_filter = $1;
    $end_filter = $2;
  }else{
    print RED, "\nProblem with range - parseMappedReads()\n\n", RESET;
    exit();
  } 

  print BLUE, "\n1.) Begin parsing input files for target chr ($target_chr: $start_filter-$end_filter) reads: $indir", RESET;
  print LOG "\n1.) Begin parsing input files for target chr ($target_chr: $start_filter-$end_filter) reads: $indir";

  my %columns;
  my @required_columns = qw(Read_ID R1_ID R1_HitType R1_GeneID R1_AlignmentLength R1_Chromosome R1_ChrStartCoords R1_ChrEndCoords R2_ID R2_HitType R2_GeneID R2_AlignmentLength R2_Chromosome R2_ChrStartCoords R2_ChrEndCoords);

  my $passing_reads = 0;

  #Get files from this directory
  print BLUE, "\n\nSearching $indir for mapped read files", RESET;
  print LOG "\n\nSearching $indir for mapped read files";

  my %mapped_read_files;
  opendir(DIRHANDLE, "$indir") || die "\nCannot open directory: $indir\n\n";
  my @test_files = readdir(DIRHANDLE);
  my $file_count = 0;

  foreach my $test_file (sort @test_files){
    my $file_path = "$indir"."$test_file";

    #Skip directories within the specified directory
    if (-e $file_path && -d $file_path){
      print YELLOW, "\n\t$file_path  is a directory - skipping", RESET;
      print LOG "\n\t$file_path  is a directory - skipping";
      next();
    }

    unless($test_file =~ /\_Lane\d+\_/){
      print YELLOW, "\n\t$test_file does not have the expected name format - skipping", RESET;
      print LOG "\n\t$test_file does not have the expected name format - skipping";
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

    $mapped_read_files{$file_count}{path} = $file_path;
  }

  my $num_files = keys %mapped_read_files;
  print BLUE, "\n\nBegin parsing $num_files mapped read files", RESET;
  print LOG "\n\nBegin parsing $num_files mapped read files";

  foreach my $file_count (sort {$mapped_read_files{$a}{path} cmp $mapped_read_files{$b}{path}} keys %mapped_read_files){
    my $infile = $mapped_read_files{$file_count}{path};
    print BLUE, "\n\nProcessing file: $infile", RESET;
    print LOG "\n\nProcessing file: $infile";

    my $line_counter = 0;
    my $block_count = 0;
    my $block_size = 1000000;
    my $header_line = 1;
    my $header;

    open (READ, "zcat $infile |") || die "\nCould not open read records summary infile: $infile\n\n";
    while(<READ>){
      $line_counter++;

      if ($line_counter == $block_size){
        $block_count++;
        
        my $message = &memoryUsage();
        $| = 1; print BLUE, "\n\tBlock $block_count (of size $block_size).  Found $passing_reads unambigous reads (cumulative) matching the target chromosome $target_chr:$range (bit score >= $min_bit_score) - $message", RESET; $| = 0;
        print LOG "\n\tBlock $block_count (of size $block_size).  Found $passing_reads unambigous reads (cumulative) matching the target chromosome $target_chr:$range (bit score >= $min_bit_score) - $message";

        $line_counter = 0;
      }

      chomp($_);
      my $current_line = $_;
      my @line = split("\t", $_);

      #Parse the column names and positions.  Check against a hard coded list of required columns before proceeding
      if ($header_line == 1){
        $header = $current_line;
        my $col_count = 0;
        foreach my $column (@line){
	  $columns{$column}{position} = $col_count;
	  $col_count++;
        }

        foreach my $req_column (@required_columns){
	  unless ($columns{$req_column}){
	    print RED, "\nRequired column: $req_column was not found in the read record file!\n\n", RESET;
	    exit();
	  }
        }
        $header_line = 0;
        next();
      }

      my $r1_id = $line[$columns{R1_ID}{position}];
      my $r2_id = $line[$columns{R2_ID}{position}];
      my $r1_gene_id = $line[$columns{R1_GeneID}{position}];
      my $r2_gene_id = $line[$columns{R2_GeneID}{position}];
      my $r1_chromosome = $line[$columns{R1_Chromosome}{position}];
      my $r2_chromosome = $line[$columns{R2_Chromosome}{position}];
      my $r1_hit_type = $line[$columns{R1_HitType}{position}];
      my $r2_hit_type = $line[$columns{R2_HitType}{position}];
      my $r1_bit_score = $line[$columns{R1_BitScore}{position}];
      my $r2_bit_score = $line[$columns{R2_BitScore}{position}];
      my $r1_chr_start_coords = $line[$columns{R1_ChrStartCoords}{position}];
      my $r2_chr_start_coords = $line[$columns{R2_ChrStartCoords}{position}];
      my $r1_chr_end_coords = $line[$columns{R1_ChrEndCoords}{position}];
      my $r2_chr_end_coords = $line[$columns{R2_ChrEndCoords}{position}];

      #Fix chromosome formats
      if ($r1_chromosome eq "chrMT"){$r1_chromosome = "chrM";}
      if ($r2_chromosome eq "chrMT"){$r2_chromosome = "chrM";}

      #Count quality, unambiguous mapped reads to ANY chromosome
      if (($r1_hit_type eq "Top_Hit") && ($r1_bit_score >= $min_bit_score)){
        $grand_total_mapped_reads++;
      }
      if (($r2_hit_type eq "Top_Hit") && ($r2_bit_score >= $min_bit_score)){
        $grand_total_mapped_reads++;
      }

      #Before proceeding, make sure at least one of these reads corresponds to the current chromosome
      unless (($r1_chromosome eq "$target_chr") || ($r2_chromosome eq "$target_chr")){
        next();
      }

      if (($r1_hit_type eq "Top_Hit") && ($r1_bit_score >= $min_bit_score) && ($r1_chromosome eq "$target_chr")){
        my @chr_starts = split(" ", $r1_chr_start_coords);
        my @chr_ends = split(" ", $r1_chr_end_coords);
        my @coords = (@chr_starts, @chr_ends);
        my @coords_sort = sort {$a <=> $b} @coords;
        my $lower = $coords_sort[0];
        my $upper = $coords_sort[scalar(@coords)-1];
        if ($lower >= $start_filter && $lower <= $end_filter && $upper >= $start_filter && $upper <= $end_filter){
          $passing_reads++;
          my $string = "$r1_gene_id\t$r1_chr_start_coords\t$r1_chr_end_coords";
          #$chr_reads{$r1_id} = $string;
          if ($chr_reads{$string}){
            $chr_reads{$string}++;
          }else{
            $chr_reads{$string}=1;
          }
        }
      }

      #Deal with Read 2
      if (($r2_hit_type eq "Top_Hit") && ($r2_bit_score >= $min_bit_score) && ($r2_chromosome eq "$target_chr")){
        my @chr_starts = split(" ", $r2_chr_start_coords);
        my @chr_ends = split(" ", $r2_chr_end_coords);
        my @coords = (@chr_starts, @chr_ends);
        my @coords_sort = sort {$a <=> $b} @coords;
        my $lower = $coords_sort[0];
        my $upper = $coords_sort[scalar(@coords)-1];
        if ($lower >= $start_filter && $lower <= $end_filter && $upper >= $start_filter && $upper <= $end_filter){
          $passing_reads++;
          my $string = "$r2_gene_id\t$r2_chr_start_coords\t$r2_chr_end_coords";
          #$chr_reads{$r2_id} = $string;
          if ($chr_reads{$string}){
            $chr_reads{$string}++;
          }else{
            $chr_reads{$string}=1;
          }
        }
      }
    }
  }

  $| = 1; print BLUE, "\n\n\tFound a total of $passing_reads unambigous reads matching the target chromosome (bit score >= $min_bit_score)", RESET; $| = 0;
  $| = 1; print BLUE, "\n\tFound a grand total of $grand_total_mapped_reads unambigous reads matching ANY chromosome (bit score >= $min_bit_score)", RESET; $| = 0;
  print LOG "\n\n\tFound a total of $passing_reads unambigous reads matching the target chromosome (bit score >= $min_bit_score)";
  print LOG "\n\tFound a grand total of $grand_total_mapped_reads unambigous reads matching ANY chromosome (bit score >= $min_bit_score)";


  return();
}


############################################################################################################################################
#Get basic info for all genes from the user specified ALEXA database                                                                       #
############################################################################################################################################
sub getBasicGeneInfo{
  my %args = @_;
  my $target_chr = $args{'-chromosome'};
  my $range = $args{'-range'};

  print BLUE, "\n\n2-a.) Getting basic gene data", RESET;
  print LOG "\n\n2-a.) Getting basic gene data";

  #Establish connection with the Alternative Splicing Expression database
  my $alexa_dbh = &connectDB('-database'=>$database, '-server'=>$server, '-user'=>$user, '-password'=>$password);

  my @gene_ids = @{&getAllGenes ('-dbh'=>$alexa_dbh, '-gene_type'=>'All', '-chromosome'=>$target_chr, '-range'=>$range)};

  #Get the gene info for all genes for which reads were found on the current chromosome
  $genes_ref = &getGeneInfo ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids, '-sequence'=>"no");

  #Get the transcript info for all transcripts of these genes
  $| = 1; print BLUE, "\n2-b.) Getting transcript data as well as exons for each transcript", RESET; $| = 0;
  print LOG "\n2-b.) Getting transcript data as well as exons for each transcript";
  $gene_transcripts_ref = &getTranscripts('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids, '-sequence'=>"no");

  #Get chromosome coordinates for all EnsEMBL transcripts and Build UCSC tracks for all ensembl transcripts
  $| = 1; print BLUE, "\n\n2-c.) Calculating chromosome coordinates for the EXONS of each gene (only for genes of the current chromosome)", RESET; $| = 0;
  print LOG "\n\n2-c.) Calculating chromosome coordinates for the EXONS of each gene (only for genes of the current chromosome)";

  foreach my $gene_id (keys %{$gene_transcripts_ref}){

    my $chromosome = $genes_ref->{$gene_id}->{chromosome};
    my $chr_strand = $genes_ref->{$gene_id}->{chr_strand};
    my $chr_start = $genes_ref->{$gene_id}->{chr_start};
    my $chr_end = $genes_ref->{$gene_id}->{chr_end};
    my $gene_start = $genes_ref->{$gene_id}->{gene_start};
    my $gene_end = $genes_ref->{$gene_id}->{gene_end};

    if ($chromosome eq "MT"){$chromosome = "M";}
    unless ($chromosome eq $target_chr){
      next();
    }

    my $transcripts_ref = $gene_transcripts_ref->{$gene_id}->{transcripts};

    foreach my $trans_id (keys %{$transcripts_ref}){
      my $exons_ref = $transcripts_ref->{$trans_id}->{exons};

      foreach my $exon_id (keys %{$exons_ref}){

	my $start = $exons_ref->{$exon_id}->{exon_start};
	my $end = $exons_ref->{$exon_id}->{exon_end};

        my $coords_ref = &convertGeneCoordinates ('-gene_object'=>$genes_ref, '-gene_id'=>$gene_id, '-start_pos'=>$start, '-end_pos'=>$end, '-ordered'=>"yes");

        $exons_ref->{$exon_id}->{chr_start} = $coords_ref->{$gene_id}->{chr_start};
	$exons_ref->{$exon_id}->{chr_end} = $coords_ref->{$gene_id}->{chr_end};
	$exons_ref->{$exon_id}->{strand} = $coords_ref->{$gene_id}->{strand};

      }

      #Now use the exon chromosome coordinates to generate custom UCSC track records for this transcript
      my $ensembl_t_id = $transcripts_ref->{$trans_id}->{ensembl_t_id};
      foreach my $exon_id (sort {$exons_ref->{$a}->{chr_start} <=> $exons_ref->{$b}->{chr_start}} keys %{$exons_ref}){
	my $start = $exons_ref->{$exon_id}->{chr_start};
	my $end = $exons_ref->{$exon_id}->{chr_end};
	my $strand = $exons_ref->{$exon_id}->{strand};

	my $ucsc_chromosome = "chr"."$target_chr";
	my $record = "\n$ucsc_chromosome\tEnsEMBL\texon\t$start\t$end\t.\t$strand\t.\t$ensembl_t_id";
	my $record_id = "$ensembl_t_id"."_"."$exon_id";
	$ensembl_transcript_tracks{$record_id} = $record;
	
      }
    }

    #Initialize gene values to be populated later
    $genes_ref->{$gene_id}->{quality_read_count} = 0;
    $genes_ref->{$gene_id}->{cumulative_coverage_raw} = 0;
    $genes_ref->{$gene_id}->{cumulative_coverage_norm1} = 0;
    $genes_ref->{$gene_id}->{bases_covered_1x} = 0;
    $genes_ref->{$gene_id}->{bases_covered_5x} = 0;
    $genes_ref->{$gene_id}->{bases_covered_10x} = 0;
    $genes_ref->{$gene_id}->{bases_covered_50x} = 0;
    $genes_ref->{$gene_id}->{bases_covered_100x} = 0;
    $genes_ref->{$gene_id}->{bases_covered_500x} = 0;
  }

  #Get exon content for genes of the chr
  $| = 1; print BLUE, "\n\n2-d.) Getting EXON CONTENT of each gene", RESET; $| = 0;
  print LOG "\n\n2-d.) Getting EXON CONTENT of each gene";

  $gene_exon_content_ref = &getExonContent ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids);

  #At this time, also get the chromosome coordinates for all exon-content coordinates
  $| = 1; print BLUE, "\n\n2-e.) Calculating chromosome coordinates for the EXON CONTENT of each gene (for only genes of the current chromosome)\n", RESET; $| = 0;
  $| = 1; print BLUE, "\n - Also initializing a COVERAGE hash for the EXON CONTENT of each gene\n", RESET; $| = 0;
  $| = 1; print BLUE, "\n - Also annotating all exons (order, acceptor and donor counts.  e.g. a_E3_a)\n", RESET; $| = 0;

  print LOG "\n\n2-e.) Calculating chromosome coordinates for the EXON CONTENT of each gene (for only genes of the current chromosome)\n";
  print LOG "\n - Also initializing a COVERAGE hash for the EXON CONTENT of each gene\n";
  print LOG "\n - Also annotating all exons (order, acceptor and donor counts.  e.g. a_E3_a)\n";

  my $counter = 0;
  my $ref_exon_counter = 0;

  foreach my $gene_id (@gene_ids){

    my %exon_content_chr;

    my $chromosome = $genes_ref->{$gene_id}->{chromosome};
    my $chr_strand = $genes_ref->{$gene_id}->{chr_strand};
    my $chr_start = $genes_ref->{$gene_id}->{chr_start};
    my $chr_end = $genes_ref->{$gene_id}->{chr_end};
    my $gene_start = $genes_ref->{$gene_id}->{gene_start};
    my $gene_end = $genes_ref->{$gene_id}->{gene_end};

    if ($chromosome eq "MT"){$chromosome = "M";}
    unless ($chromosome eq $target_chr){
      next();
    }

    $counter++;
    if ($counter == 100){
      $counter = 0;
      $| = 1; print BLUE, ".", RESET; $| = 0;
      print LOG ".";
    }

    #Calculate the size of each transcript by adding up the size of its exons
    my $size = 0;
    my $exon_content_ref = $gene_exon_content_ref->{$gene_id}->{exon_content};

    $genes_ref->{$gene_id}->{base_count} = 0;           #Total number of bases covered by exon content of this gene

    #Get the chromosome coordinates for exon content blocks
    foreach my $exon_id (sort {$exon_content_ref->{$a}->{start} <=> $exon_content_ref->{$b}->{start}} keys %{$exon_content_ref}){

      my $start = $exon_content_ref->{$exon_id}->{start};
      my $end = $exon_content_ref->{$exon_id}->{end};

      my $coords_ref = &convertGeneCoordinates ('-gene_object'=>$genes_ref, '-gene_id'=>$gene_id, '-start_pos'=>$start, '-end_pos'=>$end, '-ordered'=>"yes");
      $exon_content_chr{$exon_id}{chr_start} = $coords_ref->{$gene_id}->{chr_start};
      $exon_content_chr{$exon_id}{chr_end} = $coords_ref->{$gene_id}->{chr_end};
      $exon_content_chr{$exon_id}{strand} = $coords_ref->{$gene_id}->{strand};
      my $size = ($coords_ref->{$gene_id}->{chr_end} - $coords_ref->{$gene_id}->{chr_start})+1;
      $exon_content_chr{$exon_id}{size} = $size;
      $genes_ref->{$gene_id}->{base_count} += $size;

      #Store the exon-content coordinates
      $exon_content_ref->{$exon_id}->{chr_start} = $coords_ref->{$gene_id}->{chr_start};
      $exon_content_ref->{$exon_id}->{chr_end} = $coords_ref->{$gene_id}->{chr_end};

    }

    #Build a coverage hash for this gene (keyed on chromosome positions)
    #Use the chromosome coordinates for all exon-content coordinates
    foreach my $exon_id (sort {$exon_content_chr{$a}->{chr_start} <=> $exon_content_chr{$b}->{chr_start}} keys %exon_content_chr){

      my $start = $exon_content_chr{$exon_id}{chr_start};
      my $end = $exon_content_chr{$exon_id}{chr_end};

      #Go through each chromosome position in this exon content block and initialize that position in the hash
      for (my $i = $start; $i <= $end; $i++){
	my $gene_pos_id = "$gene_id"."_"."$i";
	$gene_coverage{$gene_pos_id} = 0;
      }
    }
    %exon_content_chr = ();

    my $transcripts_ref = $gene_transcripts_ref->{$gene_id}->{transcripts};

    #Assemble a reference set of exons (a superset of all non-redundant exons)
    my $trans_count = keys %{$transcripts_ref};
    $gene_transcripts_ref->{$gene_id}->{trans_count} = $trans_count;

    my %reference_exons;
    foreach my $trans_id (keys %{$transcripts_ref}){
      my $exons_ref = $transcripts_ref->{$trans_id}->{exons};

      foreach my $trans_exon_id (keys %{$exons_ref}){
        my $trans_exon_start = $exons_ref->{$trans_exon_id}->{exon_start};
        my $trans_exon_end = $exons_ref->{$trans_exon_id}->{exon_end};

        #Check each of the reference exons to see if one of these is the same, otherwise add it to the list
        my $redundant_exon = 0;
        foreach my $ref_exon_id (keys %reference_exons){
	  if ($trans_exon_start == $reference_exons{$ref_exon_id}{exon_start} && $trans_exon_end == $reference_exons{$ref_exon_id}{exon_end}){
	    $redundant_exon = 1;
	  }
        }
        #Unless the current transcript exon was found to be redundant, add it to the list
        unless ($redundant_exon == 1){
	  $reference_exons{$trans_exon_id}{exon_start} = $trans_exon_start;
	  $reference_exons{$trans_exon_id}{exon_end} = $trans_exon_end;
        }
      }
    }

    #Get arrays to represent the reference exons
    my $ref_exon_count = keys %reference_exons;

    my @reference_exon_starts;
    my @reference_exon_ends;
    foreach my $ref_exon_id (sort {$reference_exons{$a}->{exon_start} <=> $reference_exons{$b}->{exon_start}} keys %reference_exons){
      push (@reference_exon_starts, $reference_exons{$ref_exon_id}{exon_start});
      push (@reference_exon_ends, $reference_exons{$ref_exon_id}{exon_end});
    }

    my %starts;    #Non-redundant list of known start positions for all transcripts of the current gene
    my %ends;      #Non-redundant list of known end positions for all transcripts of the current gene
    foreach my $trans_id (keys %{$transcripts_ref}){
      my $exons_ref = $transcripts_ref->{$trans_id}->{exons};

      my @starts;
      my @ends;

      #Build arrays of the exon starts and ends for this transcript
      foreach my $trans_exon_id (sort {$exons_ref->{$a}->{exon_start} <=> $exons_ref->{$b}->{exon_start}} keys %{$exons_ref}){
        my $trans_exon_start = $exons_ref->{$trans_exon_id}->{exon_start};
        my $trans_exon_end = $exons_ref->{$trans_exon_id}->{exon_end};
        push (@starts, $trans_exon_start);
        push (@ends, $trans_exon_end);

        $starts{$trans_exon_start}{label} = "na";
        $ends{$trans_exon_end}{label} = "na";
      }
    }

    #Now assign 'exon order' to each exon content block
    my $exon_count = keys %{$exon_content_ref};

    my $exon_order = 0;

    foreach my $exon_id (sort {$exon_content_ref->{$a}->{start} <=> $exon_content_ref->{$b}->{start}} keys %{$exon_content_ref}){
      my $ec_start = $exon_content_ref->{$exon_id}->{start};
      my $ec_end = $exon_content_ref->{$exon_id}->{end};

      #Note that since we are dealing with GENE coordinates which are always stored 5' to 3', the strand does no matter at this stage
      $exon_order++;
      $exon_content_ref->{$exon_id}->{exon_order} = $exon_order;

      #For this exon content block, identify the exon starts and ends that are within it:
      #Then label the donors (exon ends) and acceptor (exon starts) from A to Z (5' to 3')

      my %ec_starts; #Starts that are actually within the current exon content block
      my %ec_ends;   #Ends that are actually within the current exon content block

      foreach my $start (keys %starts){
        if ($start >= $ec_start && $start <= $ec_end){
	  $ec_starts{$start}{tmp} = "na";
        }
      }
      foreach my $end (keys %ends){
        if ($end >= $ec_start && $end <= $ec_end){
	  $ec_ends{$end}{tmp} = "na";
        }
      }

      #Now apply the labels to each Acceptor or Donor
      my @labels = @letters;
      foreach my $ec_start (sort {$a <=> $b} keys %ec_starts){
        $starts{$ec_start}{label} = shift (@labels);
      }
      @labels = @letters;
      foreach my $ec_end (sort {$a <=> $b} keys %ec_ends){
        $ends{$ec_end}{label} = shift (@labels);
      }
    }

    #Now annotate the actual exons that fall within each exon content block
    my $acceptors_ref = \%starts;
    my $donors_ref = \%ends;
    my $ref_exons = \%reference_exons;


    #Go through all the non-redundant exons of this gene
    foreach my $ref_exon_id (sort {$ref_exons->{$a}->{exon_start} <=> $ref_exons->{$b}->{exon_start}} keys %{$ref_exons}){

      my $ref_exon_start = $ref_exons->{$ref_exon_id}->{exon_start};
      my $ref_exon_end = $ref_exons->{$ref_exon_id}->{exon_end};
      my $ref_exon_order = "na";

      #Now find the exon content block this exon belongs to
      foreach my $exon_content_id (sort {$exon_content_ref->{$a}->{start} <=> $exon_content_ref->{$b}->{start}} keys %{$exon_content_ref}){
        my $ec_exon_start = $exon_content_ref->{$exon_content_id}->{start};
        my $ec_exon_end = $exon_content_ref->{$exon_content_id}->{end};
        my $ec_exon_order = $exon_content_ref->{$exon_content_id}->{exon_order};

        #Does the current exon overlap this exon content block
        if ($ref_exon_start >= $ec_exon_start && $ref_exon_start <= $ec_exon_end && $ref_exon_end >= $ec_exon_start && $ref_exon_end <= $ec_exon_end){
	  $ref_exon_order = $ec_exon_order;
        }
      }

      #Now that the exon name (i.e. exon order from 5' to 3') has been found.  Get the acceptor and donor labels
      my $acceptor_label = $acceptors_ref->{$ref_exon_start}->{label};
      my $donor_label = $donors_ref->{$ref_exon_end}->{label};
      unless (defined($ref_exon_order)){$ref_exon_order=0;}
      unless (defined($donor_label)){$donor_label = "x";}
      unless (defined($acceptor_label)){$acceptor_label = "x";}

      my $ref_exon_name = "$acceptor_label"."_E"."$ref_exon_order"."_"."$donor_label"."__"."$gene_id";

      #Get chromosome coordinates for this reference exon
      my $coords_ref = &convertGeneCoordinates ('-gene_object'=>$genes_ref, '-gene_id'=>$gene_id, '-start_pos'=>$ref_exon_start, '-end_pos'=>$ref_exon_end, '-ordered'=>"yes");

      my $query_chr_start = $coords_ref->{$gene_id}->{chr_start};
      my $query_chr_end = $coords_ref->{$gene_id}->{chr_end};
      my $query_strand = $coords_ref->{$gene_id}->{strand};

      #NOW create the UCSC record
      my $db;
      if ($database =~ /ALEXA_(\w+)/){
        $db = $1;
      }
      $ref_exon_counter++;
      my $record = "chr$chromosome\t$db\tExon_Annotation\t$query_chr_start\t$query_chr_end\t.\t$query_strand\t.\t$ref_exon_name\n";
      $ensembl_exon_annotation_tracks{$ref_exon_counter} = $record;
    }
  }

  #Close database connection
  $alexa_dbh->disconnect();

  return();
}


#########################################################################################################################################
#build Chromosome and Gene Coverage Objects                                                                                             #
#########################################################################################################################################
sub buildCoverageObjects{
  my %args = @_;
  my $target_chr = $args{'-chromosome'};
  my $range = $args{'-range'};
  my $working_dir = $args{'-working_dir'};

  $| = 1; print BLUE, "\n\n4.) Building GENE and CHROMOSOME coverage objects for chr$target_chr:$range", RESET; $| = 0;
  print LOG "\n\n4.) Building GENE and CHROMOSOME coverage objects for chr$target_chr:$range";

  #Open the previously created berkeley DB object containing the reads of this chromosome - open in readonly mode
  $target_chr = "chr"."$target_chr";

  my $counter = 0;
  my $read_counter = 0;
  my $block_count = 0;
  my $block_size = 10000;

  #while (my ($read_id) = each %chr_reads){
  while (my ($string) = each %chr_reads){

    #Data format:  $gene_id\t $chr_start_coords\t $chr_end_coords
    #my $string = $chr_reads{$read_id};
    my $string_count = $chr_reads{$string};
    my @data = split("\t", $string);
    my $gene_id = $data[0];
    my @chr_start_coords = split(" ", $data[1]);
    my @chr_end_coords = split(" ", $data[2]);

    for (my $i = 1; $i <= $string_count; $i++){
      $counter++;
      $read_counter++;
      if ($counter == $block_size){
        $block_count++;
        $| = 1; print BLUE, "\n\tProcessed block $block_count (of size $block_size reads for $target_chr: $range) - $read_counter reads so far", RESET; $| = 0;
        print LOG "\n\tProcessed block $block_count (of size $block_size reads for $target_chr: $range) - $read_counter reads so far";
        $counter = 0;
      }

      #Add coverage for gene
      &addReadGeneCoverage('-gene_id'=>$gene_id, '-chr_starts'=>\@chr_start_coords, '-chr_ends'=>\@chr_end_coords);

      #Add coverage for chromosome
      &addReadChrCoverage('-chr_starts'=>\@chr_start_coords, '-chr_ends'=>\@chr_end_coords);
    }
  }
  return();
}


############################################################################################################################################
#Add the coverage of a read to a gene to the coverage hash for the exon content record for that gene
############################################################################################################################################
sub addReadGeneCoverage{
  my %args = @_;
  my $gene_id = $args{'-gene_id'};
  my @chr_starts = @{$args{'-chr_starts'}};
  my @chr_ends = @{$args{'-chr_ends'}};

  $genes_ref->{$gene_id}->{quality_read_count}++;

  my @ends = @chr_ends;
  foreach my $start (@chr_starts){
    my $end = shift(@ends);

    #Go through each chromosome position in this read as it is mapped to an exon and increment that position in the hash
    for (my $i = $start; $i <= $end; $i++){

      my $gene_pos_id = "$gene_id"."_"."$i";

      if ($gene_coverage{$gene_pos_id} >= 0){
	$gene_coverage{$gene_pos_id}++;
      }else{
	print RED, "\nGene_position record ($gene_pos_id) was not initialized and may not correspond to known exonic position!!", RESET;
	exit();
      }
    }
  }
  return();
}


############################################################################################################################################
#Add the coverage of a read to a gene to the coverage hash for the exon content record for the entire chromosome
############################################################################################################################################
sub addReadChrCoverage{
  my %args = @_;
  my @chr_starts = @{$args{'-chr_starts'}};
  my @chr_ends = @{$args{'-chr_ends'}};

  my @ends = @chr_ends;
  foreach my $start (@chr_starts){
    my $end = shift(@ends);

    #Go through each chromosome position in this read as it is mapped to an exon and increment that position in the hash
    for (my $i = $start; $i <= $end; $i++){

      my $pos_id = "$i";

      if ($chr_coverage{$pos_id}){
	$chr_coverage{$pos_id}++;
      }else{
	$chr_coverage{$pos_id} = 1;
      }
    }
  }
  return();
}


############################################################################################################################################
#Correct base count values at the GENE and CHROMOSOME level - correct for library size and mappability                                     #
############################################################################################################################################
sub correctBaseCounts{
  my %args = @_;
  my $normalization_value = $args{'-normalization_value'};  #User specified normalization value (# quality mapped reads from smallest library)
  my $library_size = $args{'-library_size'};                #Actual library size according to input read file

  $| = 1; print BLUE, "\n\n5.) Correct GENE and CHROMOSOME base counts for chr$chr_filter: $start_filter-$end_filter", RESET; $| = 0;
  print LOG "\n\n5.) Correct GENE and CHROMOSOME base counts for chr$chr_filter: $start_filter-$end_filter";

  my $library_correction_ratio = 1;
  unless ($library_size == 0){
    $library_correction_ratio = $normalization_value/$library_size;
  }

  #A.) CHROMOSOME BASE CORRECTION
  #Correct for library size and mapability - store values in their respective hashes
  while (my ($pos) = each %chr_coverage){

    #First correct for library size
    my $coverage_val = $chr_coverage{$pos};
    my $corrected_val = $coverage_val*$library_correction_ratio;
    my $corrected_val_f = sprintf("%.10f", $corrected_val);
    $chr_coverage_norm1{$pos} = $corrected_val_f;
  }


  #B.) GENE BASE CORRECTION
  while (my ($gene_pos) = each %gene_coverage){

    if ($gene_pos =~ /(\d+)\_(\d+)/){
      my $gene_id = $1;
      my $pos = $2;

      #First correct for library size
      my $coverage_val = $gene_coverage{$gene_pos};
      my $corrected_val = $coverage_val*$library_correction_ratio;
      my $corrected_val_f = sprintf("%.10f", $corrected_val);
      $gene_coverage_norm1{$gene_pos} = $corrected_val_f;

    }else{
      print RED, "\nGene_Pos ID not understood!\n\n", RESET;
      exit();
    }
  }
  return();
}


############################################################################################################################################
#6.) Import the genes annotation file 
############################################################################################################################################
sub importGenes{
  my %args = @_;
  my $gene_file = $args{'-infile'};
  my $current_chromosome = $args{'-chromosome'};
  my $range = $args{'-range'};

  my $start_filter;
  my $end_filter;
  if ($range =~ /(\d+)\-(\d+)/){
    $start_filter = $1;
    $end_filter = $2;
  }else{
    print RED, "\nProblem with range - importExonRegions()\n\n", RESET;
    exit();
  } 

  open (GENE, "zcat $gene_file |") || die "\nCould not open gene annotation file: $gene_file\n\n";  

  my $header = 1;
  my %columns;

  while(<GENE>){
    chomp($_);
    my $line = $_;
    my @line = split("\t", $line);

    if ($header == 1){
      $gene_header = $_;
      my $column_count = 0;
      foreach my $column (@line){
        $columns{$column}{column_pos} = $column_count;
        $column_count++;
      }
      $header = 0;
      next();
    }
    my $chr = $line[$columns{'Chromosome'}{column_pos}];
    my $gene_id = $line[$columns{'Gene_ID'}{column_pos}];
    my $start_chr = $line[$columns{'Unit1_start_chr'}{column_pos}];
    my $end_chr = $line[$columns{'Unit1_end_chr'}{column_pos}];

    if ($chr eq "MT"){
      $chr = "M";
    }

    #Only store lines corresponding to the chromosome region being processed
    unless ($chr eq $current_chromosome && ($start_chr >= $start_filter && $start_chr <= $end_filter && $end_chr >= $start_filter && $end_chr <= $end_filter)){
      next();
    }

    my $base_count = $line[$columns{'Base_Count'}{column_pos}];
    my $unmasked_base_count = $line[$columns{'UnMasked_Base_Count'}{column_pos}];
    my $coding_base_count = $line[$columns{'Coding_Base_Count'}{column_pos}];

    $genes_ref->{$gene_id}->{base_count} = $base_count;
    $genes_ref->{$gene_id}->{unmasked_base_count} = $unmasked_base_count;
    $genes_ref->{$gene_id}->{coding_base_count} = $coding_base_count;
    $genes_ref->{$gene_id}->{line} = $line;

  }
  close (GENE);

  return();
}


############################################################################################################################################
#6.) Import the exon regions database file 
############################################################################################################################################
sub importExonRegions{
  my %args = @_;
  my $exon_region_file = $args{'-infile'};
  my $current_chromosome = $args{'-chromosome'};
  my $range = $args{'-range'};

  my $start_filter;
  my $end_filter;
  if ($range =~ /(\d+)\-(\d+)/){
    $start_filter = $1;
    $end_filter = $2;
  }else{
    print RED, "\nProblem with range - importExonRegions()\n\n", RESET;
    exit();
  } 

  my %exon_regions;

  open (ER, "zcat $exon_region_file |") || die "\nCould not open exon region file: $exon_region_file\n\n";  

  my $er_count = 0;
  my $header = 1;
  my %columns;

  while(<ER>){
    chomp($_);
    my $line = $_;
    my @line = split("\t", $line);

    if ($header == 1){
      $er_header = $_;
      my $column_count = 0;
      foreach my $column (@line){
        $columns{$column}{column_pos} = $column_count;
        $column_count++;
      }
      $header = 0;
      next();
    }
    $er_count++;
    my $chr = $line[$columns{'Chromosome'}{column_pos}];
    my $gene_id = $line[$columns{'Gene_ID'}{column_pos}];
    my $start_chr = $line[$columns{'Unit1_start_chr'}{column_pos}];
    my $end_chr = $line[$columns{'Unit1_end_chr'}{column_pos}];
    my $strand = $line[$columns{'Strand'}{column_pos}];
    my $exon_region_name = $line[$columns{'Seq_Name'}{column_pos}];

    #Only store lines corresponding to the chromosome region being processed
    unless ($chr eq $current_chromosome && ($start_chr >= $start_filter && $start_chr <= $end_filter && $end_chr >= $start_filter && $end_chr <= $end_filter)){
      next();
    }


    if ($strand eq "1"){
      $strand = "+";
    }else{
      $strand = "-";
    }
    $exon_regions{$er_count}{line} = $line;
    $exon_regions{$er_count}{gene_id} = $gene_id;
    $exon_regions{$er_count}{name} = $exon_region_name;
    $exon_regions{$er_count}{chromosome} = $chr;
    $exon_regions{$er_count}{strand} = $strand;
    $exon_regions{$er_count}{start} = $start_chr;
    $exon_regions{$er_count}{end} = $end_chr;
    $exon_regions{$er_count}{expressed} = 0;

    my $ucsc_chromosome = "chr"."$chr";
    my $record_id = "$exon_region_name"."_"."$gene_id";
    my $record = "\n$ucsc_chromosome\tEnsEMBL\texon_region\t$start_chr\t$end_chr\t.\t$strand\t.\t$record_id";
    $ensembl_exon_region_tracks{$er_count} = $record;

  }
  close (ER);

  return(\%exon_regions);
}



############################################################################################################################################
#Print out Gene, Exon and Position-Bias Summary files
############################################################################################################################################
sub printGeneExonSummary{
  my %args = @_;
  my $results_dir = $args{'-results_dir'};

  $| = 1; print BLUE, "\n\n6.) Printing Gene Summary data to: $results_dir", RESET; $| = 0;
  $| = 1; print BLUE, "\n\t6-a.) Summarizing sequence coverage of each gene ...\n", RESET; $| = 0;
  print LOG "\n\n6.) Printing Gene Summary data to: $results_dir";
  print LOG "\n\t6-a.) Summarizing sequence coverage of each gene ...\n";


  #A.) CALCULATE CUMULATIVE GENE COVERAGE VALUES
  #The coverage of each base of the exon content of each gene has been stored for each gene
  #Use this information to calculate the sequence coverage of each gene
  #Summarize the % of bases with at least 1X coverage and the overall average X coverage
  while (my ($gene_pos) = each %gene_coverage){
    my $coverage_raw = $gene_coverage{$gene_pos};
    my $coverage_norm1 = $gene_coverage_norm1{$gene_pos};

    my $gene_id;
    my $pos;
    if ($gene_pos =~ /(\d+)\_(\d+)/){
      $gene_id = $1;
      $pos = $2;
    }

    $genes_ref->{$gene_id}->{cumulative_coverage_raw} += $coverage_raw;
    $genes_ref->{$gene_id}->{cumulative_coverage_norm1} += $coverage_norm1;

    if ($coverage_raw >= 1){$genes_ref->{$gene_id}->{bases_covered_1x}++};
    if ($coverage_raw >= 5){$genes_ref->{$gene_id}->{bases_covered_5x}++};
    if ($coverage_raw >= 10){$genes_ref->{$gene_id}->{bases_covered_10x}++};
    if ($coverage_raw >= 50){$genes_ref->{$gene_id}->{bases_covered_50x}++};
    if ($coverage_raw >= 100){$genes_ref->{$gene_id}->{bases_covered_100x}++};
    if ($coverage_raw >= 500){$genes_ref->{$gene_id}->{bases_covered_500x}++};
  }

  $| = 1; print BLUE, "\n\t6-b.) Now printing the actual Gene and Exon summary lines for chr: $chr_filter: $start_filter-$end_filter ...\n", RESET; $| = 0;
  print LOG "\n\t6-b.) Now printing the actual Gene and Exon summary lines for chr: $chr_filter: $start_filter-$end_filter ...\n";

  my $gene_outfile = "$results_dir"."chr"."$file_name_prefix"."_GeneExpression.txt";
  my $exon_outfile = "$results_dir"."chr"."$file_name_prefix"."_ExonRegionExpression.txt";

  #Open the output file
  open (GENE_OUT, ">$gene_outfile") || die "\nCould not open output gene summary file: $gene_outfile\n\n";
  print GENE_OUT "$gene_header\tRead_Count\tCumulative_Coverage\tAverage_Coverage_RAW\tAverage_Coverage_NORM1\tBases_Covered_1x\tPercent_Coverage_1x\tPercent_Coverage_5x\tPercent_Coverage_10x\tPercent_Coverage_50x\tPercent_Coverage_100x\tPercent_Coverage_500x\tExpressed\tLink\n";

  unless ($web_path =~ /.*\/$/){
    $web_path = "$web_path"."/";
  }

  #Go through each gene and print out basic data for it.  Also include a link to the custom track file for its chromosome
  foreach my $gene_id (sort {$a <=> $b} keys %{$genes_ref}){

    my $base_count = $genes_ref->{$gene_id}->{base_count};
    my $bases_covered_1x = $genes_ref->{$gene_id}->{bases_covered_1x};
    my $bases_covered_5x = $genes_ref->{$gene_id}->{bases_covered_5x};
    my $bases_covered_10x = $genes_ref->{$gene_id}->{bases_covered_10x};
    my $bases_covered_50x = $genes_ref->{$gene_id}->{bases_covered_50x};
    my $bases_covered_100x = $genes_ref->{$gene_id}->{bases_covered_100x};
    my $bases_covered_500x = $genes_ref->{$gene_id}->{bases_covered_500x};
    my $cumulative_coverage_raw = $genes_ref->{$gene_id}->{cumulative_coverage_raw};
    my $cumulative_coverage_norm1 = $genes_ref->{$gene_id}->{cumulative_coverage_norm1};

    my $percent_coverage_1x = sprintf("%.2f", (($bases_covered_1x/$base_count)*100));
    my $percent_coverage_5x = sprintf("%.2f", (($bases_covered_5x/$base_count)*100));
    my $percent_coverage_10x = sprintf("%.2f", (($bases_covered_10x/$base_count)*100));
    my $percent_coverage_50x = sprintf("%.2f", (($bases_covered_50x/$base_count)*100));
    my $percent_coverage_100x = sprintf("%.2f", (($bases_covered_100x/$base_count)*100));
    my $percent_coverage_500x = sprintf("%.2f", (($bases_covered_500x/$base_count)*100));
    my $average_coverage_raw = sprintf("%.10f", ($cumulative_coverage_raw/$base_count));
    my $average_coverage_norm1 = sprintf("%.10f", ($cumulative_coverage_norm1/$base_count));

    #Determine whether this element should be considered as expressed above background
    my $cutoff_test = 1;
    if ($cutoffs_file){
      my @result = @{&testExpression('-cutoffs_ref'=>$gene_cutoffs_ref, '-gene_id'=>$gene_id, '-norm_expression_value'=>$average_coverage_norm1, '-raw_expression_value'=>$average_coverage_raw, '-percent_gene_expression_cutoff'=>0)};
      $cutoff_test = $result[0];
    }

    #If this exon region meets the criteria for expression store a UCSC track record for it...
    my $expressed = 0;
    if ($percent_coverage_1x >= 50.0 && $cutoff_test == 1){
      $expressed = 1;
    }

    #Create a link to go directly to the region of this gene (+/- 100bp) and load the correct chromosome file
    my $display_start = $genes_ref->{$gene_id}->{chr_start};
    my $display_end = $genes_ref->{$gene_id}->{chr_end};
    my $temp;
    if ($display_start > $display_end){
      $temp = $display_start;
      $display_start = $display_end;
      $display_end = $temp;
    }
    $display_start -= 100;
    $display_end += 100;

    my $chromosome = "chr"."$genes_ref->{$gene_id}->{chromosome}";
    if ($chromosome eq "chrMT"){
      $chromosome = "chrM";
    }

    #Link that retains pre-existing custom tracks
    #Clean link
    my $link = "http://genome.ucsc.edu/cgi-bin/hgTracks?db="."$ucsc_build"."&position=$chromosome:$display_start-$display_end&hgt.customText=$web_path"."chr$file_name_prefix"."_geneExon.txt.gz&ctfile_"."$ucsc_build"."=";

    my $line = $genes_ref->{$gene_id}->{line};
    print GENE_OUT "$line\t$genes_ref->{$gene_id}->{quality_read_count}\t$cumulative_coverage_raw\t$average_coverage_raw\t$average_coverage_norm1\t$bases_covered_1x\t$percent_coverage_1x\t$percent_coverage_5x\t$percent_coverage_10x\t$percent_coverage_50x\t$percent_coverage_100x\t$percent_coverage_500x\t$expressed\t$link\n";
  }
  close(GENE_OUT);


  #B.) CALCULATE CUMULATIVE EXON REGION COVERAGE VALUES - And print to output file

  #Open the output file
  open (EXON_OUT, ">$exon_outfile") || die "\nCould not open output exon summary file: $exon_outfile\n\n";
  print EXON_OUT "$er_header\tCumulative_Coverage\tAverage_Coverage_RAW\tAverage_Coverage_NORM1\tBases_Covered_1x\tPercent_Coverage_1x\tPercent_Coverage_5x\tPercent_Coverage_10x\tPercent_Coverage_100x\tExpressed\tPercent_Gene_Expression\n";

  foreach my $er_id (sort {$a <=> $b} keys %{$exon_regions_ref}){

    #Traverse through the base positions of this exon content block and calculate the cumulative coverage
    #Also note the # of bases covered at 1X or greater depth, calculate percent coverage of the exon region, and if this exceeds the cutoff, store this exon region as 'expressed'
    my $bases_covered_1x = 0;
    my $bases_covered_5x = 0;
    my $bases_covered_10x = 0;
    my $bases_covered_100x = 0;
    my $cumulative_coverage_raw = 0;
    my $cumulative_coverage_norm1 = 0;

    my $gene_id = $exon_regions_ref->{$er_id}->{gene_id};
    my $start = $exon_regions_ref->{$er_id}->{start};
    my $end = $exon_regions_ref->{$er_id}->{end};
    my $exon_size = ($end - $start)+1;

    for (my $i = $start; $i <= $end; $i++){
      my $gene_pos = "$gene_id"."_"."$i";
      $cumulative_coverage_raw += $gene_coverage{$gene_pos};
      $cumulative_coverage_norm1 += $gene_coverage_norm1{$gene_pos};

      if ($gene_coverage{$gene_pos} >= 1){$bases_covered_1x++};
      if ($gene_coverage{$gene_pos} >= 5){$bases_covered_5x++};
      if ($gene_coverage{$gene_pos} >= 10){$bases_covered_10x++};
      if ($gene_coverage{$gene_pos} >= 100){$bases_covered_100x++};
    }

    my $percent_coverage_1x = sprintf("%.2f", (($bases_covered_1x/$exon_size)*100));
    my $percent_coverage_5x = sprintf("%.2f", (($bases_covered_5x/$exon_size)*100));
    my $percent_coverage_10x = sprintf("%.2f", (($bases_covered_10x/$exon_size)*100));
    my $percent_coverage_100x = sprintf("%.2f", (($bases_covered_100x/$exon_size)*100));

    my $average_coverage_raw = sprintf("%.10f", ($cumulative_coverage_raw/$exon_size));
    my $average_coverage_norm1 = sprintf("%.10f", ($cumulative_coverage_norm1/$exon_size));

    #Determine whether this element should be considered as expressed above background
    my $cutoff_test = 1;
    my $percent_gene_expression = 0;
    if ($cutoffs_file){
      my @result = @{&testExpression('-cutoffs_ref'=>$gene_cutoffs_ref, '-gene_id'=>$gene_id, '-norm_expression_value'=>$average_coverage_norm1, '-raw_expression_value'=>$average_coverage_raw, '-percent_gene_expression_cutoff'=>20)};
      $cutoff_test = $result[0];
      $percent_gene_expression = $result[1];
    }

    #If this exon region meets the criteria for expression store a UCSC track record for it...
    my $expressed = 0;
    if ($percent_coverage_1x >= $min_seq_coverage && $cutoff_test == 1){
      $expressed = 1;
      my $gene_id = $exon_regions_ref->{$er_id}->{gene_id};
      my $exon_region_name = $exon_regions_ref->{$er_id}->{name};
      my $chr = $exon_regions_ref->{$er_id}->{chromosome};
      my $strand = $exon_regions_ref->{$er_id}->{strand};
      my $start_chr = $exon_regions_ref->{$er_id}->{start};
      my $end_chr = $exon_regions_ref->{$er_id}->{end};
      my $ucsc_chromosome = "chr"."$chr";
      my $record_id = "$exon_region_name"."_"."$gene_id";
      my $record = "\n$ucsc_chromosome\tEnsEMBL\texon_region\t$start_chr\t$end_chr\t.\t$strand\t.\t$record_id";
      $expressed_exon_region_tracks{$er_id} = $record;
    }

    print EXON_OUT "$exon_regions_ref->{$er_id}->{line}\t$cumulative_coverage_raw\t$average_coverage_raw\t$average_coverage_norm1\t$bases_covered_1x\t$percent_coverage_1x\t$percent_coverage_5x\t$percent_coverage_10x\t$percent_coverage_100x\t$expressed\t$percent_gene_expression\n";
  }

  close(EXON_OUT);
  return();
}


############################################################################################################################################
#Print out UCSC track files
############################################################################################################################################
sub printUCSCTrack{
  my %args = @_;
  my $ucsc_dir = $args{'-ucsc_dir'};
  my $target_chr = $args{'-chr'};
  my $library = $args{'-library'};

  $target_chr = "chr"."$target_chr";

  unless ($ucsc_dir =~ /.*\/$/){
    $ucsc_dir = "$ucsc_dir"."/";
  }

  my $ucsc_file = "$ucsc_dir"."chr$file_name_prefix"."_geneExon.txt";
  #Print out all the UCSC track records gathered above for this chromosome
  open (UCSC, ">$ucsc_file") || die "\nCould not open ucsc file: $ucsc_file\n\n";

  $| = 1; print BLUE, "\n\n7.) Printing UCSC file for $target_chr: $ucsc_file", RESET; $| = 0;
  print LOG "\n\n7.) Printing UCSC file for $target_chr: $ucsc_file";

  #Browser line
  print UCSC "#Browser line";
  print UCSC "\nbrowser hide all";
  print UCSC "\nbrowser full knownGene";
  print UCSC "\nbrowser pack multiz28way";

  my $database_abr = $database;
  if ($database =~ /ALEXA_(\w+)/){
    $database_abr = $1;
  }

  #ANNOTATION TRACKS
  #1.) Track line for ANNOTATION of EnsEMBL transcripts of each gene
  print UCSC "\n\n#EnsEMBL transcripts";
  my $trans_track_name = "Transcripts";
  my $trans_track_description = "\"EnsEMBL Transcripts ($database_abr)\"";
  $priority_annotation++;

  print UCSC "\ntrack name=$trans_track_name description=$trans_track_description color=$colors{$color_set}{0} useScore=0 visibility=3 priority=$priority_annotation";
  print UCSC "\n\n#Begin DATA\n";

  foreach my $trans_record (keys %ensembl_transcript_tracks){
    my $record = $ensembl_transcript_tracks{$trans_record};
    print UCSC "$record";
  }

  #2.) Track line for ANNOTATION of EnsEMBL Exons and their Acceptor/Donor sites
  print UCSC "\n\n#Non redundant EnsEMBL exon annotations $database";
  my $exons_track_name = "Exons_NR";
  my $exons_track_description = "\"Annotations of Non-Redundant EnsEMBL Exons ($database_abr)\"";
  $priority_annotation++;
  print UCSC "\ntrack name=$exons_track_name description=$exons_track_description color=$colors{$color_set}{1} useScore=0 visibility=3 priority=$priority_annotation";
  print UCSC "\n\n#Begin DATA\n";

  foreach my $e_id (sort {$a <=> $b} keys %ensembl_exon_annotation_tracks){
    my $record = $ensembl_exon_annotation_tracks{$e_id};
    print UCSC "$record";
  }

  #3.) Track line for ANNOTATION of EnsEMBL Exon 'Regions'
  print UCSC "\n\n#EnsEMBL exon regions";
  my $er_track_name = "ExonRegions";
  my $er_track_description = "\"EnsEMBL Exon Regions ($database_abr)\"";
  $priority_annotation++;
  print UCSC "\ntrack name=$er_track_name description=$er_track_description color=$colors{$color_set}{0} useScore=0 visibility=3 priority=$priority_annotation";
  print UCSC "\n\n#Begin DATA\n";

  foreach my $er_id (sort {$a <=> $b} keys %ensembl_exon_region_tracks){
    my $record = $ensembl_exon_region_tracks{$er_id};
    print UCSC "$record";
  }


  #EXPRESSED REGION TRACKS
  #4.) Track line for EXPRESSED Exon 'Regions'
  # - Exon regions which are covered over X% of their bases or greater (where X is specified by the user: say 75% as --min_seq_coverage)
  print UCSC "\n\n#EnsEMBL Expressed exon regions";
  my $eer_track_name = "$library"."_Exp_ER";
  my $eer_track_description = "\"Expressed EnsEMBL Exon Regions (>= $min_seq_coverage% coverage and above intronic background cutoff) ($library_name - $database_abr)\"";
  $priority_expression++;
  print UCSC "\ntrack name=$eer_track_name description=$eer_track_description color=$colors{$color_set}{2} useScore=0 visibility=3 priority=$priority_expression";
  print UCSC "\n\n#Begin DATA\n";

  foreach my $er_id (sort {$a <=> $b} keys %expressed_exon_region_tracks){
    my $record = $expressed_exon_region_tracks{$er_id};
    print UCSC "$record";
  }


  #WIGGLE TRACKS
  my $wig_track_name;
  my $wig_track_description;
  my $current_pos;

  #5.) Create a WIG track line to display the RAW coverage level for every exonic base sequence to 1X or greater depth
  #print UCSC "\n\n#WIG TRACK: RAW Exonic Coverage calculated from Illumina Paired Reads Mapped to EnsEMBL Transcripts";
  #$wig_track_name = "$library"."_RAW_E";
  #$wig_track_description = "\"RAW Exonic Coverage. Paired Reads mapped to EnsEMBL ($library_name - $database_abr)\"";
  #$priority_wiggle++;

  #print UCSC "\ntrack name=$wig_track_name description=$wig_track_description type=wiggle_0 color=$colors{$color_set}{2} yLineMark=0.0 yLineOnOff=on visibility=hide autoScale=on graphType=bar smoothingWindow=off maxHeightPixels=120:80:10 priority=$priority_expression";
  #print UCSC "\n\n#Begin DATA\n";

  #$current_pos = -1;
  #foreach my $pos (sort {$a <=> $b} keys %chr_coverage){
    #If this is a new block of covered bases, print out a new def line
    #if ($pos == ($current_pos+1)){
      #print UCSC "\n$chr_coverage{$pos}";
    #}else{
      #my $pos_minus_1 = $pos-1;
      #print UCSC "\nfixedStep chrom=$target_chr start=$pos_minus_1 step=1";
      #print UCSC "\n0.00"; #Start each block with a 0
      #print UCSC "\n$chr_coverage{$pos}";
    #}
  #$current_pos = $pos;
  #}


  #6.) Create a WIG track line to display the NORM1 coverage level for every exonic base sequence to 1X or greater depth
  #    - NORM1 = library size adjusted values
  print UCSC "\n\n#WIG TRACK: NORM1 Exonic Coverage calculated from Illumina Paired Reads Mapped to EnsEMBL Transcripts";
  $wig_track_name = "$library"."_N1_E";
  $wig_track_description = "\"Normalized Exonic Coverage. Paired reads mapped to EnsEMBL Transcripts ($library_name - $database_abr)\"";
  $priority_wiggle++;

  print UCSC "\ntrack name=$wig_track_name description=$wig_track_description type=wiggle_0 color=$colors{$color_set}{2} yLineMark=0.0 yLineOnOff=on visibility=full autoScale=on graphType=bar smoothingWindow=off maxHeightPixels=120:80:10 priority=$priority_expression";
  print UCSC "\n\n#Begin DATA\n";

  $current_pos = -1;
  foreach my $pos (sort {$a <=> $b} keys %chr_coverage_norm1){
    my $val = $chr_coverage_norm1{$pos};

    #If this is a new block of covered bases, print out a new def line
    if ($pos == ($current_pos+1)){
      print UCSC "\n$val";
    }else{
      my $pos_minus_1 = $pos-1;
      print UCSC "\nfixedStep chrom=$target_chr start=$pos_minus_1 step=1";
      print UCSC "\n0.00"; #Start each block with a 0
      print UCSC "\n$val";
    }
    $current_pos = $pos;
  }

  close(UCSC);

  #Finally gzip the file
  my $cmd = "gzip -f $ucsc_file";
  system($cmd);

  return();
}



