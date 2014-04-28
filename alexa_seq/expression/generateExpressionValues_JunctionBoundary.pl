#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#The purpose of this script is to summarize the expression of exon-exon junctions OR exon boundaries (known and theoretical) from an EnsEMBL database
#Note that a boundary here means an extension of a known exon boundary (which in rare cases is already represented by another known, overlaping exon)

#The inputs for this script are reads mapped to individual exon-exon juctions or exon boundaries (one individual read per line)
#The outputs are a summary file describing the observed expression of all possible junctions/boundaries and, UCSC GFF custom track files depicting the junctions/boundaries observed

#Summary done for every junction/boundary:
#  - Join expression information onto the complete, annotated, junction/boundary database file:
#  - Note the number of reads hitting each junction/boundary
#  - Calculate the average coverage for each junction/boundary
#    - Even though all junctions/boundaries in the database are the same size it is still nice to have average coverage values for comparison to exons, genes, etc.
#  - Calculate the number of bases of the junction/boundary covered by one or more reads (give some sense of the diversity of reads supporting a junction/boundary)
#  - Correct the read count and coverage values to a user specified normalization value to allow direct comparison between libraries

#Create a gene level summary
#  - For each gene, note the number of known EnsEMBL junctions/boundaries, number EST/mRNA supported junctions/boundaries, total possible junctions/boundaries
#  - Note how many of these were actually expressed (>= 1 read, >= 5 reads, >= 10 reads, >= 100 reads)

#Generate a GFF style UCSC track to display all of the junctions/boundaries actually observed to be expressed for each gene
#  - Name according to level of coverage ('ReadCount'_'SeqName'_'JunctionID')

#The term 'seq' will be used to generically refer to either exon-exon junctions or exon boundaries

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

my $seq_type = '';
my $database = '';
my $server = '';
my $user = '';
my $password = '';
my $library = '';
my $library_name = '';
my $chr_filter = '';
my $seq_database = '';
my $read_record_dir = '';
my $mapped_reads_dir = '';
my $working_dir = '';
my $min_bit_score = '';
my $min_seq_coverage = '';
my $correction_factor = '';
my $results_dir = '';
my $ucsc_dir = '';
my $ucsc_build = '';
my $web_path = '';
my $color_set = '';
my $cutoffs_file = '';
my $log_file = '';

GetOptions ('database=s'=>\$database, 'server=s'=>\$server, 'user=s'=>\$user, 'password=s'=>\$password, 'seq_type=s'=>\$seq_type,
	    'library=s'=>\$library, 'library_name=s'=>\$library_name, 'chr_filter=s'=>\$chr_filter, 'ucsc_build=s'=>\$ucsc_build,
	    'seq_database=s'=>\$seq_database, 'read_record_dir=s'=>\$read_record_dir, 'mapped_reads_dir=s'=>\$mapped_reads_dir, 'working_dir=s'=>\$working_dir,
	    'min_bit_score=f'=>\$min_bit_score, 'min_seq_coverage=f'=>\$min_seq_coverage, 'correction_factor=f'=>\$correction_factor, 
	    'results_dir=s'=>\$results_dir, 'ucsc_dir=s'=>\$ucsc_dir, 'web_path=s'=>\$web_path, 'color_set=i'=>\$color_set, 'cutoffs_file=s'=>\$cutoffs_file, 'log_file=s'=>\$log_file);

#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tThis script summarizes reads mapped to exon junctions or boundaries ('seqs'), joins this info to a sequence database file and generates UCSC tracks", RESET;
print GREEN, "\n\tSpecify whether you are processing results for exon-exon junctions or exon boundaries using:  --seq_type=junction OR --seq_type=boundary", RESET;
print GREEN, "\n\tSpecify the ALEXA database and server to query using: --database and --server", RESET;
print GREEN, "\n\tSpecify the ALEXA user and password for access using: --user and --password", RESET;
print GREEN, "\n\tSpecify the library ID for read data using: --library", RESET;
print GREEN, "\n\tSpecify the library name for read data using: --library_name", RESET;
print GREEN, "\n\tSpecify the chromosome (sub-range) to be processed using: --chr_filter", RESET;
print GREEN, "\n\tSpecify the annotated sequence database corresponding to the junctions/boundaries used for mapping: --seq_database", RESET;
print GREEN, "\n\tSpecify a directory contain COMPRESSED, junction/boundary read mapping files using: --mapped_reads_dir", RESET;
print GREEN, "\n\tSpecify a directory containing read record files for this library using:  --read_record_dir", RESET;
print GREEN, "\n\tSpecify a working directory for temp files using:  --working_dir", RESET;
print GREEN, "\n\tSpecify a minimum bit score for each hit to be accepted using: --min_bit_score", RESET;
print GREEN, "\n\tSpecify the minimum percent sequence coverage for an observed sequence to be considered expressed and written to the UCSC custom track using: --min_seq_coverage", RESET;
print GREEN, "\n\t\tThis is the percentage of bases of a junction/boundary covered to 1X depth or greater", RESET;
print GREEN, "\n\t\tBy chosing a value of say 75.0 you ensure that multiple distinct reads (different start/end points) are required for each junction/boundary", RESET;
print GREEN, "\n\t\tThis prevents singletons and cases where a junction/boundary is represented by several instances of an identical read only", RESET;
print GREEN, "\n\tSpecify a correction factor for normalized junction/boundary expression values to account for bias relative to exons, etc. using: --correction_factor (1 for no correction)", RESET;
print GREEN, "\n\tSpecify a directory to write final GENE/JUNCTION-BOUNDARY results files to using: --results_dir", RESET;
print GREEN, "\n\tSpecify the target UCSC directory for custom UCSC track files using: --ucsc_dir", RESET;
print GREEN, "\n\tSpecify the UCSC build to be used for links to UCSC tracks using: --ucsc_build (e.g. hg18)", RESET;
print GREEN, "\n\tSpecify the html web path to this directory using: --web_path", RESET;
print GREEN, "\n\tSpecify which hard coded color set to use using: --color_set (1 for LibraryA, 2 for LibraryB)", RESET;
print GREEN, "\n\tSpecify the path to a file containing expression cutoffs values using:  --cutoffs_file", RESET;
print GREEN, "\n\t\tIf these have not been calculated yet, use: --cutoffs=0", RESET;
print GREEN, "\n\tSpecify a log file using: --log_file", RESET;

print GREEN, "\n\nExample: generateExpressionValues_JunctionBoundary.pl  --seq_type=junction  --database=ALEXA_hs_53_36o  --server=jango.bcgsc.ca  --user=viewer  --password=viewer  --library=HS04391  --library_name=MIP101  --chr_filter='3:10:118167002-127874001'  --mapped_reads_dir=/projects/malachig/solexa/read_records/HS04391/Junctions_v53/  --read_record_dir=/projects/malachig/solexa/read_records/HS04391/  --seq_database=/projects/malachig/sequence_databases/hs_53_36o/exonJunctions/exonJunctions_62mers_annotated.txt.gz  --working_dir=/projects/malachig/solexa/read_records/HS04391/Junctions_v53/Summary/temp/  --min_bit_score=69.9  --min_seq_coverage=75.0  --correction_factor=1  --results_dir=/projects/malachig/solexa/read_records/HS04391/Junctions_v53/Summary/junction_results/  --ucsc_dir=/home/malachig/www/public/htdocs/solexa/HS04391/  --web_path=http://www.bcgsc.ca/people/malachig/htdocs/solexa/HS04391/  --color_set=1  --cutoffs_file=/projects/malachig/solexa/figures_and_stats/HS04391/Expression_v53/HS04391_NORM1_average_coverage_cutoffs.txt  --log_file=/projects/malachig/solexa/logs/HS04391/generateExpressionValues/Junction/generateExpressionValues_Junction_chr3_10.txt\n\n", RESET;


unless ($seq_type && $database && $server && $user && $password && $library && $library_name && $chr_filter && $seq_database && $mapped_reads_dir && $read_record_dir && $working_dir && $min_bit_score && $min_seq_coverage && $correction_factor && $results_dir && $ucsc_dir && $ucsc_build && $web_path && $color_set && ($cutoffs_file || $cutoffs_file eq '0') && $log_file){
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  exit();
}

#Scale all expression values to a constant value (values will be scaled to what they would be expected to be if the library contained 10 billion mapped bases)
#If a library has less than this, the normalized expression values will be increased and vice versa
#The actual number of mapped bases is determined by parsing the 'read record' files and counting ALL bases that have been unambiguously assigned to one of the following classes:
#'ENST_U', 'INTRON_U', 'INTERGENIC_U', 'NOVEL_JUNCTION_U', 'NOVEL_BOUNDARY_U' 
my $library_normalization_value = 10000000000;
my $library_size;

#Set priority values to determine ordering of custom tracks
my $priority_annotation;
my $priority_expression;

my $file_name_prefix;
my $region_number;
my $start_filter;
my $end_filter;
my $range;
if ($chr_filter =~ /(.*)\:(\d+)\:(\d+)\-(\d+)/){
  $chr_filter = $1;
  $region_number = $2;
  $start_filter = $3;
  $end_filter = $4;
  $file_name_prefix = "$chr_filter"."_"."$region_number";
  $range = "$start_filter-$end_filter";
  unless ($end_filter > $start_filter){
    print RED, "\nStart of range must be smaller than end ($chr_filter)\n\n", RESET;
    exit();
  }
}else{
  print RED, "\nFormat of chr_filter not understood: $chr_filter (should be of the form:  Y:1:1-9999001)\n\n", RESET;
  exit();
}
if ($chr_filter eq "MT"){$chr_filter = "M"};


my $seq_id_name;
my $seq_summary_file;
my $gene_summary_file;
my $type_name;
my $type_name_b;
my $type_name_short;
my $read_class;
my %colors;
my $pge_cutoff;

if ($seq_type =~ /junction|boundary/i){
  if ($seq_type =~ /junction/i){
    $seq_id_name = "Junction_ID";
    $type_name = "Junctions";
    $type_name_b = "Junction";
    $type_name_short = "J";
    $seq_summary_file = "$results_dir"."chr"."$file_name_prefix"."_JunctionExpression.txt";
    $gene_summary_file = "$results_dir"."chr"."$file_name_prefix"."_JunctionGeneSummary.txt";
    $read_class = "ENST_U NOVEL_JUNCTION_U";
    $pge_cutoff = 2;

    $priority_annotation = 5;
    $priority_expression = 35;

    #Define some color profiles for Junctions (Library A and B colors)
    $colors{'1'}{1} = "0,102,0";       #Dark Green
    $colors{'2'}{1} = "0,204,0";     #Light Green

  }elsif($seq_type =~ /boundary/i){
    $seq_id_name = "Boundary_ID";
    $type_name = "Boundaries";
    $type_name_b = "Boundary";
    $type_name_short = "B";
    $seq_summary_file = "$results_dir"."chr"."$file_name_prefix"."_BoundaryExpression.txt";
    $gene_summary_file = "$results_dir"."chr"."$file_name_prefix"."_BoundaryGeneSummary.txt";
    $read_class = "ENST_U NOVEL_BOUNDARY_U";
    $pge_cutoff = 20;

    $priority_annotation = 10;
    $priority_expression = 40;

    #Define some color profiles for Boundaries (Library A and B colors)
    $colors{'1'}{1} = "102,0,153";       #Dark Purple
    $colors{'2'}{1} = "153,0,153";        #Light Purple
  }
}else{
  print RED, "--seq_type specified by user is not understood.  Must be 'junction' or 'boundary'", RESET;
  exit();
}
unless ($colors{$color_set}){
  print RED, "\nColor set specified by user is not understood", RESET;
  print Dumper %colors;
  exit();
}

open (LOG, ">$log_file") || die "\nCould not open log file: $log_file\n\n";

#Check input directories
$mapped_reads_dir = &checkDir('-dir'=>$mapped_reads_dir, '-clear'=>"no");
$working_dir = &checkDir('-dir'=>$working_dir, '-clear'=>"no");
$results_dir = &checkDir('-dir'=>$results_dir, '-clear'=>"no");
$ucsc_dir = &checkDir('-dir'=>$ucsc_dir, '-clear'=>"no");

unless(-e $seq_database && ($seq_database =~ /\.gz/)){
  print RED, "\nSequence annotation database file was not present or not compressed: $seq_database\n\n", RESET;
  exit();
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


#1.) Import basic gene info for all genes - limit to those within the specified chromosome region
my $genes_ref = &getBasicGeneInfo('-chromosome'=>$chr_filter, '-range'=>"$start_filter-$end_filter");


#Open global Hashes or Berkeley DB files for reading or writing
my %seq_db;
my %passing_reads;
my %seq_counts_RAW;
my %seq_counts_NORM1;
my %seq_coverage_RAW;
my %seq_percent_cov;
my %seq_coords;

my $berk = 0;
my $rm_cmd = "rm -f $working_dir"."*chr"."$file_name_prefix"."_"."*";
if ($berk == 1){
  #Delete pre-existing berkley DB files for the current chromosome
  print YELLOW, "\nDeleting Berkley DB files and creating new\n\t($rm_cmd)\n\n", RESET;
  print LOG "\nDeleting Berkley DB files and creating new\n\t($rm_cmd)\n\n";
  system($rm_cmd);

  my $seq_db_file = "$working_dir"."chr"."$file_name_prefix"."_Junction_DB.btree";
  tie(%seq_db, 'BerkeleyDB::Btree', -Cachesize => 256000000, -Filename=> $seq_db_file, -Flags => DB_CREATE) or die "can't open file $seq_db_file: $! $BerkeleyDB::Error\n";

  my $reads_file = "$working_dir"."chr"."$file_name_prefix"."_PassingReads.btree";
  tie(%passing_reads, 'BerkeleyDB::Btree', -Cachesize => 256000000, -Filename=> $reads_file , -Flags => DB_CREATE) or die "can't open file $reads_file: $! $BerkeleyDB::Error\n";

  my $seq_counts_RAW_file = "$working_dir"."chr"."$file_name_prefix"."_SeqCounts_RAW.btree";
  tie(%seq_counts_RAW, 'BerkeleyDB::Btree', -Cachesize => 256000000, -Filename=> $seq_counts_RAW_file, -Flags => DB_CREATE) or die "can't open file $seq_counts_RAW_file: $! $BerkeleyDB::Error\n";

  my $seq_coverage_RAW_file = "$working_dir"."chr"."$file_name_prefix"."_SeqCoverage_RAW.btree";
  tie(%seq_coverage_RAW, 'BerkeleyDB::Btree', -Cachesize => 256000000, -Filename=> $seq_coverage_RAW_file, -Flags => DB_CREATE) or die "can't open file $seq_coverage_RAW_file: $! $BerkeleyDB::Error\n";

  my $seq_counts_NORM1_file = "$working_dir"."chr"."$file_name_prefix"."_SeqCounts_NORM1.btree";
  tie(%seq_counts_NORM1, 'BerkeleyDB::Btree', -Cachesize => 256000000, -Filename=> $seq_counts_NORM1_file, -Flags => DB_CREATE) or die "can't open file $seq_counts_NORM1_file: $! $BerkeleyDB::Error\n";

  my $seq_percent_cov_file = "$working_dir"."chr"."$file_name_prefix"."_SeqPercentCov.btree";
  tie(%seq_percent_cov, 'BerkeleyDB::Btree', -Cachesize => 256000000, -Filename=> $seq_percent_cov_file, -Flags => DB_CREATE) or die "can't open file $seq_percent_cov_file: $! $BerkeleyDB::Error\n";

  my $seq_coords_file = "$working_dir"."chr"."$file_name_prefix"."SeqCoords.btree";
  tie(%seq_coords, 'BerkeleyDB::Btree', -Cachesize => 256000000, -Filename=> $seq_coords_file, -Flags => DB_CREATE) or die "can't open file $seq_coords_file: $! $BerkeleyDB::Error\n";
}


#2.) Import junction/boundary database and store as a berkley DB
#    - At this time, for each gene, note the number of known junctions/boundaries, number EST/mRNA supported junctions/boundaries, total possible junctions/boundaries
#    - Only grab the info actually needed for now...
#    - Only store those junctions/boundaries within the specified chromosome region
&importSeqDatabase('-seq_db_file'=>$seq_database, '-chromosome'=>$chr_filter, '-range'=>"$start_filter-$end_filter");


#3-A.) Get a list of read record files corresponding to Junction/Boundary mapped reads
my $mapped_reads_files_ref = &getReadFiles('-input_dir'=>$mapped_reads_dir);

#3-B.) Build a hash of all reads with hits to Junction/Boundary mapped reads
#      - Only store those reads that correspond to one of the Junction/Boundary IDs within the specified chromsome region
&getMappedReads('-reads_ref'=>\%passing_reads, '-reads_files'=>$mapped_reads_files_ref, '-min_bit_score'=>$min_bit_score);

#3-C.) Go through the read records files and remove reads from those stored if they have not been assigned to the correct class(es)
$library_size = &importLibrarySize ('-read_records_dir'=>$read_record_dir);


#4.) Parse mapped reads records and store junction/boundary level info as Berkley DBs
#    - Create Berkeley DB files to store the count and coverage for all junctions/boundaries with 1 or more successfully mapped reads
#    - Only allow 'Top_Hit' reads that meet the min_bit_score cutoff
my $grand_total_mapped_reads = 0;  #Count all the quality unambiguous mapped reads from ANY chromosome
&parseSeqHitRecords('-mapped_reads_files'=>$mapped_reads_files_ref, '-min_bit_score'=>$min_bit_score, '-working_dir'=>$working_dir);

#5.) Join expression values to the junction/boundary database file
my %chr_seqs;
&summarizeSeqResults('-seq_db_file'=>$seq_database, '-seq_summary_file'=>$seq_summary_file,
		     '-library_size'=>$library_size, '-library_normalization_value'=>$library_normalization_value, '-correction_factor'=>$correction_factor,
		     '-min_seq_coverage'=>$min_seq_coverage);

#6.) Summarize gene results and create UCSC tracks
&summarizeGeneResults('-gene_summary_file'=>$gene_summary_file, '-library'=>$library, '-min_seq_coverage'=>$min_seq_coverage,
		      '-ucsc_dir'=>$ucsc_dir, '-web_path'=>$web_path, '-color_set'=>$color_set);
#If neccessary, clean berkeley DB files
if ($berk == 1){
  untie(%seq_coords);
  untie(%seq_counts_RAW);
  untie(%seq_counts_NORM1);
  untie(%seq_percent_cov);
  untie(%seq_db);
  untie(%seq_coverage_RAW);
  untie(%passing_reads);

  #Delete temp files
  print BLUE, "\n\nDeleting temporary files\n\n", RESET;
  system($rm_cmd);
}

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
#getBasicGeneInfo                                                                                                                  #
####################################################################################################################################
sub getBasicGeneInfo{
  my %args = @_;
  my $target_chr = $args{'-chromosome'};
  my $range = $args{'-range'};

  #Open database connection
  my $alexa_dbh = &connectDB('-database'=>$database, '-server'=>$server, '-user'=>$user, '-password'=>$password);

  my @gene_ids = @{&getAllGenes ('-dbh'=>$alexa_dbh, '-gene_type'=>'All', '-chromosome'=>$target_chr, '-range'=>$range)};
  my $genes_ref = &getGeneInfo ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids, '-sequence'=>"no");


  #Close database connection
  $alexa_dbh->disconnect();

  foreach my $gene_id (keys %{$genes_ref}){
    my $chromosome = $genes_ref->{$gene_id}->{chromosome};
    if ($chromosome eq "MT"){$chromosome = "M";}

    $genes_ref->{$gene_id}->{total_seqs} = 0;
    $genes_ref->{$gene_id}->{ensembl_supported_seqs} = 0;
    $genes_ref->{$gene_id}->{mrna_est_supported_seqs} = 0;
    $genes_ref->{$gene_id}->{observed_known_seqs} = 0;
    $genes_ref->{$gene_id}->{observed_supported_seqs} = 0;
    $genes_ref->{$gene_id}->{observed_novel_seqs} = 0;
  }
  my $gene_count = keys %{$genes_ref};
  print BLUE, "\n\nFound $gene_count genes within chr$chr_filter:$range\n\n", RESET;
  print LOG "\n\nFound $gene_count genes within chr$chr_filter:$range\n\n";

  return($genes_ref);
}


#############################################################################################################################################
#importSeqDatabase();                                                                                                                       #
#############################################################################################################################################
sub importSeqDatabase{
  my %args = @_;
  my $seq_db_file = $args{'-seq_db_file'};
  my $target_chr = $args{'-chromosome'};
  my $range = $args{'-range'};

  $target_chr = "chr"."$target_chr";

  my $start_filter;
  my $end_filter;
  if ($range =~ /(\d+)\-(\d+)/){
    $start_filter = $1;
    $end_filter = $2;
  }else{
    print RED, "\nProblem with range - getMappedReads()\n\n", RESET;
    exit();
  } 

  print BLUE, "\n\nImporting seqs from the sequence database: $seq_db_file - limiting to $target_chr:$range\n\n", RESET;
  print LOG "\n\nImporting seqs from the sequence database: $seq_db_file - limiting to $target_chr:$range\n\n";

  open (SEQ, " zcat $seq_db_file |") || die "\nCould not open sequence db file: $seq_db_file\n\n";

  my $header = 1;
  my %columns;
  my $counter = 0;

  while(<SEQ>){
    $counter++;
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

    if ($counter == 10000){
      $| = 1; print BLUE, ".", RESET;  $| = 0;
      print LOG ".";
      $counter = 0;
    }

    #Store neccessary values to a Berkley DB file:
    #Seq_ID, Gene_ID, Supporting_EnsEMBL_Count, Supporting_mRNA_Count, Supporting_EST_Count

    my $seq_id = $line[$columns{$seq_id_name}{pos}];
    my $gene_id = $line[$columns{'Gene_ID'}{pos}];
    my $ensembl_trans_count = $line[$columns{'Supporting_EnsEMBL_Count'}{pos}];
    my $mrna_count = $line[$columns{'Supporting_mRNA_Count'}{pos}];
    my $est_count = $line[$columns{'Supporting_EST_Count'}{pos}];
    my $chromosome = $line[$columns{'Chromosome'}{pos}];

    my @coords_sort;
    if ($seq_type =~ /junction/i){
      my $unit1_start_chr = $line[$columns{'Unit1_start_chr'}{pos}];
      my $unit1_end_chr = $line[$columns{'Unit1_end_chr'}{pos}];
      my $unit2_start_chr = $line[$columns{'Unit2_start_chr'}{pos}];
      my $unit2_end_chr = $line[$columns{'Unit2_end_chr'}{pos}];
      my @coords = ($unit1_start_chr, $unit1_end_chr, $unit2_start_chr, $unit2_end_chr);
      @coords_sort = sort {$a <=> $b} @coords;
    }else{
      my $unit1_start_chr = $line[$columns{'Unit1_start_chr'}{pos}];
      my $unit1_end_chr = $line[$columns{'Unit1_end_chr'}{pos}];
      my @coords = ($unit1_start_chr, $unit1_end_chr);
      @coords_sort = sort {$a <=> $b} @coords;
    }
    my $lower = $coords_sort[0];
    my $upper = $coords_sort[scalar(@coords_sort)-1];

    #Unless the current junction/boundary is contained within the current target chromosome region, skip it 
    if (($lower >= $start_filter) && ($lower <= $end_filter) && ($upper >= $start_filter) && ($upper <= $end_filter) && ($chromosome eq $chr_filter)){
      #This seq is within the target range 
      #print YELLOW, "\nWITHIN: $seq_id\t$gene_id\t$lower\t$upper", RESET;
    }else{
      #print MAGENTA, "\nOUTSIDE: $seq_id\t$gene_id\t$lower\t$upper", RESET;
      next();
    }

    unless($genes_ref->{$gene_id}){
      print RED, "\nGene ID: $gene_id from seq DB file is not defined in the Genes object!!\n\n", RESET;
      exit();
    }

    #Store this SEQ in the list to be processed
    $seq_db{$seq_id} = "$gene_id";
 
    #Add counts for each junction/boundary to its corresponding gene
    #Add this junction/boundary to the count for its gene
    $genes_ref->{$gene_id}->{total_seqs}++;

    #Is this junction/boundary supported by one or more EnsEMBL transcripts?
    if ($ensembl_trans_count >= 1){
      $genes_ref->{$gene_id}->{ensembl_supported_seqs}++;
    }

    if ($mrna_count eq "na"){$mrna_count = 0;}
    if ($est_count eq "na"){$est_count = 0;}

    #Is this junction/boundary supported by one or more mRNAs or ESTs?
    if ($mrna_count >= 1 || $est_count >= 1){
      $genes_ref->{$gene_id}->{mrna_est_supported_seqs}++;
    }
  }

  close(SEQ);

  my $seq_count = keys %seq_db;

  print BLUE, "\n\nA total of $seq_count $seq_type sequences were imported from the sequence database for $target_chr:$range\n\n", RESET;
  print LOG "\n\nA total of $seq_count $seq_type sequences were imported from the sequence database for $target_chr:$range\n\n";

  return();
}


####################################################################################################################################
#getReadFiles                                                                                                                      #
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

  return(\%read_files);
}




####################################################################################################################################
#getMappedReads                                                                                                                    #
####################################################################################################################################
sub getMappedReads{
  my %args = @_;
  my $reads_ref = $args{'-reads_ref'};
  my $read_files_ref = $args{'-reads_files'};
  my $min_bit_score = $args{'-min_bit_score'};

  print BLUE, "\n\nImporting list of all mapped reads corresponding to chr$chr_filter:$range", RESET;
  print LOG "\n\nImporting list of all mapped reads corresponding to chr$chr_filter:$range";

  foreach my $fc (sort {$read_files_ref->{$a}->{flowcell_lane} cmp $read_files_ref->{$b}->{flowcell_lane}} keys %{$read_files_ref}){
    my $file = $read_files_ref->{$fc}{file_path};
    my $header = 1;
    my $counter = 0;

    my %columns;

    print BLUE, "\n\tProcessing: $file\n\t", RESET;
    print LOG "\n\tProcessing: $file\n\t";

    open (MAPPED, "zcat $file |") || die "\nCould not open file: $file\n\n";
    while(<MAPPED>){
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
      my $read_id = $line[$columns{'Read_ID'}{pos}];
      my $hit_type = $line[$columns{'HitType'}{pos}];
      my $bit_score = $line[$columns{'BitScore'}{pos}];
      my $seq_id = $line[$columns{$seq_id_name}{pos}];

      #Make sure this is a Top_Hit of sufficient quality AND that it corresponds to one of the SEQs in the target region
      unless ($hit_type eq "Top_Hit" && $bit_score >= $min_bit_score && $seq_db{$seq_id}){
	next();
      }
      $reads_ref->{$read_id} = 1;
    }
    close(MAPPED);
  }

  my $reads_count = keys %{$reads_ref};
  print BLUE, "\n\nStored $reads_count passing mapped reads", RESET;
  print LOG "\n\nStored $reads_count passing mapped reads";

  return();
}



#############################################################################################################################################
#parseSeqHitRecords()                                                                                                                       #
#############################################################################################################################################
sub parseSeqHitRecords{
  my %args = @_;
  my $files_ref = $args{'-mapped_reads_files'};
  my $min_bit_score = $args{'-min_bit_score'};
  my $working_dir = $args{'-working_dir'};

  my $total_read_hits = 0;
  my $passing_read_hits = 0;

  print BLUE, "\n\nProcessing $seq_type sequence read hits from the list of input files\n", RESET;
  print LOG "\n\nProcessing $seq_type sequence read hits from the list of input files\n";

  #Go through each mapped reads file and parse the count and coverage data for reads that pass quality criteria
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

      my $read_id = $line[$columns{'Read_ID'}{pos}];
      my $seq_id = $line[$columns{$seq_id_name}{pos}];
      my $hit_type = $line[$columns{'HitType'}{pos}];
      my $bit_score = $line[$columns{'BitScore'}{pos}];
      my $start = $line[$columns{'Start'}{pos}];
      my $end = $line[$columns{'End'}{pos}];

      unless ($hit_type eq "Top_Hit"){
        next();
      }

      $total_read_hits++;
      unless ($passing_reads{$read_id}){
	next();
      }
      $passing_read_hits++;
      $grand_total_mapped_reads++;

      #Add this read to the read counts record for this junction/boundary
      if ($seq_counts_RAW{$seq_id}){
	$seq_counts_RAW{$seq_id}++;
      }else{
	$seq_counts_RAW{$seq_id} = 1;
      }

      #Add the coverage for this read to the coverage record for this junction/boundary

      #Go through each position of this junction/boundary hit by this read (1-62 for a 62 mer)
      for (my $i = $start; $i <= $end; $i++){
	my $id_pos = "$seq_id"."_"."$i";

	if ($seq_coverage_RAW{$id_pos}){
	  $seq_coverage_RAW{$id_pos}++;
	}else{
	  $seq_coverage_RAW{$id_pos} = 1;
	}
      }
    }
    close(READS);
  }

  my $seqs_with_reads = keys %seq_counts_RAW;
  my $seq_bases_covered = keys %seq_coverage_RAW;

  print BLUE, "\n\nSUMMARY for chr$chr_filter:$range\n", RESET;
  print BLUE, "\n\nA total of $total_read_hits read hits were imported and $passing_read_hits of these reads passed ('TOP_HIT' and bit_score >= $min_bit_score)\n", RESET;
  print BLUE, "\n\nA total of $seqs_with_reads $seq_type sequences were hit by one more reads\n", RESET;
  print BLUE, "\n\nA total of $seq_bases_covered $seq_type bases were covered to 1X depth or greater\n", RESET;
  print BLUE, "\n\nGRAND TOTAL MAPPED READS = $grand_total_mapped_reads\n\n", RESET;

  print LOG "\n\nSUMMARY for chr$chr_filter:$range\n";
  print LOG "\n\nA total of $total_read_hits read hits were imported and $passing_read_hits of these reads passed ('TOP_HIT' and bit_score >= $min_bit_score)\n";
  print LOG "\n\nA total of $seqs_with_reads $seq_type sequences were hit by one more reads\n";
  print LOG "\n\nA total of $seq_bases_covered $seq_type bases were covered to 1X depth or greater\n";
  print LOG "\n\nGRAND TOTAL MAPPED READS = $grand_total_mapped_reads\n\n";

  return();
}


#############################################################################################################################################
#summarizeSeqResults()
#############################################################################################################################################
sub summarizeSeqResults{
  my %args = @_;
  my $infile = $args{'-seq_db_file'};
  my $outfile = $args{'-seq_summary_file'};
  my $library_size = $args{'-library_size'};
  my $library_normalization_value = $args{'-library_normalization_value'};
  my $junction_boundary_correction = $args{'-correction_factor'};
  my $min_seq_coverage = $args{'-min_seq_coverage'};

  my $library_correction_ratio = 1;
  unless ($library_size == 0){
    $library_correction_ratio = $library_normalization_value/$library_size;
  }

  print BLUE, "\n\nCreating a seq summary file: $outfile\n\n", RESET;
  print LOG "\n\nCreating a seq summary file: $outfile\n\n";

  open (SEQ, "zcat $infile |") || die "\nCould not open sequence db file: $infile\n\n";
  open (SEQ_SUM, ">$outfile") || die "\nCould not open sequence summary file: $outfile\n\n";

  my $header = 1;
  my %columns;
  my $counter = 0;

  while(<SEQ>){
    $counter++;
    chomp($_);
    my @line = split ("\t", $_);

    if ($header == 1){
      my $col_count = 0;
      foreach my $name (@line){
	$columns{$name}{pos} = $col_count;
	$col_count++;
      }
      $header = 0;

      #Write the header to a seperate file, unless it is already present
      print SEQ_SUM "$_\tRead_Count\tCumulative_Coverage\tAverage_Coverage_RAW\tAverage_Coverage_NORM1\tBases_Covered_1x\tPercent_Coverage_1x\tPercent_Coverage_5x\tPercent_Coverage_10x\tPercent_Coverage_100x\tExpressed\tPercent_Gene_Expression\n";
      next();
    }

    if ($counter == 10000){
      $| = 1; print BLUE, ".", RESET;  $| = 0;
      print LOG ".";
      $counter = 0;
    }

    my $seq_id = $line[$columns{$seq_id_name}{pos}];
    my $gene_id = $line[$columns{'Gene_ID'}{pos}];
    my $ensembl_trans_count = $line[$columns{'Supporting_EnsEMBL_Count'}{pos}];
    my $mrna_count = $line[$columns{'Supporting_mRNA_Count'}{pos}];
    my $est_count = $line[$columns{'Supporting_EST_Count'}{pos}];

    #Once again, skip this record unless it corresponds to to a SEQ within the target chromosome region
    unless ($seq_db{$seq_id}){
      next();
    }

    my $strand = $line[$columns{'Strand'}{pos}];
    my ($base_count, $unit1_start_chr, $unit1_end_chr, $unit2_start_chr, $unit2_end_chr, $seq_name, $exons_skipped);

    if ($seq_type =~ /junction/i){
      $base_count = $line[$columns{'Base_Count'}{pos}];
      $unit1_start_chr = $line[$columns{'Unit1_start_chr'}{pos}];
      $unit1_end_chr = $line[$columns{'Unit1_end_chr'}{pos}];
      $unit2_start_chr = $line[$columns{'Unit2_start_chr'}{pos}];
      $unit2_end_chr = $line[$columns{'Unit2_end_chr'}{pos}];
      $seq_name = $line[$columns{'Seq_Name'}{pos}];
      $exons_skipped = $line[$columns{'Exons_Skipped'}{pos}];
    }elsif ($seq_type =~ /boundary/i){
      $base_count = $line[$columns{'Base_Count'}{pos}];
      $unit1_start_chr = $line[$columns{'Unit1_start_chr'}{pos}];
      $unit1_end_chr = $line[$columns{'Unit1_end_chr'}{pos}];
      $unit2_start_chr = "na";
      $unit2_end_chr = "na";
      $seq_name = $line[$columns{'Seq_Name'}{pos}];
      $exons_skipped = "na";
    }

    if ($mrna_count eq "na"){$mrna_count = 0;}
    if ($est_count eq "na"){$est_count = 0;}

    #Get the read count for this junction/boundary (set to 0 if it is undefined)
    #Get the coverage values for this junction/boundary
    #Possible values are from 1 to $base_count (assuming the correct value was specified by the user)
    my $cumulative_coverage_raw = 0;
    my $bases_covered_1x = 0;
    my $bases_covered_5x = 0;
    my $bases_covered_10x = 0;
    my $bases_covered_100x = 0;
    my $read_count = 0;

    if ($seq_counts_RAW{$seq_id}){
      $read_count = $seq_counts_RAW{$seq_id};

      #If there was no read count for this junction/boundary then there would be no coverage
      for (my $i = 1; $i <= $base_count; $i++){
	my $id_pos = "$seq_id"."_"."$i";

	if ($seq_coverage_RAW{$id_pos}){

	  $cumulative_coverage_raw += $seq_coverage_RAW{$id_pos};
	  $bases_covered_1x++;

	  if ($seq_coverage_RAW{$id_pos} >= 5){$bases_covered_5x++;}
	  if ($seq_coverage_RAW{$id_pos} >= 10){$bases_covered_10x++;}
	  if ($seq_coverage_RAW{$id_pos} >= 100){$bases_covered_100x++;}
	}
      }
    }

    my $average_coverage_raw = sprintf("%.10f", ($cumulative_coverage_raw/$base_count));
    my $percent_coverage_1x = sprintf("%.2f", (($bases_covered_1x/$base_count)*100));
    my $percent_coverage_5x = sprintf("%.2f", (($bases_covered_5x/$base_count)*100));
    my $percent_coverage_10x = sprintf("%.2f", (($bases_covered_10x/$base_count)*100));
    my $percent_coverage_100x = sprintf("%.2f", (($bases_covered_100x/$base_count)*100));

    #Apply the library correction factor as well as the junction/boundary correction factor
    my $average_coverage_norm1 = (($average_coverage_raw*$library_correction_ratio)*$junction_boundary_correction);
    my $average_coverage_norm1_f = sprintf("%.10f", $average_coverage_norm1);
    my $average_coverage_norm1_ff = sprintf("%.1f", $average_coverage_norm1);

    #Count this junction/boundary towards those OBSERVED for its gene
    #Count known-expressed-junctions/boundaries AND novel-sequence-supported-expressed-junctions/boundaries AND novel-expressed-junctions/boundaries
    my $expressed = 0;
    my $percent_gene_expression = 0;
    if ($seq_counts_RAW{$seq_id}){

      $seq_counts_NORM1{$seq_id} = $average_coverage_norm1_ff;
      $seq_percent_cov{$seq_id} = $percent_coverage_1x;

      #Determine whether this element should be considered as expressed above background
      my $cutoff_test = 1;
      if ($cutoffs_file){
        my @result = @{&testExpression('-cutoffs_ref'=>$gene_cutoffs_ref, '-gene_id'=>$gene_id, '-norm_expression_value'=>$average_coverage_norm1, '-raw_expression_value'=>$average_coverage_raw, '-percent_gene_expression_cutoff'=>$pge_cutoff)};
        $cutoff_test = $result[0];
        $percent_gene_expression = $result[1];
      }

      #Unless the junction/boundary had a read it can not be counted towards a gene
      #Also do not count it towards the gene count unless it meets the coverage cutoff
      if ($percent_coverage_1x >= $min_seq_coverage && $cutoff_test == 1){
        $expressed = 1;
	if ($ensembl_trans_count >= 1){
	  $genes_ref->{$gene_id}->{observed_known_seqs}++;
	}elsif($mrna_count >= 1 || $est_count >= 1){
	  $genes_ref->{$gene_id}->{observed_supported_seqs}++;
	}else{
	  $genes_ref->{$gene_id}->{observed_novel_seqs}++;
	}

        #If this junction/boundary was expressed above the cutoff, store it in a hash keyed on chromosome so that later it can be printed to a UCSC track
        #Also store the coordinates and strand for observed exon junctions/boundaries
        my $string = "$strand\t$unit1_start_chr\t$unit1_end_chr\t$unit2_start_chr\t$unit2_end_chr\t$seq_name\t$exons_skipped\t$expressed";
        $seq_coords{$seq_id} = $string;
        my $chromosome = "$genes_ref->{$gene_id}->{chromosome}";

        if ($chr_seqs{$chromosome}){
	  my $temp_ref = $chr_seqs{$chromosome};
	  $temp_ref->{$seq_id} = $seq_counts_RAW{$seq_id};
        }else{
	  my %temp;
	  $temp{$seq_id} = $seq_counts_RAW{$seq_id};
	  $chr_seqs{$chromosome} = \%temp;
        }

      }
    }

    #Print the junction/boundary line with expression data appended
    print SEQ_SUM "$_\t$read_count\t$cumulative_coverage_raw\t$average_coverage_raw\t$average_coverage_norm1_f\t$bases_covered_1x\t$percent_coverage_1x\t$percent_coverage_5x\t$percent_coverage_10x\t$percent_coverage_100x\t$expressed\t$percent_gene_expression\n";

  }

  close (SEQ);
  close (SEQ_SUM);

  return();
}


#############################################################################################################################################
#summarizeGeneResults()
#############################################################################################################################################
sub summarizeGeneResults{
  my %args = @_;
  my $gene_summary_file = $args{'-gene_summary_file'};
  my $library = $args{'-library'};
  my $min_seq_coverage = $args{'-min_seq_coverage'};
  my $ucsc_dir = $args{'-ucsc_dir'};
  my $web_path = $args{'-web_path'};
  my $color_set = $args{'-color_set'};

  print BLUE, "\n\nCreating a gene summary file: $gene_summary_file\n\n", RESET;
  print LOG "\n\nCreating a gene summary file: $gene_summary_file\n\n";

  my $col1 = "Total_"."$type_name";
  my $col2 = "EnsEMBL_Known_"."$type_name";
  my $col3 = "mRNA-EST_Supported_"."$type_name";
  my $col4 = "Observed_Known_"."$type_name";
  my $col5 = "Observed_Supported_"."$type_name";
  my $col6 = "Observed_Novel_"."$type_name";


  open (GENE, ">$gene_summary_file") || die "\nCould not open gene summary file: $gene_summary_file\n\n";
  print GENE "ALEXA_ID\tEnsEMBL_Gene_ID\tGene_Name\t$col1\t$col2\t$col3\t$col4\t$col5\t$col6\tLink\n";
 
  foreach my $gene_id (keys %{$genes_ref}){

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

    #Link that retains pre-existing custom tracks
    #Clean link
    my $link = "http://genome.ucsc.edu/cgi-bin/hgTracks?db="."$ucsc_build"."&position=$chromosome:$display_start-$display_end&hgt.customText=$web_path"."chr$file_name_prefix"."_"."$type_name_b".".txt.gz&ctfile_"."$ucsc_build"."=";
 
    print GENE "$gene_id\t$genes_ref->{$gene_id}->{ensembl_g_id}\t$genes_ref->{$gene_id}->{gene_name}\t$genes_ref->{$gene_id}->{total_seqs}\t$genes_ref->{$gene_id}->{ensembl_supported_seqs}\t$genes_ref->{$gene_id}->{mrna_est_supported_seqs}\t$genes_ref->{$gene_id}->{observed_known_seqs}\t$genes_ref->{$gene_id}->{observed_supported_seqs}\t$genes_ref->{$gene_id}->{observed_novel_seqs}\t$link\n";

  }
  close(GENE);

  my $database_abr = $database;
  if ($database =~ /ALEXA_(\w+)/){
    $database_abr = $1;
  }

  print BLUE, "\n\nCreating a UCSC tracks files - one per chromosome", RESET;
  print LOG "\n\nCreating a UCSC tracks files - one per chromosome";

  #EXPRESSED REGION TRACKS
  #Create UCSC tracks
  foreach my $chr (sort keys %chr_seqs){

    my $seqs_ref = $chr_seqs{$chr};
    my $ucsc_file = "$ucsc_dir"."chr$file_name_prefix"."_"."$type_name_b".".txt";

    print YELLOW, "\n\twriting: $ucsc_file", RESET;
    print LOG "\n\twriting: $ucsc_file";

    open (UCSC, ">$ucsc_file") || die "\nCould not open UCSC output file: $ucsc_file\n\n";

    #Print track def lines

    #Browser line
    print UCSC "#Browser line";
    print UCSC "\nbrowser hide all";
    print UCSC "\nbrowser full knownGene";
    print UCSC "\nbrowser pack multiz28way";

    #Track line for observed exon junctions/boundaries
    print UCSC "\n\n#Observed exon $seq_type track";
    my $track_name = "$library"."_EXP_"."$type_name_short";
    my $track_description = "\"Expressed Exon $type_name (>= $min_seq_coverage% coverage and above intronic background cutoff). ($library_name - $database_abr)\"";
    $priority_expression++;

    my $color = $colors{$color_set}{1};
    print UCSC "\ntrack name=$track_name description=$track_description color=$color useScore=0 visibility=3 priority=$priority_expression";
    print UCSC "\n\n#Begin DATA\n";

    foreach my $seq_id (sort keys %{$seqs_ref}){
      #Print data lines

      my $string = $seq_coords{$seq_id};
      my @data = split("\t", $string);

      #Make sure data string has correct number of values
      unless (scalar(@data) == 8){
	print RED, "\nData object does not have expected number of values\n\n", RESET;
	exit();
      }

      my $strand = $data[0];
      my $unit1_start = $data[1];
      my $unit1_end = $data[2];
      my $unit2_start = $data[3];
      my $unit2_end = $data[4];
      my $seq_name = $data[5];
      my $exons_skipped = $data[6];
      my $expressed = $data[7];

      #See if this junction/boundary passed the percent coverage value to allow display in the UCSC tracks
      #e.g. only display junctions/boundaries that were observed at 75% (or 47 of 62 bases covered) coverage or greater
      unless ($expressed == 1){
	next();
      }

      if ($strand == 1){
	$strand = "+";
      }else{
	$strand = "-";
      }

      #Track records for Junctions
      if ($seq_type =~ /junction/i){
   
        #Name of junction: seq_name _ read_count_NORM1 _ junction_id
        my $junction_name = "$seq_name"."_S"."$exons_skipped"."_"."$seq_counts_NORM1{$seq_id}"."_"."$seq_id";

        #GFF Format: seqname  source  feature  start  end  score  strand  frame  group
        my $record_1 = "chr$chr\tEnsEMBL49\tE-E_Junction\t$unit1_start\t$unit1_end\t.\t$strand\t.\t$junction_name";
        my $record_2 = "chr$chr\tEnsEMBL49\tE-E_Junction\t$unit2_start\t$unit2_end\t.\t$strand\t.\t$junction_name";

        print UCSC "\n$record_1\n$record_2";
      }

      #Track records for Boundaries
      if ($seq_type =~ /boundary/i){
        #Name of boundary: seq_name _ read_count_NORM1 _ boundary_id
        my $boundary_name = "$seq_name"."_"."$seq_counts_NORM1{$seq_id}"."_"."$seq_id";

        #GFF Format: seqname  source  feature  start  end  score  strand  frame  group
        my $record_1 = "chr$chr\tEnsEMBL49\tE_Boundary\t$unit1_start\t$unit1_end\t.\t$strand\t.\t$boundary_name";

        print UCSC "\n$record_1";
       }
    }

    close(UCSC);

    system ("gzip -f $ucsc_file");
    print YELLOW, "\n\tCompressing: $ucsc_file\n", RESET;
    print LOG "\n\tCompressing: $ucsc_file\n";

  }

  return();
}
