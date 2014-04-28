#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#The purpose of this script is to examine coverage of intronic and intergenic bases as a function of Illumina read depth

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;
use Benchmark;
use Tie::File;
use List::Util 'shuffle';

#Load the ALEXA libraries
BEGIN {
  use Cwd 'abs_path';
  if (abs_path($0) =~ /(.*)\/.*\/.*\.pl/){
    push (@INC, $1);
  }
}
use utilities::ALEXA_DB qw(:all);
use utilities::utility qw(:all);

my $seq_database = '';
my $read_records_dir = '';
my $min_bit_score = '';          #Minimum bit score of BLAST hit for each read to be allowed for read-to-gene mappings
my $block_size = '';             #Number of reads before coverage is summarized
my $expressed_base_value = '';   #Coverage for an individual base required for it to be considered expressed/detected
my $outfile = '';
my $type = '';

GetOptions ('seq_database=s'=>\$seq_database, 'read_records_dir=s'=>\$read_records_dir, 'min_bit_score=f'=>\$min_bit_score,
	    'block_size=i'=>\$block_size, 'expressed_base_value=i'=>\$expressed_base_value,
	    'outfile=s'=>\$outfile, 'type=s'=>\$type);

#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the ALEXA database and server to query using: --database and --server", RESET;
print GREEN, "\n\tSpecify the ALEXA user and password for access using: --user and --password", RESET;
print GREEN, "\n\tspecify a path to the mapped reads files to be summarized using: --read_records_dir\n", RESET;
print GREEN, "\n\tSpecify the minimum bit score to a transcript for each read to be considered for the summary using: --min_bit_score", RESET;
print GREEN, "\n\tSpecify the number of reads to process each time before summarizing coverage using: --block_size", RESET;
print GREEN, "\n\tSpecify the number of reads required for a gene to be considered detected using: --expressed_gene_value", RESET;
print GREEN, "\n\tSpecify the X coverage level required for an exonic base to be considered detected using:  --expressed_base_value", RESET;
print GREEN, "\n\tSpecify the outfile summarizing coverage after each read block using: --outfile", RESET;
print GREEN, "\n\tSpecify the sequence type using:  --type=intronic or --type=intergenic", RESET;

print GREEN, "\n\nExample: coverageVsDepth_Gene-ExonicBase.pl  --seq_database=/projects/malachig/sequence_databases/hs_49_36k/ensembl_introns_hs_49_36k/introns_annotated.txt.gz  --read_records_dir=/projects/malachig/solexa/read_records/HS04391/Introns_v49/  --min_bit_score=48.1  --block_size=100000  --expressed_base_value=1  --outfile=/projects/malachig/solexa/read_records/HS04391/Summary/coverageVsDepth/HS04391_Lanes1-23_Intron_Coverage_10xCoverage_v49.txt  --type=intronic\n\n", RESET;

unless ($seq_database && $read_records_dir && $min_bit_score && ($block_size =~ /\d+/) && ($expressed_base_value =~ /\d+/) && $outfile && ($type =~ /intronic|intergenic/)){
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  exit();
}
unless(-e $seq_database){
  print RED, "\nSeq database: $seq_database not found\n\n", RESET;
  exit();
}

#Get list of input files containing data from the specified directory
$read_records_dir = &checkDir('-dir'=>$read_records_dir, '-clear'=>"no");
my %files;
&getDataFiles('-input_dir'=>$read_records_dir);

#Get the process ID for this script
my $pid = $$;

open (OUT, ">$outfile") || die "\nCould not open output file: $outfile\n\n";

#0.) First get all the neccessary intron/intergenic info required to perform the analysis
my $seq_ref;
my $last_i_count = 0;
my $last_bases_count = 0;

#import sequence database
&importSeqDatabase('-seq_database'=>$seq_database);

#1.) Now initialize arrays to store coverage info
my $chr_i_coverage_ref;
my $total_i_bases = 0;

print OUT "#User specified parameters:\n";
print OUT "#min_bit_score = $min_bit_score\n";
print OUT "#expressed_base_value = $expressed_base_value\n";
print OUT "#read_records_dir = $read_records_dir\n";
print OUT "#database = $seq_database\n";

$| = 1; print BLUE, "\n\n2.) Determining the non-redundant $type base coverage for ALL chromosomes\n", RESET; $| = 0;
$chr_i_coverage_ref = &getCoverageObject();

print OUT "#MappedReads\t$type"."Coverage\n";

$| = 1; print BLUE, "\n\n3-a.) Now parsing through read records files and reporting coverage after each block of $block_size reads is added\n", RESET; $| = 0;

my $total_reads_parsed = 0;
my $mapped_reads = 0;
my $grand_mapped_reads = 0;

#Report the processing time for each block of $block_size reads processed
my $t0 = new Benchmark;

#Go through each read record file and parse.
my $temp_file = "$read_records_dir"."coverage_"."$type"."_temp.txt";
my $file_io_fails = 0;
my $fc = 0;
foreach my $file_count (sort {$files{$a}->{random} <=> $files{$b}->{random}} keys %files){
  $fc++;
  print GREEN, "\nProcessing file $fc: $files{$file_count}{file_name}", RESET;

  #Access file 'in-place' as an array object using 'Tie::File'
  #Create a list of index values from line 2 to the end of the file (skip header line) [1..$line_count]
  #Then randomize the order of these index values and access the lines of the input file in this random order
  #Since 'Tie' does not seem to work on compressed file handles, the file will first have to be decompressed to a temp file
  my $cmd = "zcat $files{$file_count}{file_path} > $temp_file";
  print YELLOW, "\n\tExecuting: $cmd", RESET;
  system($cmd);
  use Fcntl 'O_RDONLY';
  my @file_array;
  tie @file_array, 'Tie::File', "$temp_file", mode => O_RDONLY;
  my $line_count = scalar(@file_array)-1;
  my $first_record = 1;
  my @record_list = (1..$line_count);
  my $record_count = scalar(@record_list)-1;
  print BLUE, "\n\tFound $line_count lines ($record_count records) in the array version of $temp_file\n", RESET;
  
  my %columns = %{$files{$file_count}{columns}};

  foreach my $i (shuffle(@record_list)){
    my $current_line = $file_array[$i];
    my @line = split("\t", $current_line);

    #Watch for sporadic failure to retrieve a line
    unless (scalar(@line) == 22){
      print MAGENTA, "\n\tFile I/O error. Retrieved:\n\t$current_line\n", RESET;
      next();
    }

    $total_reads_parsed++;

    my $read_id = $line[$columns{Read_ID}{position}];
    my $r1_id = $line[$columns{R1_ID}{position}];
    my $r2_id = $line[$columns{R2_ID}{position}];
    my $r1_hit_type = $line[$columns{R1_HitType}{position}];
    my $r2_hit_type = $line[$columns{R2_HitType}{position}];
    my $r1_bit_score = $line[$columns{R1_BitScore}{position}];
    my $r2_bit_score = $line[$columns{R2_BitScore}{position}];
    my $r1_chr = $line[$columns{R1_Chr}{position}];
    my $r2_chr = $line[$columns{R2_Chr}{position}];
    my $r1_chr_start = $line[$columns{R1_ChrStart}{position}];
    my $r1_chr_end = $line[$columns{R1_ChrEnd}{position}];
    my $r2_chr_start = $line[$columns{R2_ChrStart}{position}];
    my $r2_chr_end = $line[$columns{R2_ChrEnd}{position}];

    #change bit scores of 'NA' to 0
    if ($r1_bit_score eq "NA"){$r1_bit_score = 0;}
    if ($r2_bit_score eq "NA"){$r2_bit_score = 0;}

    #If both read alignments are too short, skip this record immediately
    unless ($r1_bit_score >= $min_bit_score || $r2_bit_score >= $min_bit_score){
      next();
    }

    #Test Read1 and Read2 to see if they pass the quality threshold individually
    my $read1_passes = 0;
    my $read2_passes = 0;

    if (($r1_hit_type eq "Top_Hit") && ($r1_bit_score >= $min_bit_score)){
      $read1_passes = 1;
    }
    if (($r2_hit_type eq "Top_Hit") && ($r2_bit_score >= $min_bit_score)){
      $read2_passes = 1;
    }

    #Fix chromosome formats
    if ($r1_chr eq "MT"){$r1_chr = "M";}
    if ($r2_chr eq "MT"){$r2_chr = "M";}

    #Add a count for the gene that this read hits.
    #Deal with READ1
    if ($read1_passes == 1){
      $mapped_reads++;
      $grand_mapped_reads++;

      #Add the coverage of this read to the gene
      &addReadCoverage('-chr'=>$r1_chr, '-chr_start'=>$r1_chr_start, '-chr_end'=>$r1_chr_end);
    }

    #Deal with READ2
    if ($read2_passes == 1){
      $mapped_reads++;
      $grand_mapped_reads++;

      #Add the coverage of this read to the gene
      &addReadCoverage('-chr'=>$r2_chr, '-chr_start'=>$r2_chr_start, '-chr_end'=>$r2_chr_end);
    }

    #After every N successfully mapped reads processed:
    # - calculate the overall coverage of the transcriptome attained
    # - calculate the number of genes with at least Y reads

    if ($mapped_reads >= $block_size){
      &summarizeCoverage();
      $mapped_reads = 0;
    }
  }
  untie @file_array;
}

my $cmd = "rm -f $temp_file";
print YELLOW, "\n\nExecuting: $cmd", RESET;
system($cmd);

if ($file_io_fails){
  print YELLOW, "\n\nEncountered $file_io_fails file I/O failures\n\n", RESET;
}

print OUT "#SCRIPT COMPLETE\n";
close (OUT);

exit();


###########################################################################################################
#Get data files and the columns of each                                                                   #
###########################################################################################################
sub getDataFiles{
  my %args = @_;
  my $dir = $args{'-input_dir'};

  my @required_columns = qw(Read_ID R1_ID R1_HitType R1_BitScore R1_Chr R1_ChrStart R1_ChrEnd R2_ID R2_HitType R2_BitScore R2_Chr R2_ChrStart R2_ChrEnd);

  my $dh = opendir(DIR, $dir) || die "\nCould not open directory: $dir\n\n";

  my @files = readdir(DIR);
  
  #Assign a random number to each file to allow random order processing of lane files
  srand();

  my $count = 0;
  foreach my $file (@files){
    my %columns;
    my $header = 1;
    
    chomp($file);
    unless ($file =~ /\.txt\.gz$/){
      next();
    }
    if (-d $file){
      next();
    }
    $count++;

    $files{$count}{file_name} = $file;
    $files{$count}{file_path} = $dir.$file;
    $files{$count}{random} = rand();

    #Get the header values for this file
    open (FILE, "zcat $dir$file |") || die "\nCould not open file: $dir$file";

    while(<FILE>){
      chomp($_);
      my @line = split("\t", $_);

      #Parse the column names and positions.  Check against a hard coded list of required columns before proceeding
      if ($header == 1){
        my $col_count = 0;
        foreach my $column (@line){
          $columns{$column}{position} = $col_count;
          $col_count++;
        }
        foreach my $req_column (@required_columns){
          unless ($columns{$req_column}){
	    print RED, "\nRequired column: $req_column was not found in the file: $file\n\n", RESET;
	    exit();
          }
        }
        last();
      }

    }
    close(FILE);
    $files{$count}{columns} = \%columns;
  }
  closedir(DIR);

  my $files_count = keys %files;
  print BLUE, "\n\nFound $files_count files to be processed (all .txt.gz files in the specified directory)", RESET;

  return();
}


############################################################################################################################################
#Get basic info from the user specified intron/intergenic region database file                                                             #
############################################################################################################################################
sub importSeqDatabase{
  my %args = @_;
  my $file = $args{'-seq_database'};

  my %seqs;

  open (IN, "zcat $file |") || die "\nCould not open intron/intergenic annotation file: $file\n\n";  

  my $i_count = 0;
  my $header = 1;
  my %columns;

  while(<IN>){
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
    $i_count++;
    my $chr = $line[$columns{'Chromosome'}{column_pos}];
    my $start_chr = $line[$columns{'Unit1_start_chr'}{column_pos}];
    my $end_chr = $line[$columns{'Unit1_end_chr'}{column_pos}];
    my $unmasked_base_count = $line[$columns{'UnMasked_Base_Count'}{column_pos}];
    my $base_count = $line[$columns{'Base_Count'}{column_pos}];

    $seqs{$i_count}{chromosome} = $chr;
    $seqs{$i_count}{start} = $start_chr;
    $seqs{$i_count}{end} = $end_chr;
    $seqs{$i_count}{unmasked_bases} = $unmasked_base_count;
    $seqs{$i_count}{base_count} = $base_count;
  }

  close(IN);

  $seq_ref = \%seqs;

  return();
}


############################################################################################################################################
#Initialize hashes to store Exon coverage for a subset of genes
############################################################################################################################################
sub getCoverageObject{
  my %args = @_;

  my %chr_i_coverage;
  my $total_i_bases = 0;
  my $total_i_unmasked_bases = 0;

  #Group the intron/intergenic regions by chromosome
  foreach my $i (keys %{$seq_ref}){
    my $chr = $seq_ref->{$i}->{chromosome};

    if ($chr_i_coverage{$chr}){
      push(@{$chr_i_coverage{$chr}{i_list}}, $i);
      $chr_i_coverage{$chr}{unmasked_count} += $seq_ref->{$i}->{unmasked_bases};
      $total_i_bases += $seq_ref->{$i}->{base_count};
      $total_i_unmasked_bases += $seq_ref->{$i}->{unmasked_bases};
    }else{
      my @tmp;
      push(@tmp, $i);
      $chr_i_coverage{$chr}{i_list} = \@tmp;
      $chr_i_coverage{$chr}{unmasked_count} = $seq_ref->{$i}->{unmasked_bases};
      $total_i_bases += $seq_ref->{$i}->{base_count};
      $total_i_unmasked_bases += $seq_ref->{$i}->{unmasked_bases};
    }
  }

  #Build a coverage hash for chromosome (keyed on chromosome positions)
  #At this time, also get the chromosome coordinates for all exon-content coordinates
  $| = 1; print BLUE, "\na.) Initializing a COVERAGE hash for the $type CONTENT of each chromosome\n", RESET; $| = 0;

  foreach my $current_chr (sort keys %chr_i_coverage){

    my @chr_i_list = @{$chr_i_coverage{$current_chr}{i_list}};
    my $chr_i_count = scalar(@chr_i_list);
    $| = 1; print BLUE, "\n\tProcessing chr$current_chr (has $chr_i_count $type regions)", RESET; $| = 0;

    $chr_i_coverage{$current_chr}{i_base_count} = &sumCoverage('-current_chr'=>$current_chr, '-i_list'=>\@chr_i_list);
  
    #initialize an empty coverage hash
    my %temp;
    $chr_i_coverage{$current_chr}{coverage} = \%temp;
  }

  my $chr_count = keys %chr_i_coverage;
  my $i_count = keys %{$seq_ref};

  $| = 1; print BLUE, "\n\n\tFound $total_i_bases non-redundant $type bases ($total_i_unmasked_bases UnMasked) in $i_count $type regions on $chr_count chromosomes/contigs", RESET; $| = 0;
  print OUT "#Found $total_i_bases non-redundant $type bases ($total_i_unmasked_bases UnMasked) in $i_count $type regions on $chr_count chromosomes/contigs\n";

  return(\%chr_i_coverage);
}


############################################################################################################################################
#sumCoverage
############################################################################################################################################
sub sumCoverage{
  my %args = @_;
  my @chr_i_list = @{$args{'-i_list'}}; 
  my $current_chr = $args{'-current_chr'};

  my $counter = 0;
  my $highest_coord = 0;
  my $lowest_coord = 100000000000000000000000000000000000000000000;

  my %coverage;

  foreach my $i (@chr_i_list){

    $counter++;
    if ($counter == 100){
      $counter = 0;
      $| = 1; print BLUE, ".", RESET; $| = 0;
    }

    my $start = $seq_ref->{$i}->{start};
    my $end = $seq_ref->{$i}->{end};
    if ($end > $highest_coord){$highest_coord = $end;}
    if ($start < $lowest_coord){$lowest_coord = $start;}

    #Go through each chromosome position in this exon content block and initialize that position in the hash
    for (my $i = $start; $i <= $end; $i++){
      $coverage{$i} = 0;
    }
  }

  my $chr_i_bases = keys %coverage;
  $total_i_bases += $chr_i_bases;
  my $coord_span = ($highest_coord - $lowest_coord)+1;

  $| = 1; print BLUE, "\n\t\tFound $chr_i_bases $type bases covering a total of $coord_span chromosome bases", RESET; $| = 0;

  undef %coverage;
  %coverage = ();
  return($chr_i_bases);
}



############################################################################################################################################
#Add the coverage of a read to a gene to the coverage hash for the exon content record for that gene
############################################################################################################################################
sub addReadCoverage{
  my %args = @_;
  my $chr = $args{'-chr'};
  my $chr_start = $args{'-chr_start'};
  my $chr_end = $args{'-chr_end'};

  #A.) First store coverage at the level of entire chromosomes
  my $coverage_ref = $chr_i_coverage_ref->{$chr}->{coverage};

  #Go through each chromosome position in this read as it is mapped to an exon and increment that position in the hash
  for (my $i = $chr_start; $i <= $chr_end; $i++){
    if ($coverage_ref->{$i}){
      $coverage_ref->{$i}++;
    }else{
      $coverage_ref->{$i} = 1;
    }
  }

  return();
}


################################################################################################################################
#summarize read coverage of transcriptome and number of genes identified after given number of reads
################################################################################################################################
sub summarizeCoverage{

  my $bases_covered = 0;

  #Calculate the exonic base coverage for the whole transcriptome
  foreach my $chr (sort keys %{$chr_i_coverage_ref}){

    my $coverage_ref = $chr_i_coverage_ref->{$chr}->{coverage};

    foreach my $pos (keys %{$coverage_ref}){

      my $cov = $coverage_ref->{$pos};

      if ($cov >= $expressed_base_value){
	$bases_covered++;
      }
    }
  }

  my $bases_diff = $bases_covered - $last_bases_count;
  $last_bases_count = $bases_covered;

  #Determine elapsed time since last report:
  my $t1 = new Benchmark;
  my $td = timediff($t1, $t0);
  $t0 = $t1;

  #Determine current %memory usage
  my $ps_query = `ps u -p $pid`;
  my @process_info = split ("\n", $ps_query);
  my $memory_usage = '';
  if ($process_info[1] =~ /\S+\s+\S+\s+\S+\s+(\S+)\s+/){
    $memory_usage = $1;
  } 

  my $string = timestr($td);
  my $seconds = '';
  if ($string =~ /(\d+)\s+wallclock/){
    $seconds = $1;
  } 

  $| = 1; print BLUE, "\n[READS: $grand_mapped_reads]\t[BASES ($type): $bases_covered ($bases_diff new)]", RESET; $| = 0; 
  $| = 1; print BLUE, "\n\t[Elapsed Time = $seconds seconds]\t[% memory usage = $memory_usage]", RESET; $| = 0; 
  print OUT "$grand_mapped_reads\t$bases_covered\n";

  return();
}
