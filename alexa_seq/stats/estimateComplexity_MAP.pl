#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#Take a directory containing read record files, merge these into one massive temp file, randomly sample this file and assess complexity
#Complexity will be assessed as: 
#  a.) number of unique reads per n reads total
#  b.) number of distinct non-unique reads per n reads total

#NOTE on 'SEQ' vs. 'MAP' approaches to assessing library complexity
#The 'MAP' approach uses uniqueness of the outer mapping positions of each read pair as a measure of library complexity

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;
use Tie::File;
use List::Util 'shuffle';
use Fcntl 'O_RDONLY';

#Load the ALEXA libraries
BEGIN {
  use Cwd 'abs_path';
  if (abs_path($0) =~ /(.*)\/.*\/.*\.pl/){
    push (@INC, $1);
  }
}
use utilities::utility qw(:all);

#Initialize command line options
my $read_records_dir = '';
my $outfile = '';
my $iterations = '';
my $temp_dir = '';
my $ensembl_version = '';
my $block_size = '';

GetOptions ('read_records_dir=s'=>\$read_records_dir, 'temp_dir=s'=>\$temp_dir, 'outfile=s'=>\$outfile, 'iterations=i'=>\$iterations, 'ensembl_version=i'=>\$ensembl_version, 'block_size=i'=>\$block_size);

#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the directory containing read record files using: --read_record_dir", RESET;
print GREEN, "\n\tSpecify a temp directory using:  --temp_dir", RESET;
print GREEN, "\n\tSpecify an output file for results using: --outfile", RESET;
print GREEN, "\n\tSpecify the number of sampling iterations to perform using: --iterations", RESET;
print GREEN, "\n\tSpecify the ensembl_version used for the analysis using:  --ensembl_version", RESET;
print GREEN, "\n\tSpecify the read sample size (e.g. 1,000,000 reads) using: --block_size", RESET;
print GREEN, "\n\nExample: estimateComplexity_MAP.pl  --read_records_dir=/projects/malachig/alexa_seq/read_records/HS04391/  --temp_dir=/projects/malachig/alexa_seq/temp/complexity/  --outfile=/projects/malachig/alexa_seq/figures_and_stats/HS04391/LibraryQuality/LibraryComplexity_MAP.txt  --iterations=3  --ensembl_version=53  --block_size=1000000\n\n", RESET;

unless ($read_records_dir && $temp_dir && $outfile && $iterations && $ensembl_version && $block_size){
  print RED, "\nRequired input parameter(s) missing!\n\n", RESET;
  exit();
}
$| = 1;

if (-e $outfile){
  print YELLOW, "\n\nOutfile: $outfile already exists, aborting...\n", RESET;
  exit();
}

#Global vars
my $pid = $$;
my $c = 0;
my $n_reads = 0;

#Check input directories
$read_records_dir = &checkDir('-dir'=>$read_records_dir, '-clear'=>"no");
$temp_dir = &checkDir('-dir'=>$temp_dir, '-clear'=>"no");

#Get the mapped reads dir
my $mapped_reads_dir = "$read_records_dir"."ENST_v$ensembl_version/";
$mapped_reads_dir = &checkDir('-dir'=>$mapped_reads_dir, '-clear'=>"no");

#Get all files in the mapped reads dir
print BLUE, "\n\n1.) Getting mapped reads files from: $mapped_reads_dir", RESET;
my $dh = opendir(DIR, $mapped_reads_dir) || die "\nCould not open directory: $mapped_reads_dir\n\n";
my @files = readdir(DIR);
my %files;
my $fc = 0;
my $record_count = 0;
foreach my $file (@files){
  my $reads_stored = 0;
  if ($file =~ /(\w+\.txt)\.gz$/){
    $fc++;
    my $name = $1;
    my $file_path = "$mapped_reads_dir"."$file";
    my $temp_file_path = "$temp_dir"."MAP_"."$name";

    #Create a temp copy of the file in the temp dir
    open (READS, "zcat $file_path |") || die "\n\nCould not open input file: $file_path\n\n";
    open(TEMP, ">$temp_file_path") || die "\n\nCould not open temp file for output: $temp_file_path\n\n";
    my $header = 1;
    my %columns;
    while(<READS>){
      chomp($_);
      my @line = split("\t", $_);
      if ($header == 1){
        my $c = 0;
        foreach my $col (@line){
          $columns{$col}{position} = $c;
          $c++;
        }
        $header = 0; 
        next();
      }

      #Get the gene, chr, chr starts, and chr ends for each record
      my $r1_gene_id = $line[$columns{'R1_GeneID'}{position}];
      my $r2_gene_id = $line[$columns{'R2_GeneID'}{position}];
      my $r1_chr = $line[$columns{'R1_Chromosome'}{position}];
      my $r2_chr = $line[$columns{'R2_Chromosome'}{position}];
      my $r1_chr_starts = $line[$columns{'R1_ChrStartCoords'}{position}];
      my $r2_chr_starts = $line[$columns{'R2_ChrStartCoords'}{position}];
      my $r1_chr_ends = $line[$columns{'R1_ChrEndCoords'}{position}];
      my $r2_chr_ends = $line[$columns{'R2_ChrEndCoords'}{position}];

      #Limit the analysis to correctly paired reads
      unless ($r1_gene_id eq $r2_gene_id && $r1_gene_id =~ /\d+/ && $r2_gene_id =~ /\d+/){
        next();
      }
      my @c1 = split(" ", $r1_chr_starts);
      my @c2 = split(" ", $r2_chr_starts);
      my @c3 = split(" ", $r1_chr_ends);
      my @c4 = split(" ", $r2_chr_ends);
      my @coords = (@c1, @c2, @c3, @c4);
      my @sort_coords = sort {$a <=> $b} @coords;
      my $grand_start = $sort_coords[0];
      my $grand_end = @sort_coords[scalar(@sort_coords)-1];

      my $read_record = "$r1_chr"."_"."$grand_start"."-"."$grand_end";
      print TEMP "$read_record\n";
      $record_count++;
      $reads_stored++;
    }
    close(READS);
    close(TEMP);

    if ($reads_stored > $block_size){
      $files{$fc}{name} = $name;
      $files{$fc}{file} = $file;
      $files{$fc}{file_path} = $file_path;
      $files{$fc}{temp_file_path} = $temp_file_path;
      $files{$fc}{reads_stored} = $reads_stored;
      print BLUE, "\n\t$file", RESET;
    }else{
      print YELLOW, "\n\tNot enough reads for $file ($reads_stored reads stored)", RESET;
    }
  }
}
closedir(DIR);


#How many reads should be drawn from each lane?
my $lane_count = keys %files;
unless ($lane_count > 0){
  print RED, "\nNo single lane had more than $block_size valid reads stored - aborting\n\n", RESET;
  exit();
}
my $reads_per_lane = sprintf("%.0f", ($block_size/$lane_count))+1;
print BLUE, "\n\tPrinted a total of $record_count reads to a total of $lane_count files", RESET;

print BLUE, "\n\n2.) Sampling $reads_per_lane reads randomly from each of $lane_count lanes of data until $block_size total reads have been sampled ";

#Now for each iteration, sample from the temp lane files
my @unique_reads;
my @unique_reads_p;
my @non_unique_reads;
my @non_unique_reads_p;
my @distinct_non_unique_reads;
my @distinct_non_unique_reads_p;

print BLUE, "\nIteration\tTotal_Records\tSample_Size\tUnique_Reads\tNon_Unique_Reads\tDistinct_Non_Unique_Reads ", RESET;
for (my $i = 1; $i <= $iterations; $i++){
  my %h = ();
  my $sample_size = 0;

  foreach my $fc (sort {$a <=> $b} keys %files){
    print YELLOW, ".", RESET;
    my $temp_file_path = $files{$fc}{temp_file_path};
    my @file_array;
    tie @file_array, 'Tie::File', "$temp_file_path", mode => O_RDONLY;
    my $line_count = scalar(@file_array)-1;
    my @record_list = (1..$line_count);
    my @shuffled = shuffle(@record_list);
    my @shuffled_slice = @shuffled[1..$reads_per_lane];

    foreach my $i (@shuffled_slice){
      my $current_read = $file_array[$i];
      chomp($current_read);
      #Watch for sporadic failure to retrieve a line
      unless ($current_read){
        print MAGENTA, "\n\tFile I/O error. Retrieved:\n\t$current_read\n", RESET;
        next();
      }

      #Consider reads as a pair
      if ($h{$current_read}){
        $h{$current_read}++;
      }else{
        $h{$current_read}=1;
      }
      $sample_size++;

      #Once the target sample size is achieved, discontinue sampling
      if ($sample_size >= $block_size){
        last();
      }
    }
    untie @file_array;
  }

  #Count the unique reads, non-unique reads, and DISTINCT non-unique reads
  #The distinct non-unique reads give us some sense whether all the non-unique reads are the same thing over and over, or different things.
  my $distinct_non_unique_reads = 0;
  my $unique_reads = 0;
  while (my ($r) = each %h){
    if ($h{$r} > 1){
      $distinct_non_unique_reads++;
    }else{
      $unique_reads++;
    }
  }
  my $unique_reads_p = sprintf("%.2f", ($unique_reads/$sample_size)*100);
  my $non_unique_reads = $sample_size - $unique_reads;
  my $non_unique_reads_p = sprintf("%.2f", ($non_unique_reads/$sample_size)*100);
  my $distinct_non_unique_reads_p = sprintf("%.2f", ($distinct_non_unique_reads/$non_unique_reads)*100);

  push(@unique_reads, $unique_reads);
  push(@unique_reads_p, $unique_reads_p);
  push(@non_unique_reads, $non_unique_reads);
  push(@non_unique_reads_p, $non_unique_reads_p);
  push(@distinct_non_unique_reads, $distinct_non_unique_reads);
  push(@distinct_non_unique_reads_p, $distinct_non_unique_reads_p);

  print BLUE, "\n$i\t$record_count\t$sample_size\t$unique_reads($unique_reads_p%)\t$non_unique_reads($non_unique_reads_p%)\t$distinct_non_unique_reads($distinct_non_unique_reads_p%)", RESET;
}

#Display the averages across all iterations
my $avg_unique_reads = sprintf("%.2f", &avg('-array'=>\@unique_reads));
my $avg_unique_reads_p = sprintf("%.2f", &avg('-array'=>\@unique_reads_p));
my $avg_non_unique_reads = sprintf("%.2f", &avg('-array'=>\@non_unique_reads));
my $avg_non_unique_reads_p = sprintf("%.2f", &avg('-array'=>\@non_unique_reads_p));
my $avg_distinct_non_unique_reads = sprintf("%.2f", &avg('-array'=>\@distinct_non_unique_reads));
my $avg_distinct_non_unique_reads_p = sprintf("%.2f", &avg('-array'=>\@distinct_non_unique_reads_p));

open (OUT, ">$outfile") || die "\n\nCould not open output file: $outfile\n\n";
print BLUE, "\n\nAVG\t$record_count\t$block_size\t$avg_unique_reads($avg_unique_reads_p%)\t$avg_non_unique_reads($avg_non_unique_reads_p%)\t$avg_distinct_non_unique_reads($avg_distinct_non_unique_reads_p%)", RESET;
print OUT "record_count\tblock_size\tavg_unique_reads(%)\tavg_non_unique_reads(%)\tavg_distinct_non_unique_reads(%)\n";
print OUT "$record_count\t$block_size\t$avg_unique_reads($avg_unique_reads_p%)\t$avg_non_unique_reads($avg_non_unique_reads_p%)\t$avg_distinct_non_unique_reads($avg_distinct_non_unique_reads_p%)\n";
close(OUT);


#Clean up the temp dir
foreach my $fc (sort {$a <=> $b} keys %files){
  my $temp_file_path = $files{$fc}{temp_file_path};
  my $cmd = "rm -f $temp_file_path";
  system($cmd);
}


#Summarize memory usage
my $message = &getMemUsage;
print BLUE, "\n\n$message\n\n", RESET;

$| = 0;

exit();


###########################################################################################################################################
#Get average of an array of numbers                                                                                                       #
###########################################################################################################################################
sub avg{
  my %args = @_;
  my @array = @{$args{'-array'}};
  my $n = scalar(@array);
  my $sum = 0;
  foreach my $x (@array){
    $sum += $x;
  }
  my $avg = $sum/$n;

  return($avg);
}


###########################################################################################################################################
#getMemoryUsage                                                                                                                           #
###########################################################################################################################################
sub getMemUsage{
  my $ps_query = `ps -p $pid -o pmem,rss`;
  my @process_info = split ("\n", $ps_query);
  my $memory_usage = '';
  my $memory_usage_p = '';
  if ($process_info[1] =~ /(\S+)\s+(\S+)/){
    $memory_usage_p = $1;
    $memory_usage = $2;
  }
  my $memory_usage_m = sprintf("%.1f", ($memory_usage/1024));
  my $message = "Memory usage: $memory_usage_m Mb ($memory_usage_p%)"; 

  return($message)
}
