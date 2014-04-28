#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#The purpose of this script is to examine coverage of exon junctions as a function of Illumina read depth

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

my $junction_database_file = '';     
my $read_records_dir = '';
my $min_bit_score = '';              #Minimum bit score of BLAST hit for each read to be allowed for read-to-gene mappings
my $block_size = '';                 #Number of reads before coverage is summarized
my $expressed_junction_value = '';   #Number of reads (read equivalents) required for an junction to be considered expressed/detected
my $outfile = '';

GetOptions ('junction_database_file=s'=>\$junction_database_file,
	    'read_records_dir=s'=>\$read_records_dir, 'min_bit_score=f'=>\$min_bit_score,
	    'block_size=i'=>\$block_size, 
	    'expressed_junction_value=i'=>\$expressed_junction_value,
	    'outfile=s'=>\$outfile);

#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify a path to an annotated junction database file using: --junction_database_file\n", RESET;
print GREEN, "\n\tspecify a path to the read record files to be summarized using: --read_records_dir", RESET;
print GREEN, "\n\tSpecify the minimum bit score to a transcript for each read to be considered for the summary using: --min_bit_score", RESET;
print GREEN, "\n\tSpecify the number of reads to process each time before summarizing coverage using: --block_size", RESET;
print GREEN, "\n\tSpecify the number of reads required for an exon to be considered detected using: --expressed_exon_value", RESET;
print GREEN, "\n\tSpecify the outfile summarizing coverage after each read block using: --outfile", RESET;

print GREEN, "\n\nExample: coverageVsDepth_Junction.pl  --junction_database_file=/projects/malachig/sequence_databases/hs_53_36o/exonJunctions/exonJunctions_62mers_annotated.txt.gz  --read_records_dir=/projects/malachig/solexa/read_records/HS0499/Junctions_v53/  --min_bit_score=48.1  --block_size=100000  --expressed_junction_value=10  --outfile=/projects/malachig/solexa/read_records/HS0499/Summary/HS0499_Junction_Coverage_10xCoverage.txt\n\n", RESET;

unless ($junction_database_file && $read_records_dir && $min_bit_score && ($block_size =~ /\d+/) && ($expressed_junction_value =~ /\d+/) && $outfile){
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  exit();
}

#Get list of input files containing data from the specified directory
$read_records_dir = &checkDir('-dir'=>$read_records_dir, '-clear'=>"no");
my %files;
&getDataFiles('-input_dir'=>$read_records_dir);

#Get the process ID for this script
my $pid = $$;

open (OUT, ">$outfile") || die "\nCould not open output file: $outfile\n\n";

#0.) First get all the neccessary junction database info required to perform the analysis
my $last_known_junction_count = 0;
my $last_novel_junction_count = 0;
my %known_junctions;
my %novel_junctions;

&getJunctionInfo('-infile'=>$junction_database_file);


#Print out header info for output file
print OUT "#User specified parameters:\n";
print OUT "#Junction database: $junction_database_file\n";
print OUT "#min_bit_score = $min_bit_score\n";
print OUT "#expressed_junction_value = $expressed_junction_value\n";
print OUT "#read_records_dir = $read_records_dir\n";

print OUT "#MappedReads\tExpressedKnownJunctions\tExpressedNovelJunctions\n";

$| = 1; print BLUE, "\n\n3-a.) Now parsing through read records files and reporting coverage after each block of $block_size reads is added\n", RESET; $| = 0;

my $total_reads_parsed = 0;
my $mapped_reads = 0;
my $grand_mapped_reads = 0;

#Report the processing time for each block of $block_size reads processed
my $t0 = new Benchmark;

#Go through each read record file and parse.
my $temp_file = "$read_records_dir"."coverage_junction_temp.txt";
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
    unless (scalar(@line) == 9){
      print MAGENTA, "\n\tFile I/O error. Retrieved:\n\t$current_line\n", RESET;
      next();
    }

    $total_reads_parsed++;

    my $read_id = $line[$columns{Read_ID}{position}];
    my $hit_type = $line[$columns{HitType}{position}];
    my $junction_id = $line[$columns{Junction_ID}{position}];
    my $bit_score = $line[$columns{BitScore}{position}];

    #change bit scores of 'NA' to 0
    if ($bit_score eq "NA"){$bit_score = 0;}

    #If both read alignments are too short, skip this record immediately
    unless ($bit_score >= $min_bit_score){
      next();
    }

    #Test Read to see if it passes the quality threshold individually
    my $read_passes = 0;

    if (($hit_type eq "Top_Hit") && ($bit_score >= $min_bit_score)){
      $read_passes = 1;
    }

    #Add a count for the junction that this read hits.
    #Deal with READ1
    if ($read_passes == 1){
      $mapped_reads++;
      $grand_mapped_reads++;

      #Count known and novel junctions seperately
      if ($known_junctions{$junction_id}){
        $known_junctions{$junction_id}{c}++;
      }elsif($novel_junctions{$junction_id}){
        $novel_junctions{$junction_id}{c}++;
      }else{
        print RED, "\nJunction ID ($junction_id) not recognized as a known or novel junction\n\n", RESET;
        exit();
      }
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

  my @required_columns = qw(Read_ID HitType Junction_ID AlignmentLength PercentIdentity BitScore Start End Strand);

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
#Get basic info for all known and predicted junctions from the user specified junction database file                                       #
############################################################################################################################################
sub getJunctionInfo{
  my %args = @_;
  my $infile = $args{'-infile'};

  print BLUE, "\n\nGetting junction info from database file: $infile\n", RESET;
  open(JUNCTIONS, "zcat $infile |") || die "\nCould not open junction database file: $infile\n\n";

  #populate the following two data structures
  # %known_junctions;
  # %novel_junctions;

  my @required_columns = qw (Junction_ID Supporting_EnsEMBL_Count);

  my %columns;
  my $header = 1;
  
  while(<JUNCTIONS>){
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
        print RED, "\nRequired column: $req_column was not found in the file: $infile\n\n", RESET;
        exit();
        }
      }
      $header = 0;
      next();
    }

    my $junction_id = $line[$columns{Junction_ID}{position}];
    my $supporting_ensembl_count = $line[$columns{Supporting_EnsEMBL_Count}{position}];

    if ($supporting_ensembl_count > 0){
      $known_junctions{$junction_id}{c} = 0;
    }else{
      $novel_junctions{$junction_id}{c} = 0;
    }
  }

  close(JUNCTIONS);

  #Print out the total number of known and predicted exon junctions from the database file
  my $known_junction_count = keys %known_junctions;
  my $novel_junction_count = keys %novel_junctions;

  $| = 1; print BLUE, "\n\n\tFound $known_junction_count known junctions and $novel_junction_count novel junctions in junction database", RESET; $| = 0;
  print OUT "#Found $known_junction_count known junctions and $novel_junction_count novel junctions in junction database\n";

  return();
}


################################################################################################################################
#summarize read coverage of transcriptome and number of genes identified after given number of reads
################################################################################################################################
sub summarizeCoverage{
  
  my $expressed_known_junction_count = 0;
  my $expressed_novel_junction_count = 0;

  foreach my $junction_id (keys %known_junctions){
    if ($known_junctions{$junction_id}{c} >= $expressed_junction_value){
      $expressed_known_junction_count++;
    }
  }

  foreach my $junction_id (keys %novel_junctions){
    if ($novel_junctions{$junction_id}{c} >= $expressed_junction_value){
      $expressed_novel_junction_count++;
    }
  }

  my $known_junction_diff = $expressed_known_junction_count - $last_known_junction_count;
  my $novel_junction_diff = $expressed_novel_junction_count - $last_novel_junction_count;
  $last_known_junction_count = $expressed_known_junction_count;
  $last_novel_junction_count = $expressed_novel_junction_count;

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
  }else{
    print YELLOW, "\n$string\n", RESET;
  }

  $| = 1; print BLUE, "\n[READS: $grand_mapped_reads]\t[KNOWN JUNCTIONS: $expressed_known_junction_count ($known_junction_diff new)]\t[NOVEL JUNCTIONS: $expressed_novel_junction_count ($novel_junction_diff new)]", RESET; $| = 0; 
  $| = 1; print BLUE, "\n\t[Elapsed Time = $seconds seconds]\t[% memory usage = $memory_usage]", RESET; $| = 0; 
  print OUT "$grand_mapped_reads\t$expressed_known_junction_count\t$expressed_novel_junction_count\n";

  return();
}
