#!/usr/bin/perl -w
#Written by Malachi Griffith and Obi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#The purpose of this script is to process raw solexa read files for a single lane and create a read record file
#It will get _seq.txt files from the specified source directory and concatenate them to the specified output file
#It will also check to make sure the analysis is actually finished (checks for presence of s_*_finished.txt file)
#It will also check to see if there were any warnings recorded for this lane

#Several file types are supported (see elsewhere for formatting of these file types):
#seq - expected naming: s_1_0001_seq.txt (i.e., s_lane_tile_seq.txt)
#qseq - expected naming: s_1_1_0001_seq.txt (i.e., s_lane_read_tile_qseq.txt)
#fastq - expected naming: s_1_1_sequence.txt (i.e., s_lane_read_sequence.txt)

#Compressed files (.bz or .gz) will also be recognized (i.e., s_1_0001_seq.txt.gz; s_1_1_0001_seq.txt.gz; s_1_1_sequence.txt.gz)

#Note that the second read of a read pair is reverse complemented by this script

#The number of ambiguous bases per read will be noted and ambiguous bases will be converted from '.'s to 'N's

#For qseq or fastq files, the mean phred quality score of the read will be calculated.
#There are three types of quality score schemes. See the following urls to determine which applies to your data.
#http://en.wikipedia.org/wiki/FASTQ_format
#http://maq.sourceforge.net/fastq.shtml
#http://maq.sourceforge.net/qual.shtml

#To skip phred score filtering, set phred_cutoff=0.

#It will also use mdust to identify reads with low complexity regions and note the number of mdust bases
#- The default mdust cutoff value of 27 will be used.  This will result in only very low complexity sequences being masked
#- A lower value than this may eliminate some usable reads

#A simple regex will determine how many reads appear to correspond to polyA tails
#  - e.g. look for a stretch of at least 10 consecutive A's or T's and 1/2 or more of the read's overall length is comprised of A's or T's
#  - These reads will simply be summarized, not actually removed or masked

#If the user specifies, the reads will also be trimmed by the specified amount
#- This is done by removing N bases from the end of each read (before reverse complementing R2)

#Basic stats provided for the lane:
#Total number of reads found
#Total number of poor quality reads (N > max ambiguous bases OR P < min phred cutoff)
#Total number of low complexity reads (>= $low_complexity_cutoff mdust bases)
#- The subset of these reads which are actually polyA/T reads (>= $low_complexity_cutoff A's or T's and a stretch of at least 12 consecutive A' or T's)

#Total number of bases found
#Total number (and percent) of ambiguous bases found
#Total number (and percent) of mdust bases found

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;
#use BerkeleyDB;
use Tie::File;
use Fcntl 'O_RDONLY';

#Load the ALEXA modules
BEGIN {
  use Cwd 'abs_path';
  if (abs_path($0) =~ /(.*)\/.*\/.*\.pl/){
    push (@INC, $1);
  }
}
use utilities::utility qw(:all);

#Initialize command line options
my $input_dir = '';
my $flowcell_name = '';
my $lane_number = '';
my $raw_read_dir = '';
my $read_record_dir = '';
my $mdust_bin = '';
my $temp_dir = '';
my $trim = '';
my $low_complexity_cutoff = '';
my $log_dir = '';
my $repair = '';
my $file_type = '';
my $force = '';
my $n_cutoff = '';
my $qual_type = '';
my $phred_cutoff = '';

GetOptions ('input_dir=s'=>\$input_dir, 'flowcell_name=s'=>\$flowcell_name, 'lane_number=i'=>\$lane_number, 'raw_read_dir=s'=>\$raw_read_dir, 
	    'read_record_dir=s'=>\$read_record_dir, 'mdust_bin=s'=>\$mdust_bin, 'temp_dir=s'=>\$temp_dir,
	    'trim=i'=>\$trim, 'low_complexity_cutoff=i'=>\$low_complexity_cutoff,
	    'log_dir=s'=>\$log_dir, 'repair=s'=>\$repair, 'file_type=s'=>\$file_type, 'force=s'=>\$force, 'n_cutoff=i'=>\$n_cutoff, 'phred_cutoff=i'=>\$phred_cutoff, 'qual_type=s'=>\$qual_type);

#Provide instruction to the user
print GREEN, "\n\nNOTE: This script assumes you are processing PAIRED Solexa reads.  It will divide each concatenated read into two seperate reads", RESET;
print GREEN, "\n\tEach read of the pair is presumed to be 1/2 of the read sequence", RESET;
print GREEN, "\n\tThe second read will be reverse complemented", RESET;
print GREEN, "\n\tAll bases reported by Solexa as '.' will be converted to N's\n\n", RESET;

print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the directory containing Solexa seq files using: --input_dir", RESET;
print GREEN, "\n\tSpecify the flowcell name for this input file using: --flowcell_name (this will be appended to read names)", RESET;
print GREEN, "\n\tSpecify the lane number for this input file using: --lane_number (this will be verified against files in the input_dir and used for the log file)", RESET;
print GREEN, "\n\tSpecify a raw read output directory using: --raw_read_dir", RESET;
print GREEN, "\n\tSpecify a read record output directory using: --read_record_dir", RESET;
print GREEN, "\n\t\tThis read record file will have each paired read as well as the number of ambiguous bases for each read, etc.", RESET;
print GREEN, "\n\tSpecify the location of the mdust binary using: --mdust_bin", RESET;
print GREEN, "\n\tSpecify a temp directory for mdust using: --temp_dir", RESET;
print GREEN, "\n\tIf you wish to trim the ends of each read of a pair use the trim option.  (e.g. --trim=1 will remove 1 base from the end of both reads)", RESET;
print GREEN, "\n\tSpecify the maximum number of low complexity bases allowed in a read (reads exceeding this will be classified as low complexity reads) using: --low_complexity_cutoff", RESET;
print GREEN, "\n\tSpecify the maximum number of N (uncalled) bases allowed in a read (reads exceeding this will be classified as low quality reads, default=1) using: --n_cutoff", RESET;
print GREEN, "\n\tSpecify the minimum average phred score allowed for a read (reads exceeding this will be classified as low quality reads) using: --phred_cutoff", RESET;
print GREEN, "\n\t\tSet to 1 (default) for no phred filtering (phred scores still calculated). Set to 0 to skip phred calculations altogether", RESET;
print GREEN, "\n\tSpecify the file type using: --file_type (currently supported: seq | qseq | fastq)", RESET;
print GREEN, "\n\tSpecify the quality score type (for fastq or qseq only) using: --qual_type (illumina_phred | sanger_phred | solexa_illumina)", RESET;
print GREEN, "\n\tAlso specify the directory for a log file using: --log_dir", RESET;
print GREEN, "\n\tIf the mdust code failed but the raw seq file was created successfully you can try the repair option: --repair=yes", RESET;
print GREEN, "\n\tTo prevent user prompts use: --force=yes", RESET;

print GREEN, "\n\nExample: processRawSolexaReads.pl  --input_dir=/archive/solexa1_2/data2/080325_SOLEXA4_0028_20821AAXX/HS04391..L1/Data/C1-74_Firecrest1.9.1_03-04-2008_aldente/Bustard1.9.1_03-04-2008_aldente/  --flowcell_name=20821AAXX  --lane_number=1  --raw_read_dir=/projects/malachig/solexa/raw_seq_data/HS04391/  --read_record_dir=/projects/malachig/solexa/read_records/HS04391/  --mdust_bin=/home/malachig/tools/dust/mdust/mdust  --temp_dir=/projects/malachig/solexa/temp  --trim=0  --low_complexity_cutoff=24  --n_cutoff=1  --phred_cutoff=10  --file_type=qseq  --qual_type=illumina_phred  --log_dir=/projects/malachig/solexa/logs/\n\n", RESET;

#Make sure all options were specified
unless ($input_dir && $flowcell_name && $lane_number && $raw_read_dir && $read_record_dir && $mdust_bin && $temp_dir && ($low_complexity_cutoff =~ /\d+/) && $log_dir && $file_type){
  print RED, "\nRequired input parameter(s) missing or incorrect format!\n\n", RESET;
  exit();
}

#Specify default values for optional arguments
unless($repair){
  $repair = "no";
}

unless($n_cutoff){
  $n_cutoff = 1;
}

#By default do not filter on phred. Phred of 1 is the lowest possible value (corresponds to 100% probability of incorrect base call)
unless($phred_cutoff){
  $phred_cutoff = 1;
}

unless($qual_type){
  $qual_type = "illumina_phred";
}

#Check working dirs before proceeding
#Add trailing '/' to directory paths if they were forgotten
$input_dir = &checkDir('-dir'=>$input_dir, '-clear'=>"no");
$raw_read_dir = &checkDir('-dir'=>$raw_read_dir, '-clear'=>"no");
$read_record_dir = &checkDir('-dir'=>$read_record_dir, '-clear'=>"no");
$log_dir = &checkDir('-dir'=>$log_dir, '-clear'=>"no");
$temp_dir = &checkDir('-dir'=>$temp_dir, '-clear'=>"no");

#Counters for final stats
my $total_readpairs = 0;
my $duplicate_reads = 0;
my $duplicate_reads_rc = 0;
my $total_reads = 0;
my $R1_low_quality_reads = 0;
my $R2_low_quality_reads = 0;
my $R1_low_phred_reads = 0;
my $R2_low_phred_reads = 0;
my $R1_grand_total_phred = 0;
my $R2_grand_total_phred = 0;
my $R1_N_filtered_reads = 0;
my $R2_N_filtered_reads = 0;
my $R1_singleN_reads = 0;
my $R2_singleN_reads = 0;
my $R1_low_complexity_reads = 0;
my $R2_low_complexity_reads = 0;
my $R1_polyAT_reads = 0;
my $R2_polyAT_reads = 0;
my $total_bases = 0;
my $R1_ambiguous_bases = 0;
my $R2_ambiguous_bases = 0;
my $R1_mdust_bases = 0;
my $R2_mdust_bases = 0;

my $read1_a_count = 0;
my $read2_a_count = 0;
my $read1_c_count = 0;
my $read2_c_count = 0;
my $read1_g_count = 0;
my $read2_g_count = 0;
my $read1_t_count = 0;
my $read2_t_count = 0;

my %phred_dist;

#open log file to record the output of this script
my $log_file = "$log_dir"."processSolexaReads_LOG_"."$flowcell_name"."_"."$lane_number".".txt";
open (LOG, ">$log_file") || die "\nCould not open log_file: $log_file\n\n";
print LOG "\nFollowing is a log for the script processRawSolexaRun.pl\n\n";
print LOG "\nUser specified the following options:\ninput_dir = $input_dir\nflowcell_name = $flowcell_name\nlane_number=$lane_number\nraw_read_dir = $raw_read_dir\nread_record_dir = $read_record_dir\nmdust_bin = $mdust_bin\ntemp_fir = $temp_dir\nlow_complexity_cutoff = $low_complexity_cutoff\ntrim = $trim\nlog_dir = $log_dir\n";

if ($trim =~ /\d+/){
  print YELLOW, "\nBoth R1 and R2 will be trimmed by $trim bases\n\n", RESET;
  print LOG "\nBoth R1 and R2 will be trimmed by $trim bases\n\n";
}else{
  $trim = 0;
  print YELLOW, "\nNo --trim option specified, read length will remain unchanged\n\n", RESET;
  print LOG "\nNo --trim option specified, read length will remain unchanged\n\n";
}
my $message = &memoryUsage(); $| = 1; print YELLOW, "\n$message\n", RESET; $| = 0; 
 
#Check source directory for raw files and create a concatenated raw read file and save this to the specified path
my $raw_read_file;
if ($repair =~ /y|yes/i){
  $raw_read_file = "$raw_read_dir"."$flowcell_name"."_Lane"."$lane_number".".txt";
}else{
  #Check the input source directory to ensure that the analysis was complete and look for warnings - also determine the lane number
  my $lane_number_check = &checkSourceDirectory('-input_dir'=>$input_dir, '-flowcell_name'=>$flowcell_name, '-raw_read_dir'=>$raw_read_dir);
  if ($lane_number_check==$lane_number){
    print BLUE, "\n\tThe lane number found matches the lane number specified at runtime", RESET;
    print BLUE, "\n\n", RESET;
    print LOG "\n\tThe lane number found matches the lane number specified at runtime";
    print LOG "\n\n";
  }else{
    print RED, "\nSpecified lane number ($lane_number) does not match lane number found ($lane_number_check)\n\n", RESET;
    exit();
  }
  $raw_read_file = &createRawReadFile('-input_dir'=>$input_dir, '-flowcell_name'=>$flowcell_name, '-raw_read_dir'=>$raw_read_dir, '-lane_number'=>$lane_number);
}

#Create a read record file from the raw read file
$message = &memoryUsage(); $| = 1; print YELLOW, "\n$message\n", RESET; $| = 0; 
my $read_record_file = &createReadRecords('-input_file'=>$raw_read_file, '-flowcell_name'=>$flowcell_name, '-read_record_dir'=>$read_record_dir, '-lane_number'=>$lane_number);

#Use mdust and regular expressions to identify low complexity bases in the read record file
#The read record file created above will be updated
$message = &memoryUsage(); $| = 1; print YELLOW, "\n$message\n", RESET; $| = 0; 
&identifyLowComplexityBases('-input_file'=>$read_record_file, '-temp_dir'=>$temp_dir);

#Print out statistics summarizing this lane
$message = &memoryUsage(); $| = 1; print YELLOW, "\n$message\n", RESET; $| = 0; 
&printStats('-flowcell_name'=>$flowcell_name, '-lane_number'=>$lane_number, '-read_record_dir'=>$read_record_dir);

#Perform a final check of file integrity to ensure that none of the columns are undefined (i.e. an mdust batch may have failed)
$message = &memoryUsage(); $| = 1; print YELLOW, "\n$message\n", RESET; $| = 0; 
open (READS, "$read_record_file") || die "\nCould not open read_record_file: $read_record_file\n\n", RESET;
while(<READS>){
  chomp($_);
  my @line = split("\t", $_);

  #Make sure all values are defined!!  Missing values might indicate file corruption at some point
  unless ($line[0] && $line[1] && $line[2] && $line[3] && $line[4] && $line[5] && $line[6] && ($line[7] =~ /\d+/) && ($line[8] =~ /\d+/) && ($line[9] =~ /\d+/) && ($line[10] =~ /\d+/)){
    print RED, "\n\nFound an undefined value in the input read records file!!!\n\n", RESET;
    print RED, "RECORD: $_\n\n", RESET;
    exit();
  }
}
close (READS);

#Compress the read records file immediately
$message = &memoryUsage(); $| = 1; print YELLOW, "\n$message\n", RESET; $| = 0; 

print BLUE, "\n\nCompressing read record file", RESET;
my $gzip_cmd = "gzip -f $read_record_file";
system ($gzip_cmd);

$message = &memoryUsage(); $| = 1; print YELLOW, "\n$message\n", RESET; $| = 0; 
 
print "\nSCRIPT COMPLETE\n\n";
print LOG "\nSCRIPT COMPLETE\n\n";

close(LOG);

exit();



############################################################################################################
#Perform basic checks on the input source directory
############################################################################################################
sub checkSourceDirectory{
  my %args = @_;
  my $input_dir = $args{'-input_dir'};
  my $flowcell_name = $args{'-flowcell_name'};
  my $raw_read_dir = $args{'-raw_read_dir'};

  my $lane_found;

  #Check the validity of required directories specified by the user
  unless (-d $input_dir && -e $input_dir){
    print RED, "\nInput directory does not appear to be a valid directory:\n\tinput_dir = $input_dir\n\n", RESET;
    exit();
  }

  print BLUE, "\nPerforming basic checks on the input source directory\n\t$input_dir", RESET;
  print LOG "\nPerforming basic checks on the input source directory\n\t$input_dir";

  #See if the specified input directory seems to correspond to the specified flowcell name
  if ($input_dir =~ /$flowcell_name/){
    print BLUE, "\n\n\tThe specified flowcell name ($flowcell_name) appears to match the specified source directory", RESET;
    print LOG "\n\n\tThe specified flowcell name ($flowcell_name) appears to match the specified source directory";
  }else{
    print RED, "\n\n\tThe specified flowcell name ($flowcell_name) doesn't seem to correspond to the specified source directory\n\n", RESET;
    exit();
  }

  #Check for '_finished.txt' and 'finished.txt' files which indicate that the analysis run completed
  my $cmd_a = "ls $input_dir"."finished.txt 2>/dev/null";
  my $cmd_b = "ls $input_dir"."s_"."$lane_number"."_finished.txt 2>/dev/null";
  my $result_a = `$cmd_a`;
  my $result_b = `$cmd_b`;

  my @finished_files_a = split(" ", $result_a);
  my @finished_files_b = split(" ", $result_b);

  if ((scalar(@finished_files_a) + scalar(@finished_files_b)) == 2){
    print BLUE, "\n\tFound two finished.txt files indicating that the analysis of this lane was completed", RESET;
    print LOG "\n\tFound two finished.txt files indicating that the analysis of this lane was completed";
  }else{
    unless ($force){
      print RED, "\n\tDid not find *finished.txt files indicating that the analysis of this lane may not be complete\n\tCheck with LIMS before proceeding!\n\n", RESET;
      print YELLOW, "Proceed anyway? ", RESET;
      my $answer = <>;
      chomp($answer);
      unless($answer =~ /yes|y/i){
        exit();
      }
    }
  }

  opendir(DIRHANDLE, "$input_dir") || die "\nCannot open directory: $input_dir\n\n";
  my @files = readdir(DIRHANDLE);
  my @seq_files;

  #seq files
  if($file_type eq "seq"){
    foreach my $file (@files){
      #Skip all files except _seq.txt files
      my $lane_found;
      if ($file =~ /^s_(\d+)_\d+_seq\.txt$/ || $file =~ /^s_(\d+)_\d+_seq\.txt\.bz2$/ || $file =~ /(\d+)\.seq\.txt\.gz$/ || $file =~ /^s_(\d+)_\d+_seq\.txt\.gz$/){
	$lane_found = $1;
      }else{
	next();
      }
      #skip files corresponding to a lane other than that specified by the user (for directories where all lanes are together)
      unless ($lane_found == $lane_number){
	next();
      }
      push (@seq_files, $file);
    }
  }

  #qseq files
  if($file_type eq "qseq"){
    foreach my $file (@files){
      #Skip all files except _qseq.txt files
      my $lane_found;
      if ($file =~ /^s_(\d+)_\d+_\d+_qseq\.txt$/ || $file =~ /^s_(\d+)_\d+_\d+_qseq\.txt\.bz2$/ || $file =~ /^s_(\d+)_\d+_\d+_qseq\.txt\.gz$/){
        $lane_found = $1;
      }else{
        next();
      }
      #skip files corresponding to a lane other than that specified by the user (for directories where all lanes are together)
      unless ($lane_found == $lane_number){
        next();
      }
      push (@seq_files, $file);
    }
  }

  #fastq files
  if($file_type eq "fastq"){
    foreach my $file (@files){
      #Skip all files except _sequence.txt files
      my $lane_found;
      if ($file =~ /^s_(\d+)_\d+_sequence\.txt$/ || $file =~ /^s_(\d+)_\d+_sequence\.txt\.bz2$/ || $file =~ /^s_(\d+)_\d+_sequence\.txt\.gz$/){
        $lane_found = $1;
      }else{
        next();
      }
      #skip files corresponding to a lane other than that specified by the user (for directories where all lanes are together)
      unless ($lane_found == $lane_number){
        next();
      }
      push (@seq_files, $file);
    }
  }

  closedir(DIRHANDLE);


  #Determine the lane number for the specified source directory
  my $test_file = $seq_files[0];
  if ($test_file =~ /^s_(\d+)_\d+_seq.txt/ || $test_file =~ /^s_(\d+)_\d+_\d+_qseq.txt/ || $test_file =~ /^s_(\d+)_\d+_sequence.txt/ || $test_file =~ /(\d+)\.seq\.txt\.gz$/){
    $lane_found = $1;
  }else{
    print RED, "\nseq file format not understood - could not determine lane number for this source directory\n\n", RESET;
    exit();
  }

  unless ($lane_found == $lane_number){
    print RED, "\nLane number found does not match that supplied by user\n\n", RESET;
    exit();
  }

  foreach my $file (@files){

    #Copy any warning files over to the raw_read_dir
    if ($file =~ /warning/i){

      print YELLOW, "\n\tFound a warning file: $file", RESET;
      print LOG "\n\tFound a warning file: $file";
      my $cmd = "cp $input_dir"."$file". " $raw_read_dir"."$flowcell_name"."_"."Lane$lane_number"."_"."$file";
      system($cmd);
    }
  }

  print BLUE, "\n\tThe lane number use to get files from this source directory was $lane_number", RESET;
  print LOG "\n\tThe lane number use to get files from this source directory was $lane_number";
  return($lane_number);
}


############################################################################################################
#Create a raw read file by concatenating raw files in the source directory
############################################################################################################
sub createRawReadFile{
  my %args = @_;
  my $input_dir = $args{'-input_dir'};
  my $flowcell_name = $args{'-flowcell_name'};
  my $raw_read_dir = $args{'-raw_read_dir'};
  my $lane_number = $args{'-lane_number'};

  my $raw_read_file = "$raw_read_dir"."$flowcell_name"."_Lane"."$lane_number".".txt";

  my %seq_files;

  print BLUE, "\nCreating raw_read_file: $raw_read_file\n\n", RESET;
  print LOG "\nCreating raw_read_file: $raw_read_file\n\n";

  #If the file is already there, delete it
  if (-e $raw_read_file && -s $raw_read_file){

    if ($force){
      my $cmd = "rm -f $raw_read_file";
      print YELLOW, "\nDeleting $raw_read_file\n\n", RESET;
      system($cmd);
    }else{
      print YELLOW, "\nRaw read file already exists.  Overwrite (y/n)? ", RESET;
      my $answer = <>;
      chomp($answer);

      if ($answer =~ /^y$|^yes$/i){
        my $cmd = "rm -f $raw_read_file";
        print YELLOW, "\nDeleting $raw_read_file\n\n", RESET;
        system($cmd);
      }else{
        print RED, "\nNo point in continuing... exiting\n\n", RESET;
        exit();
      }
    }
  }

  opendir(DIRHANDLE, "$input_dir") || die "\nCannot open directory: $input_dir\n\n";
  my @files = readdir(DIRHANDLE);

  print BLUE, "\nCreating concatenated raw sequence file:\n\t$raw_read_file\n\n", RESET;
  print LOG "\nCreating concatenated raw sequence file:\n\t$raw_read_file\n\n";

  my $gr = 0;
  my $r = 0;
  my $r_block = 1000000;

  ####################################################################################
  #For seq files                                                                     #
  ####################################################################################
  if ($file_type eq "seq"){
    foreach my $file (@files){
      #Skip all files except _seq.txt files
      my $lane_found;
      my $tile_found;
      if ($file =~ /^s_(\d+)_(\d+)_seq\.txt$/ || $file =~ /^s_(\d+)_(\d+)_seq\.txt\.bz2$/ || $file =~ /^s_(\d+)_(\d+)_seq\.txt\.gz$/){
     	$lane_found = $1;
	$tile_found = $2;
      }elsif($file =~ /(\d+)\.seq\.txt\.gz$/){
	$lane_found = $1;
	$tile_found = 1;
      }else{
	next();
      }
      #skip files corresponding to a lane other than that specified by the user (for directories where all lanes are together)
      unless ($lane_found == $lane_number){
	next();
      }
      $seq_files{$tile_found}{file} = $file;
    }

    foreach my $tile (sort {$a cmp $b} keys %seq_files){
      my $file = $seq_files{$tile}{file};
      print YELLOW, "\n\tCat: $file", RESET;

      #If file is a simple _seq.txt file just cat it to the raw_read_file
      if ($file =~ /_seq\.txt$/){
	my $cmd = "cat $input_dir"."$file >> $raw_read_file";
	system ($cmd);
      }

      #If file is a _seq.txt.bz2 file, copy it to the raw_read_dir, bunzip2 it, cat it to the raw_read_file, and then delete it
      if ($file =~ /_seq\.txt\.bz2$/){
	my $cmd = "bzcat $input_dir"."$file >> $raw_read_file";
	system ($cmd);

      }
      if ($file =~ /(\d+)\.seq\.txt\.gz$/ || $file =~ /^s_(\d+)_(\d+)_seq\.txt\.gz$/){
	my $cmd = "zcat $input_dir"."$file >> $raw_read_file";
	system ($cmd);
      }
    }
  }


  ####################################################################################
  #For qseq files - need to make temp copies of these to join paired reads together  #
  ####################################################################################
  if ($file_type eq "qseq"){

    my $temp_raw_read_dir = "$raw_read_dir"."$flowcell_name"."_"."Lane$lane_number/";
    my $rm_cmd = "rm -fr $temp_raw_read_dir";
    system($rm_cmd);
    mkdir ($temp_raw_read_dir);

    my @tempfiles;
    open (RAWREADS, ">$raw_read_file") or die "can't open $raw_read_file for write\n";
    foreach my $file (@files){
      #Skip all files except _qseq.txt files
      my $lane_found;
      my $tile_found;
      my $read1or2_found;
      if ($file =~ /^s_(\d+)_(\d+)_(\d+)_qseq\.txt$/ || $file =~ /^s_(\d+)_(\d+)_(\d+)_qseq\.txt\.bz2$/ || $file =~ /^s_(\d+)_(\d+)_(\d+)_qseq\.txt\.gz$/){
	$lane_found = $1;
	$read1or2_found=$2;
	$tile_found = $3;
      }else{
	next();
      }
      #skip files corresponding to a lane other than that specified by the user (for directories where all lanes are together)
      unless ($lane_found == $lane_number){
	next();
      }

      my $filepath = "$input_dir"."$file";

      if ($file =~ /_qseq\.txt\.bz2$/){
	my $cp_cmd = "cp "."$input_dir"."$file"." $temp_raw_read_dir";
	my $bunzip2_cmd = "bunzip2 "."$temp_raw_read_dir"."$file";
	system($cp_cmd);
	system($bunzip2_cmd);
	$file=~s/\.bz2$//; #the bunzip2'd file will have a new name
	$filepath = "$temp_raw_read_dir"."$file"; #Change filepath to use bunzipped version in temp dir
	push(@tempfiles,$filepath);
      }

      if ($file =~ /_qseq\.txt\.gz$/){
	my $cp_cmd = "cp "."$input_dir"."$file"." $temp_raw_read_dir";
	my $gunzip_cmd = "gunzip "."$temp_raw_read_dir"."$file";
	system($cp_cmd);
	system($gunzip_cmd);
	$file=~s/\.gz$//; #the gunzip'd file will have a new name
	$filepath = "$temp_raw_read_dir"."$file"; #Change filepath to use bunzipped version in temp dir
	push(@tempfiles,$filepath);
      }

      $seq_files{$tile_found}{$read1or2_found} = $filepath;
    }

    foreach my $tile (sort {$a cmp $b} keys %seq_files){
      my $read1file = $seq_files{$tile}{'1'};
      my $read2file = $seq_files{$tile}{'2'};
      print YELLOW, "Combining $read1file and $read2file\n", RESET;

      my %seq_data;
      open (FILE1, $read1file) or die "can't open $read1file\n";
      while (<FILE1>){
	my @data = split ("\t",$_);
	my $read = "$data[2]\t$data[3]\t$data[4]\t$data[5]";
	my $sequence = $data[8];
	my $quality_string = $data[9];
	my $mean_phred = &calculateMeanPhredScore('-quality_string'=>$quality_string, '-read'=>"R1");
	$seq_data{$read}{R1}=$sequence;
	$seq_data{$read}{Q1}=$mean_phred;
      }
      close FILE1;

      open (FILE2, $read2file) or die "can't open $read2file\n";
      while (<FILE2>){
	my @data = split ("\t",$_);
	my $read = "$data[2]\t$data[3]\t$data[4]\t$data[5]";
	my $sequence2 = $data[8];
	my $quality_string2 = $data[9];
	my $mean_phred2 = &calculateMeanPhredScore('-quality_string'=>$quality_string2, '-read'=>"R2");

        #Print combined reads to raw reads file
	print RAWREADS "$read\t$seq_data{$read}{R1}"."$sequence2\t"."$seq_data{$read}{Q1}\t$mean_phred2\n";
      }
      close FILE2;
    }
    close RAWREADS;

    #Delete the temp dir
    system($rm_cmd);

  }

  ####################################################################################
  #For fastq files - need to make temp copies of these to join paired reads together #
  ####################################################################################
  if ($file_type eq "fastq"){

    my $temp_raw_read_dir = "$raw_read_dir"."$flowcell_name"."_"."Lane$lane_number/";
    my $rm_cmd = "rm -fr $temp_raw_read_dir";
    system($rm_cmd);
    mkdir ($temp_raw_read_dir);

    my @tempfiles;
    open (RAWREADS, ">$raw_read_file") or die "can't open $raw_read_file for write\n";
    foreach my $file (@files){
      #Skip all files except _sequence.txt files
      my $lane_found;
      my $read1or2_found;
      if ($file =~ /^s_(\d+)_(\d+)_sequence\.txt$/ || $file =~ /^s_(\d+)_(\d+)_sequence\.txt\.bz2$/ || $file =~ /^s_(\d+)_(\d+)_sequence\.txt\.gz$/){
	$lane_found = $1;
	$read1or2_found=$2;
      }else{
	next();
      }
      #skip files corresponding to a lane other than that specified by the user (for directories where all lanes are together)
      unless ($lane_found == $lane_number){
	next();
      }

      my $filepath = "$input_dir"."$file";

      if ($file =~ /_sequence\.txt\.bz2$/){
        #my $cp_cmd = "cp "."$input_dir"."$file"." $temp_raw_read_dir";
        #my $bunzip2_cmd = "bunzip2 "."$temp_raw_read_dir"."$file";
        #system($cp_cmd);
        #system($bunzip2_cmd);
        #$file=~s/\.bz2$//; #the bunzip2'd file will have a new name
        #$filepath = "$temp_raw_read_dir"."$file"; #Change filepath to use bunzipped version in temp dir
        my $new_file;
        if ($file =~ /(.*)\.bz2/){
          $new_file = $1;
        }
	my $bzcat_cmd = "bzcat "."$input_dir"."$file"." > $temp_raw_read_dir"."$new_file";
        system($bzcat_cmd);
	$filepath = "$temp_raw_read_dir"."$new_file"; #Change filepath to use bunzipped version in temp dir

	push(@tempfiles,$filepath);
      }

      if ($file =~ /_sequence\.txt\.gz$/){
	my $cp_cmd = "cp "."$input_dir"."$file"." $temp_raw_read_dir";
	my $gunzip_cmd = "gunzip "."$temp_raw_read_dir"."$file";
	system($cp_cmd);
	system($gunzip_cmd);
	$file=~s/\.gz$//; #the gunzip'd file will have a new name
	$filepath = "$temp_raw_read_dir"."$file"; #Change filepath to use bunzipped version in temp dir
	push(@tempfiles,$filepath);
      }
      $seq_files{$read1or2_found} = $filepath;
    }

    my $read1file = $seq_files{'1'};
    my $read2file = $seq_files{'2'};
    print YELLOW, "Combining $read1file and $read2file\n", RESET;
    $message = &memoryUsage(); $| = 1; print YELLOW, "\n$message\n", RESET; $| = 0; 
    
    #Assume the fastq files are ordered and of equal length
    my @file_array1;
    tie @file_array1, 'Tie::File', "$read1file", mode => O_RDONLY, memory => 0; #Since each record is being used once, caching is pointless, thus set memory to '0'
    #my $line_count1 = scalar(@file_array1)-1;
    #my $line_count1 = $#file_array1;
    my $line_count1 = `wc -l $read1file`; chomp($line_count1);
    if ($line_count1 =~ /^(\d+)/){$line_count1=$1;}else{print RED, "\n\nLine count on file 1 failed\n\n", RESET; exit();}
    
    my @file_array2;
    tie @file_array2, 'Tie::File', "$read2file", mode => O_RDONLY, memory => 0; #Since each record is being used once, caching is pointless, thus set memory to '0'
    #my $line_count2 = scalar(@file_array2)-1;
    #my $line_count2 = $#file_array2;
    my $line_count2 = `wc -l $read2file`; chomp($line_count2);
    if ($line_count2 =~ /^(\d+)/){$line_count2=$1;}else{print RED, "\n\nLine count on file 2 failed\n\n", RESET; exit();}

    if ($line_count1 == $line_count2){
      print YELLOW, "\tFound $line_count1 lines in both files\n", RESET;
    }else{
      print RED, "\n\nR1 and R2 file line counts are not equal!  Aborting\n\n", RESET;
      exit();
    }

    $message = &memoryUsage(); $| = 1; print YELLOW, "\n$message\n", RESET; $| = 0; 
    for (my $i = 0; $i < $line_count1; $i+=4){
      $r++;
      $gr++;
      if ($r == $r_block){
        $r = 0;
        $message = &memoryUsage(); $| = 1; print YELLOW, "\n\tProcessed $gr reads so far - $message", RESET; $| = 0; 
      }
      my $read1;
      my $read1_line = $file_array1[$i];
      my $sequence1 = $file_array1[$i+1];
      my $quality_string1 = $file_array1[$i+3];
      my $read2;
      my $read2_line = $file_array2[$i];
      my $sequence2 = $file_array2[$i+1];
      my $quality_string2 = $file_array2[$i+3];
      chomp($read1_line); chomp($sequence1); chomp($quality_string1);
      chomp($read2_line); chomp($sequence2); chomp($quality_string2);
      if($read1_line =~ /\@\S+\:(\d+)\:(\d+)\:(\d+)\:(\d+)\#\w+\/(\d+)/){
         $read1 = "$1\t$2\t$3\t$4";
      }elsif($read1_line =~ /\@\S+\:(\d+)\:(\d+)\:(\d+)\:(\d+)\/(\d+)/){
         $read1 = "$1\t$2\t$3\t$4";
      }elsif($read1_line =~ /\@\S+\:\d+\:\w+\:(\d+)\:(\d+)\:(\d+)\:(\d+)\s(\d+)\:\S+/){
         $read1 = "$1\t$2\t$3\t$4";
      }else{
        print "Unexpected fastq read 1 name format: $read1_line\n"; 
        exit();
      }
      if($read2_line =~ /\@\S+\:(\d+)\:(\d+)\:(\d+)\:(\d+)\#\w+\/(\d+)/){
         $read2 = "$1\t$2\t$3\t$4";
      }elsif($read2_line =~ /\@\S+\:(\d+)\:(\d+)\:(\d+)\:(\d+)\/(\d+)/){
         $read2 = "$1\t$2\t$3\t$4";
      }elsif($read1_line =~ /\@\S+\:\d+\:\w+\:(\d+)\:(\d+)\:(\d+)\:(\d+)\s(\d+)\:\S+/){
         $read1 = "$1\t$2\t$3\t$4";
      }else{
        print "Unexpected fastq read 2 name format: $read2_line\n";
        exit();
      }
      if($sequence1=~/[^ACTGN\.]/){print "Unexpected sequence: $sequence1\n"; exit;}
      if($sequence2=~/[^ACTGN\.]/){print "Unexpected sequence: $sequence2\n"; exit;}
      my $mean_phred1 = &calculateMeanPhredScore('-quality_string'=>$quality_string1, '-read'=>"R1");
      my $mean_phred2 = &calculateMeanPhredScore('-quality_string'=>$quality_string2, '-read'=>"R2");
      #Print combined reads to raw reads file
      print RAWREADS "$read1\t"."$sequence1"."$sequence2\t"."$mean_phred1\t$mean_phred2\n";
 
    }
    untie @file_array1;
    untie @file_array2;
    close RAWREADS;

    #Delete the temp dir
    system($rm_cmd);
  }

  ####################################################################################
  #Export files                                                                      #
  ####################################################################################
  if ($file_type eq "export"){



  }



  return($raw_read_file);
}

############################################################################################################
#Open input file and write to read output file
############################################################################################################
sub createReadRecords{
  my %args = @_;
  my $input_file = $args{'-input_file'};
  my $flowcell_name = $args{'-flowcell_name'};
  my $read_record_dir = $args{'-read_record_dir'};
  my $lane_number = $args{'-lane_number'};

  my $read_record_file = "$read_record_dir"."$flowcell_name"."_Lane"."$lane_number".".txt";;

  print BLUE, "\n\nCreating a read record file (converting ambiguous bases from Solexa to Ns and counting them): $read_record_file\n\n", RESET;
  print LOG "\n\nCreating a read record file (converting ambiguous bases from Solexa to Ns and counting them): $read_record_file\n\n";

  open (IN_FILE, "$input_file") || die "\nCould not open input_file: $input_file\n\n";
  open (READ_RECORD_FILE, ">$read_record_file") || die "\nCould not open read_record_file: $read_record_file\n\n";

  print READ_RECORD_FILE "Read_ID\tRead1_ID\tRead2_ID\tRead1_Status\tRead2_Status\tRead1_Seq\tRead2_seq\tRead1_Ambig_Count\tRead2_Ambig_Count\n";
  my $line_count = 0;
  while(<IN_FILE>){

    $line_count++;
    chomp($_);
    my @line = split("\t", $_);

    my $lane = $line[0];
    my $tile = $line[1];
    my $x_coord = $line[2];
    my $y_coord = $line[3];
    my $full_seq = $line[4];

    #Make sure each of the 4 components of the ID from above are integers!
    if ($lane =~ /(\d+)/){
      $lane = $1;
    }else{
      print RED, "\nProblem with form of record: $_\n\n", RESET;
      exit();
    }
    if ($tile =~ /(\d+)/){
      $tile = $1;
    }else{
      print RED, "\nProblem with form of record: $_\n\n", RESET;
      exit();
    }
    if ($x_coord =~ /(\d+)/){
      $x_coord = $1;
    }else{
      print RED, "\nProblem with form of record: $_\n\n", RESET;
      exit();
    }
    if ($y_coord =~ /(\d+)/){
      $y_coord = $1;
    }else{
      print RED, "\nProblem with form of record: $_\n\n", RESET;
      exit();
    }

    my $read_base = "$flowcell_name"."_"."$lane"."_"."$tile"."_"."$x_coord"."_"."$y_coord";
    my $read1_name = "$read_base"."_R1";
    my $read2_name = "$read_base"."_R2";

    my $seq_length = length($full_seq);
    my $read_length = ($seq_length/2);

    #Change all '.' values to 'N'
    $full_seq =~ tr/./N/;

    #Get each read of the pair and apply the trim
    my $read1 = substr($full_seq, 0, $read_length-$trim);
    my $read2 = substr($full_seq, $read_length, $read_length-$trim);
    my $original_read2 = $read2;

    #reverse complement the second read
    #$read2 =~ tr/gatcnGATCN/ctagnCTAGN/;
    #$read2 = reverse $read2;

    #Sanity check
    unless ((length($read1) == length($read2)) && (length($read2) == ($read_length - $trim))){
      print RED, "\nSome problem with dividing reads in half\n\n", RESET;
      exit();
    }

    my $read1_status = "Unassigned";
    my $read2_status = "Unassigned";

    #Count the N's in each read
    my $read1_n_count = ($read1 =~ tr/N/N/);
    my $read2_n_count = ($read2 =~ tr/N/N/);

    #Update stats
    $total_readpairs++;
    $total_reads += 2;

    #Consider any read with more than max ambiguous base to be low quality
    if ($read1_n_count > $n_cutoff){
      $R1_low_quality_reads++;
      $R1_N_filtered_reads++;
      $read1_status = "Low_Quality";
    }
    if ($read2_n_count > $n_cutoff){
      $R2_low_quality_reads++;
      $R2_N_filtered_reads++;
      $read2_status = "Low_Quality";
    }

    #For fastq and qseq files, consider any read with less than min mean phred score to be low quality. 
    #Only bother if not already set determined as low quality above
    my ($mean_phred1, $mean_phred2);
    if ($file_type eq "qseq" || $file_type eq "fastq"){
      $mean_phred1 = $line[5];
      $mean_phred2 = $line[6];
      if ($mean_phred1 =~ /(\d+\.?\d*)/){
        $mean_phred1 = $1;
      }else{
        print RED, "\nProblem with form of record: $_\n\n", RESET;
        exit();
      }
      if ($mean_phred2 =~ /(\d+\.?\d*)/){
        $mean_phred2 = $1;
      }else{
        print RED, "\nProblem with form of record: $_\n\n", RESET;
        exit();
      }
      #keep running sum of all phred scores so that grand mean phred score can be determined
      $R1_grand_total_phred += $mean_phred1;
      $R2_grand_total_phred += $mean_phred2;

      if ($mean_phred1 < $phred_cutoff){
        $R1_low_phred_reads++;
        unless ($read1_status eq "Low_Quality"){ #Only count as low-quality read if not already assigned as such by number of Ns
          $R1_low_quality_reads++;
          $read1_status = "Low_Quality";
        }
      }
      if ($mean_phred2 < $phred_cutoff){
        $R2_low_phred_reads++;
        unless ($read2_status eq "Low_Quality"){ #Only count as low-quality read if not already assigned as such by number of Ns
          $R2_low_quality_reads++;
          $read2_status = "Low_Quality";
        }
      }
    } 

    #Print out the actual read record
    print READ_RECORD_FILE "$read_base\t$read1_name\t$read2_name\t$read1_status\t$read2_status\t$read1\t$read2\t$read1_n_count\t$read2_n_count\n";

    $total_bases += (length($read1) + length($read2));
    $R1_ambiguous_bases += $read1_n_count;
    $R2_ambiguous_bases += $read2_n_count;

    #Count the occurence of each base in each read
    my $r1 = $read1;
    my $r2 = $read2;
    $read1_a_count += ($r1 =~ tr/A/A/);
    $read2_a_count += ($r2 =~ tr/A/A/);
    $read1_c_count += ($r1 =~ tr/C/C/);
    $read2_c_count += ($r2 =~ tr/C/C/);
    $read1_g_count += ($r1 =~ tr/G/G/);
    $read2_g_count += ($r2 =~ tr/G/G/);
    $read1_t_count += ($r1 =~ tr/T/T/);
    $read2_t_count += ($r2 =~ tr/T/T/);
  }
  close (IN_FILE);
  close (READ_RECORD_FILE);

  return($read_record_file);
}

###################################################################################################################
#Calculate mean phred score
###################################################################################################################
sub calculateMeanPhredScore{
  my %args = @_;
  my $quality_string = $args{'-quality_string'};
  my $read = $args{'-read'};
  my @qual_chrs;
  my @phred_scores;
  chomp $quality_string;

  #If phred cutoff is set to 0, skip mean_phred calculation
  if ($phred_cutoff==0){
    my $mean_phred=0;
    return($mean_phred);
  }

  #Extract quality characters from quality string
  for (my $i = 0; $i < length($quality_string); $i++) {
    my $q = substr ($quality_string, $i, 1)."\n";
    push (@qual_chrs, $q);
  }

  #Check option for what kind of quality data you have. Then calculate phred score
  #Illumina/phred (aka phred64) (Default assumed by script)
  if ($qual_type eq "illumina_phred"){
    foreach my $q (@qual_chrs){
      my $Q = ord($q) - 64;
      push (@phred_scores, $Q);
    }
  }

  #Sanger/phred (aka phred33)
  if ($qual_type eq "sanger_phred"){
    foreach my $q (@qual_chrs){
      my $Q = ord($q) - 33;
      push (@phred_scores, $Q);
    }
  }

    #Solexa/Illumina (Calculation from: http://maq.sourceforge.net/qual.shtml)
  if ($qual_type eq "solexa_illumina"){
    foreach my $q (@qual_chrs){
      my $Q = 10 * log(1 + 10 ** ((ord($q) - 64) / 10.0)) / log(10);
      push (@phred_scores, $Q);
    }
  }
  
  #Calculate mean phred score
  my $mean_phred = eval(join("+", @phred_scores)) / @phred_scores;

  #Store mean phred score counts (for R1 and R2) for later plotting of distribution
  my $rounded = sprintf ("%d", $mean_phred);
  $phred_dist{$rounded}{$read}++;

  return($mean_phred);
}


###################################################################################################################
#Use mdust to determine the number of low complexity bases in each read in the newly created read file
###################################################################################################################
sub identifyLowComplexityBases{
  my %args = @_;
  my $read_record_file = $args{'-input_file'};
  my $temp_dir = $args{'-temp_dir'};

  print BLUE, "\nRunning mdust on each sequence and updating the read record file with this info\n\n", RESET;
  print LOG "\nRunning mdust on each sequence and updating the read record file with this info\n\n";

  my $line_count = 0;
  my $block_size = 50000;
  my $block_count = 0;

  #Create a $temp_dir called 'mdust_flowcell_lane'
  my $new_dir_name = "mdust_"."$flowcell_name"."_"."$lane_number";
  my $new_dir_path = "$temp_dir"."$new_dir_name";
  my $working_dir = &createNewDir('-path'=>$temp_dir, '-new_dir_name'=>$new_dir_name, '-force'=>"yes");

  #Specify a directory where probe files can be temporarily be created and fed into mdust
  my $mdust_tempfile = "$working_dir"."temp.fa";
  my $results_file = "$working_dir"."reads.txt";

  open (READ_RECORD_FILE, "$read_record_file") || die "\nCould not open read file: $read_record_file\n\n";
  open (RESULTS_FILE, ">$results_file") || die "\nCould not open temp results file: $results_file\n\n";

  #Go through the readfile and get the probe sequences to be processed in blocks
  my %reads;
  while (<READ_RECORD_FILE>){
    $line_count++;
    #Skip the header line
    if ($line_count == 1){
      chomp($_);
      print RESULTS_FILE "$_\tRead1_mdust_bases\tRead2_mdust_bases\n";
      next();
    }
    chomp($_);
    my @line = split("\t", $_);

    #Get the read IDs and sequences from each line
    my $id = $line[0];
    my $id1 = $line[1];
    my $id2 = $line[2];
    my $status1 = $line[3];
    my $status2 = $line[4];
    my $read1 = $line[5];
    my $read2 = $line[6];
    my $ambig1 = $line[7];
    my $ambig2 = $line[8];

    $reads{$id}{id1} = $id1;
    $reads{$id}{id2} = $id2;
    $reads{$id}{status1} = $status1;
    $reads{$id}{status2} = $status2;
    $reads{$id}{sequence1} = $read1;
    $reads{$id}{sequence2} = $read2;
    $reads{$id}{ambig1} = $ambig1;
    $reads{$id}{ambig2} = $ambig2;
    $reads{$id}{line_count} = $line_count;

    #Once a full set of probes has been collected process them with mdust
    if ((keys %reads) == $block_size){
      $block_count++;
      $message = &memoryUsage(); $| = 1; print YELLOW, "\n\tProcessed $block_count blocks of $block_size reads so far - $message", RESET; $| = 0; 
      print LOG "\n\tProcessed $block_count blocks of $block_size reads so far - $message";
      &processReads('-read_hash'=>\%reads, '-mdust_tempfile'=>$mdust_tempfile);
      %reads = ();
    }
  }
  #Process the final incomplete block of probes (as long as there was at least 1 read left)
  my $count = keys %reads;
  if ($count > 1){
    &processReads('-read_hash'=>\%reads, '-mdust_tempfile'=>$mdust_tempfile);
  }

  close(RESULTS_FILE);
  close(READ_RECORD_FILE);

  #Overwrite the original READFILE with the new RESULTS file
  my $cmd = "mv $results_file $read_record_file";
  print BLUE, "\nExecuting: $cmd\n\n", RESET;
  print LOG "\nExecuting: $cmd\n\n";
  system("$cmd");

  #Deleted the $temp_dir called 'mdust_flowcell_lane'
  my $rm_cmd = "rm -fr $new_dir_path";
  print BLUE, "\nExecuting: $rm_cmd\n\n", RESET;
  print LOG "\nExecuting: $rm_cmd\n\n";
  system("$rm_cmd");

  return();
}


###################################################################################################
#Process Probes with mdust
###################################################################################################
sub processReads{
  my %args = @_;
  my $read_list_ref = $args{'-read_hash'};
  my $mdust_tempfile = $args{'-mdust_tempfile'};

  my $mdust_cutoff = 27;

  #Create a temp fasta file with the probe sequences in this block
  open (FASTA_FILE, ">$mdust_tempfile") || die "\nCould not open $mdust_tempfile\n\n";
  foreach my $read_id (sort keys %{$read_list_ref}){
    print FASTA_FILE ">$read_list_ref->{$read_id}->{id1}\n$read_list_ref->{$read_id}->{sequence1}\n";
    print FASTA_FILE ">$read_list_ref->{$read_id}->{id2}\n$read_list_ref->{$read_id}->{sequence2}\n";
  }
  close (FASTA_FILE);

  my $result = `$mdust_bin $mdust_tempfile -v $mdust_cutoff`;

  while (!$result){
    sleep(1);
    print YELLOW, "\nmdust exec failure ... retrying", RESET;
    $result = `$mdust_bin $mdust_tempfile -v $mdust_cutoff`;
  }


  #Split the resulting output back into single probe entries
  my @dusted_results = split (">", $result);

  #Remove the first entry which is empty
  shift(@dusted_results);

  #Grab the probe ID and sequence again
  foreach my $entry (@dusted_results){
#    if ($entry =~ /^([a-zA-Z0-9]+_\d+_\d+_\d+_\d+)_R(\d+)(.*)/s){
    if ($entry =~ /^(\S+_\d+_\d+_\d+_\d+)_R(\d+)(.*)/s){
      my $read_id = $1;
      my $read_num = $2;
      my $seq = $3;
      chomp($seq);

      #clean-up sequence data by removing spaces, >,  and making uppercase
      $seq =~ s/[\s\>]//g;
      my $upper_seq = uc($seq);

      #Count the number of N's
      my @dusted_n_count = $upper_seq =~ m/N/g;
      my $dusted_n_count = @dusted_n_count;

      if ($read_num eq "1"){
	$read_list_ref->{$read_id}->{mdust_count1} = $dusted_n_count;
      }elsif ($read_num eq "2"){
	$read_list_ref->{$read_id}->{mdust_count2} = $dusted_n_count;

      }else{
	print RED, "\nRead num not understood\n\n", RESET;
	exit();
      }

    }else{
      print RED, "\nCould not understand output format of entry from mdust!\n\tENTRY: $entry\n\n", RESET;
      close (LOG);
      exit();
    }
  }

  #Print out the probe info and append the mdust cutoff and number of dust masked N's (representing simple repeats) found
  foreach my $read_id (sort {$read_list_ref->{$a}->{line_count} <=> $read_list_ref->{$b}->{line_count}} keys %{$read_list_ref}){
    my $read1 = $read_list_ref->{$read_id}->{sequence1};
    my $read2 = $read_list_ref->{$read_id}->{sequence2};
    my $mdust_count1 = $read_list_ref->{$read_id}->{mdust_count1};
    my $mdust_count2 = $read_list_ref->{$read_id}->{mdust_count2};
    my $read1_n_count = $read_list_ref->{$read_id}->{ambig1};
    my $read2_n_count = $read_list_ref->{$read_id}->{ambig2};

    #1.) Low complexity check
    #Count those reads with more than $low_complexity_cutoff mdust bases that would not already be eliminated because of ambiguous bases

    #2.) Determine the subset of low complexity sequences that are actually polyAT sequences
    #Look for polyA/T reads by searching for reads with a consecutive stretch of 12 A's or T's and at least $low_complexity_cutoff A's or T's overall

    if ($mdust_count1 > $low_complexity_cutoff && $read1_n_count <= 1){
      $R1_low_complexity_reads++;
      $read_list_ref->{$read_id}->{status1} = "Low_Complexity";

      my $R1_polyA_check = 0;
      my $R1_polyT_check = 0;
      if ($read1 =~ /AAAAAAAAAAAA/){
	my $A_count = ($read1 =~ tr/A/A/);
	if ($A_count > $low_complexity_cutoff){
	  $R1_polyA_check = 1;
	}
      }
      if ($read1 =~ /TTTTTTTTTTTT/){
	my $T_count = ($read1 =~ tr/T/T/);
	if ($T_count > $low_complexity_cutoff){
	  $R1_polyT_check = 1;
	}
      }
      if ($R1_polyA_check == 1 || $R1_polyT_check == 1){
	$R1_polyAT_reads++;
      }

    }
    if ($mdust_count2 > $low_complexity_cutoff && $read2_n_count <= 1){
      $R2_low_complexity_reads++;
      $read_list_ref->{$read_id}->{status2} = "Low_Complexity";

      my $R2_polyA_check = 0;
      my $R2_polyT_check = 0;
      if ($read2 =~ /AAAAAAAAAAAA/){
	my $A_count = ($read2 =~ tr/A/A/);
	if ($A_count > $low_complexity_cutoff){
	  $R2_polyA_check = 1;
	}
      }
      if ($read2 =~ /TTTTTTTTTTTT/){
	my $T_count = ($read2 =~ tr/T/T/);
	if ($T_count > $low_complexity_cutoff){
	  $R2_polyT_check = 1;
	}
      }
      if ($R2_polyA_check == 1 || $R2_polyT_check == 1){
	$R2_polyAT_reads++;
      }
    }

    $R1_mdust_bases += ($read_list_ref->{$read_id}->{mdust_count1} - $read_list_ref->{$read_id}->{ambig1});
    $R2_mdust_bases += ($read_list_ref->{$read_id}->{mdust_count2} - $read_list_ref->{$read_id}->{ambig2});


    #3.) Duplicate reads of a pair check
    #Count the cases where R1 and R2 are identical - will need to watch out for these later, for now just note them to see if a library has an unusual number
    #Do not count those that have already been defined as low quality reads or low complexity reads
    my $original_read2 = $read2;
    $original_read2 = reverse $original_read2;
    $original_read2 =~ tr/gatcnGATCN/ctagnCTAGN/;
    my $duplicate_found = 0;

    if (($read1 eq $original_read2) && ($read1_n_count <= 1) && ($read2_n_count <= 1) && ($mdust_count1 <= $low_complexity_cutoff) && ($mdust_count2 <= $low_complexity_cutoff)){
      $duplicate_reads++;
      $duplicate_found = 1;
      $read_list_ref->{$read_id}->{status1} = "Duplicate";
      $read_list_ref->{$read_id}->{status2} = "Duplicate";
      #print YELLOW, "\nBEFORE reverse complement\n\tR1: $read1\n\tR2: $original_read2\n", RESET;
    }
    if (($read1 eq $read2) && ($read1_n_count <= 1) && ($read2_n_count <= 1) && ($mdust_count1 <= $low_complexity_cutoff) && ($mdust_count2 <= $low_complexity_cutoff)){
      $duplicate_reads_rc++;
      $duplicate_found = 1;
      $read_list_ref->{$read_id}->{status1} = "Duplicate";
      $read_list_ref->{$read_id}->{status2} = "Duplicate";
      #print YELLOW, "\nAFTER reverse complement\n\tR1: $read1\n\tR2: $read2\n", RESET;
    }

    #4.) Single N check
    #Note the number of reads with a single N.  These will be allowed but will be counted for interest's sake
    #Only count those that are not failed for some reason (low complexity or duplicate)
    if ($read1_n_count == 1 && $mdust_count1 <= $low_complexity_cutoff && $duplicate_found == 0){
      $R1_singleN_reads++;
    }
    if ($read2_n_count == 1 && $mdust_count2 <= $low_complexity_cutoff && $duplicate_found == 0){
      $R2_singleN_reads++;
    }

    print RESULTS_FILE "$read_id\t$read_list_ref->{$read_id}->{id1}\t$read_list_ref->{$read_id}->{id2}\t$read_list_ref->{$read_id}->{status1}\t$read_list_ref->{$read_id}->{status2}\t$read_list_ref->{$read_id}->{sequence1}\t$read_list_ref->{$read_id}->{sequence2}\t$read_list_ref->{$read_id}->{ambig1}\t$read_list_ref->{$read_id}->{ambig2}\t$read_list_ref->{$read_id}->{mdust_count1}\t$read_list_ref->{$read_id}->{mdust_count2}\n";

  }
  return();
}


###################################################################################################
#Print out statistics summarizing this lane
###################################################################################################
sub printStats{
  my %args = @_;
  my $flowcell_name = $args{'-flowcell_name'};
  my $lane_number = $args{'-lane_number'};
  my $read_record_dir = $args{'-read_record_dir'};

  my $duplicate_reads_p = sprintf("%.2f", (($duplicate_reads/$total_readpairs)*100));
  my $duplicate_reads_rc_p = sprintf("%.2f", (($duplicate_reads_rc/$total_readpairs)*100));

  my $R12_low_quality_reads = $R1_low_quality_reads + $R2_low_quality_reads;
  my $R12_low_quality_reads_p = sprintf("%.2f", (($R12_low_quality_reads/$total_reads)*100));
  my $R1_low_quality_reads_p = sprintf("%.2f", (($R1_low_quality_reads/$total_reads)*100));
  my $R2_low_quality_reads_p = sprintf("%.2f", (($R2_low_quality_reads/$total_reads)*100));

  my $R12_low_phred_reads = $R1_low_phred_reads + $R2_low_phred_reads;
  my $R12_low_phred_reads_p = sprintf("%.2f", (($R12_low_phred_reads/$total_reads)*100));
  my $R1_low_phred_reads_p = sprintf("%.2f", (($R1_low_phred_reads/$total_reads)*100));
  my $R2_low_phred_reads_p = sprintf("%.2f", (($R2_low_phred_reads/$total_reads)*100));

  my $R12_grand_total_phred = $R1_grand_total_phred + $R2_grand_total_phred;
  my $R1_grand_mean_phred = sprintf("%.2f", ($R1_grand_total_phred/$total_readpairs));
  my $R2_grand_mean_phred = sprintf("%.2f", ($R2_grand_total_phred/$total_readpairs));
  my $R12_grand_mean_phred = sprintf("%.2f", ($R12_grand_total_phred/$total_reads));

  my $R12_N_filtered_reads = $R1_N_filtered_reads + $R2_N_filtered_reads;
  my $R12_N_filtered_reads_p = sprintf("%.2f", (($R12_N_filtered_reads/$total_reads)*100));
  my $R1_N_filtered_reads_p = sprintf("%.2f", (($R1_N_filtered_reads/$total_reads)*100));
  my $R2_N_filtered_reads_p = sprintf("%.2f", (($R2_N_filtered_reads/$total_reads)*100));

  my $R12_singleN_reads = $R1_singleN_reads + $R2_singleN_reads;
  my $R12_singleN_reads_p = sprintf("%.2f", (($R12_singleN_reads/$total_reads)*100));
  my $R1_singleN_reads_p = sprintf("%.2f", (($R1_singleN_reads/$total_reads)*100));
  my $R2_singleN_reads_p = sprintf("%.2f", (($R2_singleN_reads/$total_reads)*100));

  my $R12_low_complexity_reads = $R1_low_complexity_reads + $R2_low_complexity_reads;
  my $R12_low_complexity_reads_p = sprintf("%.2f", (($R12_low_complexity_reads/$total_reads)*100));
  my $R1_low_complexity_reads_p = sprintf("%.2f", (($R1_low_complexity_reads/$total_reads)*100));
  my $R2_low_complexity_reads_p = sprintf("%.2f", (($R2_low_complexity_reads/$total_reads)*100));

  my $R12_polyAT_reads = $R1_polyAT_reads + $R2_polyAT_reads;
  my $R12_polyAT_reads_p = sprintf("%.2f", (($R12_polyAT_reads/$total_reads)*100));
  my $R1_polyAT_reads_p = sprintf("%.2f", (($R1_polyAT_reads/$total_reads)*100));
  my $R2_polyAT_reads_p = sprintf("%.2f", (($R2_polyAT_reads/$total_reads)*100));

  my $total_readpairs_m = sprintf("%.2f", ($total_readpairs/1000000));
  my $total_reads_m = sprintf("%.2f", ($total_reads/1000000));

  my $total_bases_m = sprintf("%.2f", ($total_bases/1000000));

  my $R12_ambiguous_bases = $R1_ambiguous_bases + $R2_ambiguous_bases;
  my $R12_ambiguous_bases_p = sprintf("%.2f", (($R12_ambiguous_bases/$total_bases)*100));
  my $R1_ambiguous_bases_p = sprintf("%.2f", (($R1_ambiguous_bases/$total_bases)*100));
  my $R2_ambiguous_bases_p = sprintf("%.2f", (($R2_ambiguous_bases/$total_bases)*100));

  my $R12_mdust_bases = $R1_mdust_bases + $R2_mdust_bases;
  my $R12_mdust_bases_p = sprintf("%.2f", (($R12_mdust_bases/$total_bases)*100));
  my $R1_mdust_bases_p = sprintf("%.2f", (($R1_mdust_bases/$total_bases)*100));
  my $R2_mdust_bases_p = sprintf("%.2f", (($R2_mdust_bases/$total_bases)*100));

  my $R12_A_count_p = sprintf("%.3f", ((($read1_a_count + $read2_a_count)/$total_bases)*100));
  my $R12_C_count_p = sprintf("%.3f", ((($read1_c_count + $read2_c_count)/$total_bases)*100));
  my $R12_G_count_p = sprintf("%.3f", ((($read1_g_count + $read2_g_count)/$total_bases)*100));
  my $R12_T_count_p = sprintf("%.3f", ((($read1_t_count + $read2_t_count)/$total_bases)*100));

  my $R1_A_count_p = sprintf("%.3f", ((($read1_a_count)/($total_bases/2))*100));
  my $R1_C_count_p = sprintf("%.3f", ((($read1_c_count)/($total_bases/2))*100));
  my $R1_G_count_p = sprintf("%.3f", ((($read1_g_count)/($total_bases/2))*100));
  my $R1_T_count_p = sprintf("%.3f", ((($read1_t_count)/($total_bases/2))*100));

  my $R2_A_count_p = sprintf("%.3f", ((($read2_a_count)/($total_bases/2))*100));
  my $R2_C_count_p = sprintf("%.3f", ((($read2_c_count)/($total_bases/2))*100));
  my $R2_G_count_p = sprintf("%.3f", ((($read2_g_count)/($total_bases/2))*100));
  my $R2_T_count_p = sprintf("%.3f", ((($read2_t_count)/($total_bases/2))*100));


  print BLUE, "\n\nThe following are statistics describing this flowcell: $flowcell_name Lane: $lane_number\n", RESET;
  print BLUE, "\n\tTotal ReadPairs = $total_readpairs ($total_readpairs_m m)", RESET;
  print BLUE, "\n\t\tNumber of ReadPairs where both reads are identical (and each read has no more than 1 ambiguous base) = $duplicate_reads ($duplicate_reads_p%)", RESET;
  print BLUE, "\n\t\tNumber of ReadPairs where both reads are identical AFTER reverse complementing Read2 = $duplicate_reads_rc ($duplicate_reads_rc_p%)", RESET;
  print BLUE, "\n\n\tTotal Individual Reads = $total_reads ($total_reads_m m)", RESET;

  print BLUE, "\n\n\tLow Quality Reads (> $n_cutoff ambiguous bases) = $R12_N_filtered_reads ($R12_N_filtered_reads_p%)", RESET;
  print BLUE, "\n\t\tRead1 = $R1_N_filtered_reads ($R1_N_filtered_reads_p%)\tRead2 = $R2_N_filtered_reads ($R2_N_filtered_reads_p%)", RESET;

  print BLUE, "\n\n\tLow Quality Reads (< $phred_cutoff phred score) = $R12_low_phred_reads ($R12_low_phred_reads_p%)", RESET;
  print BLUE, "\n\t\tRead1 = $R1_low_phred_reads ($R1_low_phred_reads_p%)\tRead2 = $R2_low_phred_reads ($R2_low_phred_reads_p%)", RESET;

  print BLUE, "\n\n\tGrand Mean phred score (all reads) = $R12_grand_mean_phred", RESET;
  print BLUE, "\n\t\tRead1 = $R1_grand_mean_phred\tRead2 = $R2_grand_mean_phred", RESET;

  print BLUE, "\n\n\tLow Quality Reads (> $n_cutoff ambiguous bases OR < $phred_cutoff phred score) = $R12_low_quality_reads ($R12_low_quality_reads_p%)", RESET;
  print BLUE, "\n\t\tRead1 = $R1_low_quality_reads ($R1_low_quality_reads_p%)\tRead2 = $R2_low_quality_reads ($R2_low_quality_reads_p%)", RESET;

  print BLUE, "\n\n\tSingle N Reads (ambiguous bases == 1) = $R12_singleN_reads ($R12_singleN_reads_p%)", RESET;
  print BLUE, "\n\t\tRead1 = $R1_singleN_reads ($R1_singleN_reads_p%)\tRead2 = $R2_singleN_reads ($R2_singleN_reads_p%)", RESET;

  print BLUE, "\n\n\tLow Complexity Reads (> $low_complexity_cutoff mdust masked bases & <= 1 ambiguous bases) = $R12_low_complexity_reads ($R12_low_complexity_reads_p%)", RESET;
  print BLUE, "\n\t\tRead1 = $R1_low_complexity_reads ($R1_low_complexity_reads_p%)\tRead2 = $R2_low_complexity_reads ($R2_low_complexity_reads_p%)", RESET;

  print BLUE, "\n\n\tPolyAT Reads (> $low_complexity_cutoff A's or T's and a stretch of 12 consecutive) = $R12_polyAT_reads ($R12_polyAT_reads_p%)", RESET;
  print BLUE, "\n\t\tRead1 = $R1_polyAT_reads ($R1_polyAT_reads_p%)\tRead2 = $R2_polyAT_reads ($R2_polyAT_reads_p%)", RESET;

  print BLUE, "\n\n\tTotal Bases generated = $total_bases ($total_bases_m mB)", RESET;
  print BLUE, "\n\n\tAmbiguous Bases = $R12_ambiguous_bases ($R12_ambiguous_bases_p%)", RESET;
  print BLUE, "\n\t\tRead1 = $R1_ambiguous_bases ($R1_ambiguous_bases_p%)\tRead2 = $R2_ambiguous_bases ($R2_ambiguous_bases_p%)", RESET;

  print BLUE, "\n\n\tLow Complexity Bases = $R12_mdust_bases ($R12_mdust_bases_p%)", RESET;
  print BLUE, "\n\t\tRead1 = $R1_mdust_bases ($R1_mdust_bases_p%)\tRead2 = $R2_mdust_bases ($R2_mdust_bases_p%)", RESET;

  print BLUE, "\n\n\tPercent base content:", RESET;
  print BLUE, "\n\t\tBoth Reads:\tA = $R12_A_count_p%\tC = $R12_C_count_p%\tG = $R12_G_count_p%\tT = $R12_T_count_p%", RESET;
  print BLUE, "\n\t\tRead1:\tA = $R1_A_count_p%\tC = $R1_C_count_p%\tG = $R1_G_count_p%\tT = $R1_T_count_p%", RESET;
  print BLUE, "\n\t\tRead2:\tA = $R2_A_count_p%\tC = $R2_C_count_p%\tG = $R2_G_count_p%\tT = $R2_T_count_p%\n\n", RESET;

  print LOG "\n\nThe following are statistics describing this flowcell: $flowcell_name Lane: $lane_number\n";
  print LOG "\n\tTotal ReadPairs = $total_readpairs ($total_readpairs_m m)";
  print LOG "\n\t\tNumber of ReadPairs where both reads are identical (and each read has no more than 1 ambiguous base) = $duplicate_reads ($duplicate_reads_p%)";
  print LOG "\n\t\tNumber of ReadPairs where both reads are identical AFTER reverse complementing Read2 = $duplicate_reads_rc ($duplicate_reads_rc_p%)";
  print LOG "\n\n\tTotal Individual Reads = $total_reads ($total_reads_m m)";

  print LOG "\n\n\tLow Quality Reads (> $n_cutoff ambiguous bases) = $R12_N_filtered_reads ($R12_N_filtered_reads_p%)";
  print LOG "\n\t\tRead1 = $R1_N_filtered_reads ($R1_N_filtered_reads_p%)\tRead2 = $R2_N_filtered_reads ($R2_N_filtered_reads_p%)";

  print LOG "\n\n\tLow Quality Reads (< $phred_cutoff phred score) = $R12_low_phred_reads ($R12_low_phred_reads_p%)";
  print LOG "\n\t\tRead1 = $R1_low_phred_reads ($R1_low_phred_reads_p%)\tRead2 = $R2_low_phred_reads ($R2_low_phred_reads_p%)";

  print LOG "\n\n\tGrand Mean phred score (all reads) = $R12_grand_mean_phred";
  print LOG "\n\t\tRead1 = $R1_grand_mean_phred\tRead2 = $R2_grand_mean_phred";

  print LOG "\n\n\tLow Quality Reads (> $n_cutoff ambiguous bases OR < $phred_cutoff phred score) = $R12_low_quality_reads ($R12_low_quality_reads_p%)";
  print LOG "\n\t\tRead1 = $R1_low_quality_reads ($R1_low_quality_reads_p%)\tRead2 = $R2_low_quality_reads ($R2_low_quality_reads_p%)";

  print LOG "\n\n\tSingle N Reads (ambiguous bases == 1) = $R12_singleN_reads ($R12_singleN_reads_p%)";
  print LOG "\n\t\tRead1 = $R1_singleN_reads ($R1_singleN_reads_p%)\tRead2 = $R2_singleN_reads ($R2_singleN_reads_p%)";

  print LOG "\n\n\tLow Complexity Reads (> $low_complexity_cutoff mdust masked bases & <= 1 ambiguous bases) = $R12_low_complexity_reads ($R12_low_complexity_reads_p%)";
  print LOG "\n\t\tRead1 = $R1_low_complexity_reads ($R1_low_complexity_reads_p%)\tRead2 = $R2_low_complexity_reads ($R2_low_complexity_reads_p%)";

  print LOG "\n\n\tPolyAT Reads (> $low_complexity_cutoff A's or T's and a stretch of 12 consecutive) = $R12_polyAT_reads ($R12_polyAT_reads_p%)";
  print LOG "\n\t\tRead1 = $R1_polyAT_reads ($R1_polyAT_reads_p%)\tRead2 = $R2_polyAT_reads ($R2_polyAT_reads_p%)";

  print LOG "\n\n\tTotal Bases generated = $total_bases ($total_bases_m mB)";
  print LOG "\n\n\tAmbiguous Bases = $R12_ambiguous_bases ($R12_ambiguous_bases_p%)";
  print LOG "\n\t\tRead1 = $R1_ambiguous_bases ($R1_ambiguous_bases_p%)\tRead2 = $R2_ambiguous_bases ($R2_ambiguous_bases_p%)";

  print LOG "\n\n\tLow Complexity Bases = $R12_mdust_bases ($R12_mdust_bases_p%)";
  print LOG "\n\t\tRead1 = $R1_mdust_bases ($R1_mdust_bases_p%)\tRead2 = $R2_mdust_bases ($R2_mdust_bases_p%)";

  print LOG "\n\n\tPercent base content:";
  print LOG "\n\t\tBoth Reads:\tA = $R12_A_count_p%\tC = $R12_C_count_p%\tG = $R12_G_count_p%\tT = $R12_T_count_p%";
  print LOG "\n\t\tRead1:\tA = $R1_A_count_p%\tC = $R1_C_count_p%\tG = $R1_G_count_p%\tT = $R1_T_count_p%";
  print LOG "\n\t\tRead2:\tA = $R2_A_count_p%\tC = $R2_C_count_p%\tG = $R2_G_count_p%\tT = $R2_T_count_p%\n\n";

  #Also dump these stats to a tab delimited file stored in the read_records_dir
  my $summary_stats_file = "$read_record_dir"."Summary/SummaryStats.csv";

  unless (-e $summary_stats_file && -s $summary_stats_file){
    my $cmd1 = "echo 'flowcell_name,lane_number,total_readpairs,duplicate_readpair_reads,duplicate_readpair_reads_afterRC,total_reads,R12_low_quality_reads,R1_low_quality_reads,R2_low_quality_reads,R12_singleN_reads,R1_singleN_reads,R2_singleN_reads,R12_low_complexity_reads,R1_low_complexity_reads,R2_low_complexity_reads,R12_polyAT_reads,R1_polyAT_reads,R2_polyAT_reads,total_bases,R12_ambiguous_bases,R1_ambiguous_bases,R2_ambiguous_bases,R12_mdust_bases,R1_mdust_bases,R2_mdust_bases,R12_A_percent,R12_C_percent,R12_G_percent,R12_T_percent,R1_A_percent,R1_C_percent,R1_G_percent,R1_T_percent,R2_A_percent,R2_C_percent,R2_G_percent,R2_T_percent,R12_N_filtered_reads,R1_N_filtered_reads,R2_N_filtered_reads,R12_low_phred_reads,R1_low_phred_reads,R2_low_phred_reads,R12_grand_mean_phred,R1_grand_mean_phred,R2_grand_mean_phred' >> $summary_stats_file";
    system($cmd1);
  }
  my $cmd2 = "echo '$flowcell_name,$lane_number,$total_readpairs,$duplicate_reads,$duplicate_reads_rc,$total_reads,$R12_low_quality_reads,$R1_low_quality_reads,$R2_low_quality_reads,$R12_singleN_reads,$R1_singleN_reads,$R2_singleN_reads,$R12_low_complexity_reads,$R1_low_complexity_reads,$R2_low_complexity_reads,$R12_polyAT_reads,$R1_polyAT_reads,$R2_polyAT_reads,$total_bases,$R12_ambiguous_bases,$R1_ambiguous_bases,$R2_ambiguous_bases,$R12_mdust_bases,$R1_mdust_bases,$R2_mdust_bases,$R12_A_count_p,$R12_C_count_p,$R12_G_count_p,$R12_T_count_p,$R1_A_count_p,$R1_C_count_p,$R1_G_count_p,$R1_T_count_p,$R2_A_count_p,$R2_C_count_p,$R2_G_count_p,$R2_T_count_p,$R12_N_filtered_reads,$R1_N_filtered_reads,$R2_N_filtered_reads,$R12_low_phred_reads,$R1_low_phred_reads,$R2_low_phred_reads,$R12_grand_mean_phred,$R1_grand_mean_phred,$R2_grand_mean_phred' >> $summary_stats_file";

  system($cmd2);

  #Print a file with the phred score distribution
  my $phred_dist_file = "$read_record_dir"."Summary/"."$flowcell_name"."_"."$lane_number"."_"."phred.txt";
  open (PHRED, ">$phred_dist_file") || die "\nCould not open: $phred_dist_file\n\n";
  print PHRED "phred_score\tR1_count\tR2_count\n";
  foreach my $phred (sort{$a<=>$b} keys %phred_dist){
    #Some phred scores may not have instances for both R1 and R2 - check and set to 'NA' if no value found
    unless($phred_dist{$phred}{R1}){$phred_dist{$phred}{R1}="NA";}
    unless($phred_dist{$phred}{R2}){$phred_dist{$phred}{R2}="NA";}
    print PHRED "$phred\t$phred_dist{$phred}{R1}\t$phred_dist{$phred}{R2}\n";
  }
  close PHRED;

  return();
}

