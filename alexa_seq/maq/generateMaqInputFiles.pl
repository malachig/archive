#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#Given a source directory for a single lane of raw solexa data (_seq.txt and _prb.txt files) create the files neccessary to run MAQ

#This involves the following file conversions:
#Seq/Prb  -->  .fq (fastq)  -->  .bfq (binary fastq)

#a.) Concatenate *_seq.txt and *_prb.txt
#b.) Convert *_seq.txt and *_prb.txt files to a .fastq file 
#     - e.g. /home/malachig/tools/MAQ/maq-0.6.6_i686-linux/fq_all2std.pl seqprb2std lane1_seq.txt lane1_prb.txt > lane1.fastq
#c.) Convert .fastq file to MAQ .bfq (binary fastq) file
#     - e.g. /home/malachig/tools/MAQ/maq-0.6.6_i686-linux/maq fastq2bfq [-n nreads] lane1.fastq lane1.bfq

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);

#Initialize command line options
my $maq_bin_dir = '';
my $input_dir = '';
my $flowcell_name = '';
my $lane_number = '';
my $output_dir = '';
my $paired = '';
my $trim = '';
my $tile_count = '';
my $log_file = '';

GetOptions ('maq_bin_dir=s'=>\$maq_bin_dir, 'input_dir=s'=>\$input_dir, 'flowcell_name=s'=>\$flowcell_name, 'lane_number=i'=>\$lane_number,
	    'output_dir=s'=>\$output_dir, 'paired=s'=>\$paired, 'trim=i'=>\$trim, 'tile_count=i'=>\$tile_count,
	    'log_file=s'=>\$log_file);

#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the complete path to your MAQ installation directory using: --maq_bin_dir", RESET;
print GREEN, "\n\tSpecify the directory containing Solexa seq files using: --input_dir", RESET;
print GREEN, "\n\tSpecify the flowcell name for this input file file using: --flowcell_name (this will be appended to read names)", RESET;
print GREEN, "\n\tSpecify the lane number for this input file using: --lane_number (this will be verified against files in the input_dir and used for the log file)", RESET;
print GREEN, "\n\tSpecify an output directory for MAQ formated files using: --output_dir", RESET;
print GREEN, "\n\tIf the input data is paired, specify this so that R1 and R2 will be split by using: --paired=yes", RESET;
print GREEN, "\n\tIf you wish to trim the ends of each read of a pair use the trim option.  (e.g. --trim=1 will remove 1 base from the end of both reads)", RESET;
print GREEN, "\n\tSpecify the number of tiles per flow cell lane (usually 100 or 200) using: --tile_count", RESET;
print GREEN, "\n\tAlso specify the name of a log file to be used to store a record of this script using: --log_file", RESET;

print GREEN, "\n\nExample: generateMaqInputFiles.pl  --maq_bin_dir=/home/malachig/tools/MAQ/maq-0.6.6_x86_64-linux/  --input_dir=/archive/solexa1_2/data2/080325_SOLEXA4_0028_20821AAXX/HS04391..L1/Data/C1-74_Firecrest1.9.1_03-04-2008_aldente/Bustard1.9.1_03-04-2008_aldente/  --flowcell_name=20821AAXX  --lane_number=1  --output_dir=/projects/malachig/solexa/maq_analysis  --paired=yes  --tile_count=100  --log_file=/projects/malachig/solexa/logs/HS04391/20821AAXX_Lane1/generateMaqInputFiles_LOG.txt\n\n", RESET;


#Make sure all options were specified
unless ($maq_bin_dir && $input_dir && $flowcell_name && $lane_number && $tile_count && $output_dir && $log_file){
  print RED, "\nRequired input parameter(s) missing or incorrect format!\n\n", RESET;
  exit();
}
unless ($paired =~ /^yes|y/i){
  print YELLOW, "\nPaired option not specified ...  Assuming this data consists of single end reads!!!\n\n", RESET;
}

#open log file to record the output of this script
open (LOG, ">$log_file") || die "\nCould not open log_file: $log_file\n\n";
print LOG "\nFollowing is a log for the script generateMaqInputFiles.pl\n\n";
print LOG "\nUser specified the following options:\ninput_dir = $input_dir\nflowcell_name = $flowcell_name\noutput_dir = $output_dir\npaired = $paired\ntrim = $trim\nlog_file = $log_file\n";

if ($trim =~ /\d+/){
  print YELLOW, "\nBoth R1 and R2 will be trimmed by $trim bases\n\n", RESET;
  print LOG "\nBoth R1 and R2 will be trimmed by $trim bases\n\n";
}else{
  $trim = 0;
  print YELLOW, "\nNo --trim option specified, read length will remain unchanged\n\n", RESET;
  print LOG "\nNo --trim option specified, read length will remain unchanged\n\n";
}

#Add trailing '/' to directory paths if they were forgotten
unless ($input_dir =~ /.*\/$/){
  $input_dir = "$input_dir"."/";
}
unless ($output_dir =~ /.*\/$/){
  $output_dir = "$output_dir"."/";
}
unless ($maq_bin_dir =~ /.*\/$/){
  $maq_bin_dir = "$maq_bin_dir"."/";
}

#Check the validity of required directories and binaries specified by the user
unless (-d $input_dir && -e $input_dir){
  print RED, "\nInput directory does not appear to be a valid directory:\n\tinput_dir = $input_dir\n\n", RESET;
  exit();
}
unless (-d $output_dir && -e $output_dir){
  print RED, "\nOutput directory does not appear to be a valid directory:\n\toutput_dir = $output_dir\n\n", RESET;
  exit();
}
unless (-d $maq_bin_dir && -e $maq_bin_dir){
  print RED, "\nMAQ binary directory does not appear to be a valid directory:\n\tmaq_bin_dir = $maq_bin_dir\n\n", RESET;
  exit();
}

my $fq_all2std_bin = "$maq_bin_dir"."fq_all2std.pl";
my $maq_bin = "$maq_bin_dir"."maq";

#Check the input source directory to ensure that the analysis was complete and look for warnings - also determine the lane number
my $lane_number_check = &checkSourceDirectory('-input_dir'=>$input_dir, '-flowcell_name'=>$flowcell_name);

if ($lane_number_check==$lane_number){
  print BLUE, "\n\tThe lane number found matches the lane number specified at runtime", RESET;
  print BLUE, "\n\n", RESET;
  print LOG "\n\tThe lane number found matches the lane number specified at runtime";
  print LOG "\n\n";
}else{
  print RED, "\nSpecified lane number ($lane_number) does not match lane number found ($lane_number_check)\n\n", RESET;
  exit();
}


#a.) Concatenate *_seq.txt and *_prb.txt

#Check source directory for raw files and create a concatenated raw seq file and save this to the specified path
my @raw_seq_files = &createRawSeqFile('-input_dir'=>$input_dir, '-flowcell_name'=>$flowcell_name, '-raw_read_dir'=>$output_dir, '-lane_number'=>$lane_number);

#get the number of reads by performing a wordcount on the _seq.txt file
my $wc_cmd = "wc -l $raw_seq_files[0]";
my $result = `$wc_cmd`;
my $read_count;
if ($result =~ /^(\d+)/){
  $read_count = $1;
}else{
  print RED, "\nCould not figure out the number of reads in the _seq.txt file\n\n", RESET;
  exit();
}
print YELLOW, "\nFound $read_count reads in the raw _seq.txt file\n\n", RESET;

#Check source directory for raw files and create a concatenated raw prb file and save this to the specified path
my @raw_prb_files = &createRawPrbFile('-input_dir'=>$input_dir, '-flowcell_name'=>$flowcell_name, '-raw_read_dir'=>$output_dir, '-lane_number'=>$lane_number);


#b.) Convert *_seq.txt and *_prb.txt files to a .fastq file 
#     - e.g. /home/malachig/tools/MAQ/maq-0.6.6_i686-linux/fq_all2std.pl seqprb2std lane1_seq.txt lane1_prb.txt > lane1.fastq
my $fastq_file;
my $fastq_file_R1;
my $fastq_file_R2;

if ($paired =~ /^yes|y/i){

  $fastq_file_R1 = "$output_dir"."$flowcell_name"."_Lane"."$lane_number"."_R1".".fastq";
  my $fq_cmd_R1 = "$fq_all2std_bin seqprb2std $raw_seq_files[0] $raw_prb_files[0] > $fastq_file_R1";

  $fastq_file_R2 = "$output_dir"."$flowcell_name"."_Lane"."$lane_number"."_R2".".fastq";
  my $fq_cmd_R2 = "$fq_all2std_bin seqprb2std $raw_seq_files[1] $raw_prb_files[1] > $fastq_file_R2";

  print BLUE, "\nConverting R1 _seq.txt and _prb.txt files to a single FASTQ file with the following command\n\n$fq_cmd_R1\n", RESET;
  print LOG "\nConverting R1 _seq.txt and _prb.txt files to a single FASTQ file with the following command\n\n$fq_cmd_R1\n";
  system ($fq_cmd_R1);

  print BLUE, "\nConverting R2 _seq.txt and _prb.txt files to a single FASTQ file with the following command\n\n$fq_cmd_R2\n", RESET;
  print LOG "\nConverting R2 _seq.txt and _prb.txt files to a single FASTQ file with the following command\n\n$fq_cmd_R2\n";
  system ($fq_cmd_R2);

}else{
  $fastq_file = "$output_dir"."$flowcell_name"."_Lane"."$lane_number".".fastq";
  my $fq_cmd = "$fq_all2std_bin seqprb2std $raw_seq_files[0] $raw_prb_files[0] > $fastq_file";

  print BLUE, "\nConverting _seq.txt and _prb.txt files to a single FASTQ file with the following command\n\n$fq_cmd\n", RESET;
  print LOG "\nConverting _seq.txt and _prb.txt files to a single FASTQ file with the following command\n\n$fq_cmd\n";
  system ($fq_cmd);
}

#c.) Convert .fastq file to MAQ .bfq (binary fastq) file
#     - e.g. /home/malachig/tools/MAQ/maq-0.6.6_i686-linux/maq fastq2bfq [-n nreads] lane1.fastq lane1.bfq

if ($paired =~ /^yes|y/i){
  my $binary_fastq_file_R1 = "$output_dir"."$flowcell_name"."_Lane"."$lane_number"."_R1".".bfq";
  my $bfq_cmd_R1 = "$maq_bin fastq2bfq -n $read_count $fastq_file_R1 $binary_fastq_file_R1";

  my $binary_fastq_file_R2 = "$output_dir"."$flowcell_name"."_Lane"."$lane_number"."_R2".".bfq";
  my $bfq_cmd_R2 = "$maq_bin fastq2bfq -n $read_count $fastq_file_R2 $binary_fastq_file_R2";

  print BLUE, "\nConverting R1 FASTQ file to a R1 Binary FASTQ file with the following command\n\n$bfq_cmd_R1\n", RESET;
  print LOG "\nConverting R1 FASTQ file to a R1 Binary FASTQ file with the following command\n\n$bfq_cmd_R1\n";
  system ($bfq_cmd_R1);

  print BLUE, "\nConverting R2 FASTQ file to a R2 Binary FASTQ file with the following command\n\n$bfq_cmd_R2\n", RESET;
  print LOG "\nConverting R2 FASTQ file to a R2 Binary FASTQ file with the following command\n\n$bfq_cmd_R2\n";
  system ($bfq_cmd_R2);

  #Clean up R1 intermediary files
  my $del_cmd_seq_R1 = "rm -f $raw_seq_files[0]";
  my $del_cmd_prb_R1 = "rm -f $raw_prb_files[0]";
  #my $del_cmd_fq_R1 = "rm -f $fastq_file_R1";  #Leave the fastq file since it is used with Novoalign

  print BLUE, "\nDeleting intermediary R1 files with the following commands\n\n\t$del_cmd_seq_R1\n\t$del_cmd_prb_R1\n\n", RESET;
  print LOG "\nDeleting intermediary R1 files with the following commands\n\n\t$del_cmd_seq_R1\n\t$del_cmd_prb_R1\n\n";
  system ($del_cmd_seq_R1);
  system ($del_cmd_prb_R1);
  #system ($del_cmd_fq_R1);

  #Clean up R2 intermediary files
  my $del_cmd_seq_R2 = "rm -f $raw_seq_files[1]";
  my $del_cmd_prb_R2 = "rm -f $raw_prb_files[1]";
  #my $del_cmd_fq_R2 = "rm -f $fastq_file_R2";  #Leave the fastq file since it is used with Novoalign

  print BLUE, "\nDeleting intermediary R2 files with the following commands\n\n\t$del_cmd_seq_R2\n\t$del_cmd_prb_R2\n\n", RESET;
  print LOG "\nDeleting intermediary R2 files with the following commands\n\n\t$del_cmd_seq_R2\n\t$del_cmd_prb_R2\n\n";
  system ($del_cmd_seq_R2);
  system ($del_cmd_prb_R2);
  #system ($del_cmd_fq_R2);

}else{
  my $binary_fastq_file = "$output_dir"."$flowcell_name"."_Lane"."$lane_number".".bfq";
  my $bfq_cmd = "$maq_bin fastq2bfq -n $read_count $fastq_file $binary_fastq_file";

  print BLUE, "\nConverting FASTQ file to a Binary FASTQ file with the following command\n\n$bfq_cmd\n", RESET;
  print LOG "\nConverting FASTQ file to a Binary FASTQ file with the following command\n\n$bfq_cmd\n";
  system ($bfq_cmd);

  #Clean up intermediary files
  my $del_cmd_seq = "rm -f $raw_seq_files[0]";
  my $del_cmd_prb = "rm -f $raw_prb_files[0]";
  #my $del_cmd_fq = "rm -f $fastq_file";  #Leave the fastq file since it is used with Novoalign

  print BLUE, "\nDeleting intermediary files with the following commands\n\n\t$del_cmd_seq\n\t$del_cmd_prb\n\n", RESET;
  print LOG "\nDeleting intermediary files with the following commands\n\n\t$del_cmd_seq\n\t$del_cmd_prb\n\n";
  system ($del_cmd_seq);
  system ($del_cmd_prb);
  #system ($del_cmd_fq);
}

exit();



############################################################################################################
#Perform basic checks on the input source directory
############################################################################################################
sub checkSourceDirectory{
  my %args = @_;
  my $input_dir = $args{'-input_dir'};
  my $flowcell_name = $args{'-flowcell_name'};

  my $lane_found;

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

  #See if the specified input directory seems to correspond to the specified flowcell name
  if ($input_dir =~ /$flowcell_name/){
    print BLUE, "\n\n\tThe specified flowcell name ($flowcell_name) appears to match the specified source directory", RESET;
    print LOG "\n\n\tThe specified flowcell name ($flowcell_name) appears to match the specified source directory";
  }else{
    print RED, "\n\n\tThe specified flowcell name ($flowcell_name) doesn't seem to correspond to the specified source directory\n\n", RESET;
    exit();
  }

  #Check for '_finished.txt' and 'finished.txt' files which indicate that the analysis run completed
  my $cmd_a = "ls $input_dir"."finished.txt";
  my $cmd_b = "ls $input_dir"."s_"."$lane_number"."_finished.txt";
  my $result_a = `$cmd_a`;
  my $result_b = `$cmd_b`;

  my @finished_files_a = split(" ", $result_a);
  my @finished_files_b = split(" ", $result_b);

  if ((scalar(@finished_files_a) + scalar(@finished_files_b)) == 2){
    print BLUE, "\n\tFound two finished.txt files indicating that the analysis of this lane was completed", RESET;
    print LOG "\n\tFound two finished.txt files indicating that the analysis of this lane was completed";
  }else{
    print RED, "\n\tDid not find *finished.txt files indicating that the analysis of this lane may not be complete\n\tCheck with LIMS before proceeding!\n\n", RESET;
    exit();
  }

  opendir(DIRHANDLE, "$input_dir") || die "\nCannot open directory: $input_dir\n\n";
  my @files = readdir(DIRHANDLE);
  my @seq_files;

  foreach my $file (@files){

    #Skip all files except _seq.txt files
    my $lane_found;
    if ($file =~ /^s_(\d+)_\d+_seq\.txt$/ || $file =~ /^s_(\d+)_\d+_seq\.txt\.bz2$/){
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

  #Check to make sure seq files count matches that expected by user
  my $tile_files_found = scalar(@seq_files);
  unless ($tile_files_found == $tile_count){
    print RED, "\nFound file count ($tile_files_found) other than $tile_count seq files ... something is wrong?\n\n", RESET;
    exit();
  }
  closedir(DIRHANDLE);

  #Determine the lane number for the specified source directory
  my $test_file = $seq_files[0];
  if ($test_file =~ /^s_(\d+)_\d+_seq.txt/){
    $lane_found = $1;
  }else{
    print RED, "\nseq file format not understood - could not determine lane number for this source directory\n\n", RESET;
    exit();
  }

  unless ($lane_found == $lane_number){
    print RED, "\nLane number found does not match that supplied by user\n\n", RESET;
    exit();
  }

  print BLUE, "\n\tThe lane number for this source directory appears to be $lane_number", RESET;
  print BLUE, "\n\n", RESET;
  print LOG "\n\tThe lane number for this source directory appears to be $lane_number";
  print LOG "\n\n";

  return($lane_number);
}


############################################################################################################
#Create a raw Seq file by concatenating raw seq files in the source directory
############################################################################################################
sub createRawSeqFile{
  my %args = @_;
  my $input_dir = $args{'-input_dir'};
  my $flowcell_name = $args{'-flowcell_name'};
  my $raw_seq_dir = $args{'-raw_read_dir'};
  my $lane_number = $args{'-lane_number'};

  my @return_files;
  my $raw_seq_file;
  my $raw_seq_file_R1;
  my $raw_seq_file_R2;

  my %seq_files;

  #Figure out file names depending on the paired or unpaired cases
  if ($paired =~ /^yes|y/i){

    $raw_seq_file_R1 = "$raw_seq_dir"."$flowcell_name"."_Lane"."$lane_number"."_R1"."_seq.txt";
    $raw_seq_file_R2 = "$raw_seq_dir"."$flowcell_name"."_Lane"."$lane_number"."_R2"."_seq.txt";
    push (@return_files, $raw_seq_file_R1);
    push (@return_files, $raw_seq_file_R2);

    &checkFile('-file'=>$raw_seq_file_R1);
    &checkFile('-file'=>$raw_seq_file_R2);

  }else{

    $raw_seq_file = "$raw_seq_dir"."$flowcell_name"."_Lane"."$lane_number"."_seq.txt";
    push (@return_files, $raw_seq_file);

    print BLUE, "\nCreating raw_seq_file: $raw_seq_file\n\n", RESET;
    print LOG "\nCreating raw_seq_file: $raw_seq_file\n\n";

    &checkFile('-file'=>$raw_seq_file);
  }

  #Get the actual seq files for processing
  opendir(DIRHANDLE, "$input_dir") || die "\nCannot open directory: $input_dir\n\n";
  my @files = readdir(DIRHANDLE);

  foreach my $file (@files){

    #Skip all files except _seq.txt files that correspond to the specified lane
    my $lane_found;
    my $tile_found;

    if ($file =~ /^s_(\d+)_(\d+)_seq\.txt$/ || $file =~ /^s_(\d+)_(\d+)_seq\.txt\.bz2$/){
      $lane_found = $1;
      $tile_found = $2;
    }else{
      next();
    }

    #skip files corresponding to a lane other than that specified by the user (for directories where all lanes are together)
    unless ($lane_found == $lane_number){
      next();
    }
    $seq_files{$tile_found}{file} = "$input_dir"."$file";
    $seq_files{$tile_found}{name} = $file;
  }

  if ($paired =~ /^yes|y/i){
    print BLUE, "\nCreating concatenated raw sequence files:\n\t$raw_seq_file_R1\n\t$raw_seq_file_R2\n\n", RESET;
    print LOG "\nCreating concatenated raw sequence files:\n\t$raw_seq_file_R1\n\t$raw_seq_file_R2\n\n";

    foreach my $tile (sort {$a cmp $b} keys %seq_files){
      my $file = $seq_files{$tile}{file};
      my $file_name = $seq_files{$tile}{name};


      #If file is a _seq.txt.bz2 file, copy it to the working dir, bunzip2 it, cat it to the raw_read_file, and then delete it
      my $bzipped = 0;
      if ($file =~ /_seq\.txt\.bz2$/){
	$bzipped = 1;
	my $cp_cmd = "cp $file $raw_seq_dir";
	my $bunzip2_cmd = "bunzip2 "."$raw_seq_dir"."$file_name";
	system($cp_cmd);
	system($bunzip2_cmd);
	$file = "$raw_seq_dir"."$file_name";
	$file=~s/\.bz2$//; #the bunzip2'd file will have a new name
      }

      print YELLOW, "\n\tPrint: $file", RESET;

      open (SEQ, "$file") || die "\nCould not open SEQ file: $file\n\n";
      open (OUT_R1, ">>$raw_seq_file_R1") || die "\nCould not open R1 output seq file: $raw_seq_file_R1";
      open (OUT_R2, ">>$raw_seq_file_R2") || die "\nCould not open R1 output seq file: $raw_seq_file_R2";

      while (<SEQ>){
	my $line = $_;
	chomp($line);
	my @data = split("\t", $line);
	my $full_seq = $data[4];

	my $seq_length = length($full_seq);
	my $read_length = ($seq_length/2);

	my $read1 = substr($full_seq, 0, $read_length-$trim);
	my $read2 = substr($full_seq, $read_length, $read_length-$trim);

	print OUT_R1 "$data[0]\t$data[1]\t$data[2]\t$data[3]\t$read1\n";
	print OUT_R2 "$data[0]\t$data[1]\t$data[2]\t$data[3]\t$read2\n";

      }

      close (SEQ);
      close (OUT_R1);
      close (OUT_R2);

      #If the file was bzipped remove the temp file created
      if ($bzipped == 1){
	my $rm_cmd = "rm $file";
	system($rm_cmd);
      }

    }

  }else{
    print BLUE, "\nCreating concatenated raw sequence file:\n\t$raw_seq_file\n\n", RESET;
    print LOG "\nCreating concatenated raw sequence file:\n\t$raw_seq_file\n\n";

    foreach my $tile (sort {$a cmp $b} keys %seq_files){
      my $file = $seq_files{$tile}{file};
      my $file_name = $seq_files{$tile}{name};

      print YELLOW, "\n\tCat: $file", RESET;

      #If file is a _seq.txt.bz2 file, copy it to the working dir, bunzip2 it, cat it to the raw_read_file, and then delete it
      my $bzipped = 0;
      if ($file =~ /_seq\.txt\.bz2$/){
	$bzipped = 1;
	my $cp_cmd = "cp $file $raw_seq_dir";
	my $bunzip2_cmd = "bunzip2 "."$raw_seq_dir"."$file_name";
	system($cp_cmd);
	system($bunzip2_cmd);
	$file = "$raw_seq_dir"."$file_name";
	$file=~s/\.bz2$//; #the bunzip2'd file will have a new name
      }

      my $cmd = "cat $input_dir"."$file >> $raw_seq_file";
      system ($cmd);

      #If the file was bzipped remove the temp file created
      if ($bzipped == 1){
	my $rm_cmd = "rm $file";
	system($rm_cmd);
      }

    }
  }

  print "\n\n";

  return(@return_files);
}


############################################################################################################
#Create a raw Prb file by concatenating raw prb files in the source directory
############################################################################################################
sub createRawPrbFile{
  my %args = @_;
  my $input_dir = $args{'-input_dir'};
  my $flowcell_name = $args{'-flowcell_name'};
  my $raw_prb_dir = $args{'-raw_read_dir'};
  my $lane_number = $args{'-lane_number'};

  my %prb_files;

  my @return_files;
  my $raw_prb_file;
  my $raw_prb_file_R1;
  my $raw_prb_file_R2;

  #Figure out file names depending on the paired or unpaired cases
  if ($paired =~ /^yes|y/i){

    $raw_prb_file_R1 = "$raw_prb_dir"."$flowcell_name"."_Lane"."$lane_number"."_R1"."_prb.txt";
    $raw_prb_file_R2 = "$raw_prb_dir"."$flowcell_name"."_Lane"."$lane_number"."_R2"."_prb.txt";
    push (@return_files, $raw_prb_file_R1);
    push (@return_files, $raw_prb_file_R2);

    &checkFile('-file'=>$raw_prb_file_R1);
    &checkFile('-file'=>$raw_prb_file_R2);

  }else{

    $raw_prb_file = "$raw_prb_dir"."$flowcell_name"."_Lane"."$lane_number"."_prb.txt";
    push (@return_files, $raw_prb_file);

    print BLUE, "\nCreating raw_prb_file: $raw_prb_file\n\n", RESET;
    print LOG "\nCreating raw_prb_file: $raw_prb_file\n\n";

    &checkFile('-file'=>$raw_prb_file);
  }


  opendir(DIRHANDLE, "$input_dir") || die "\nCannot open directory: $input_dir\n\n";
  my @files = readdir(DIRHANDLE);
  my @prb_files;

  foreach my $file (@files){
    #Skip all files except _seq.txt files
    my $lane_found;
    my $tile_found;

    if ($file =~ /^s_(\d+)_(\d+)_prb\.txt$/ || $file =~ /^s_(\d+)_(\d+)_prb\.txt\.bz2$/){
      $lane_found = $1;
      $tile_found = $2;
    }else{
      next();
    }

    #skip files corresponding to a lane other than that specified by the user (for directories where all lanes are together)
    unless ($lane_found == $lane_number){
      next();
    }
    $prb_files{$tile_found}{file} = "$input_dir"."$file";
    $prb_files{$tile_found}{name} = "$file";

  }

  if ($paired =~ /^yes|y/i){
    print BLUE, "\nCreating concatenated raw prb files:\n\t$raw_prb_file_R1\n\t$raw_prb_file_R2\n\n", RESET;
    print LOG "\nCreating concatenated raw prb files:\n\t$raw_prb_file_R1\n\t$raw_prb_file_R2\n\n";

    foreach my $tile (sort {$a cmp $b} keys %prb_files){
      my $file = $prb_files{$tile}{file};
      my $file_name = $prb_files{$tile}{name};

      print YELLOW, "\n\tPrint: $file", RESET;

      #If file is a _prb.txt.bz2 file, copy it to the working dir, bunzip2 it, cat it to the raw_read_file, and then delete it
      my $bzipped = 0;
      if ($file =~ /_prb\.txt\.bz2$/){
	$bzipped = 1;
	my $cp_cmd = "cp $file $raw_prb_dir";
	my $bunzip2_cmd = "bunzip2 "."$raw_prb_dir"."$file_name";
	system($cp_cmd);
	system($bunzip2_cmd);
	$file = "$raw_prb_dir"."$file_name";
	$file=~s/\.bz2$//; #the bunzip2'd file will have a new name
      }

      open (PRB, "$file") || die "\nCould not open PRB file: $file\n\n";
      open (OUT_R1, ">>$raw_prb_file_R1") || die "\nCould not open R1 output prb file: $raw_prb_file_R1";
      open (OUT_R2, ">>$raw_prb_file_R2") || die "\nCould not open R1 output prb file: $raw_prb_file_R2";

      while (<PRB>){
	my $line = $_;
	chomp($line);


	my @prbs = split("\t", $line);
	my $seq_length = scalar(@prbs);
	my $read_length = ($seq_length/2);

	my @prbs_R1 = @prbs[0 .. (($read_length-1)-$trim)];
	my @prbs_R2 = @prbs[($read_length) .. (($seq_length-1)-$trim)];

	#print "\n\nDEBUG: 0 .. (($read_length-1)-$trim)";
	#print "\nDEBUG: ($read_length) .. (($seq_length-1)-$trim)\n";
	#exit();

	unless (scalar(@prbs_R1) == scalar(@prbs_R2)){
	  print RED, "\nProbability arrays not equal length!!\n\n", RESET;
	  exit();
	}

	my $temp_sep = $";
	$" = "\t";
	print OUT_R1 "@prbs_R1\n";
	print OUT_R2 "@prbs_R2\n";
	$" = $temp_sep;

      }

      close (PRB);
      close (OUT_R1);
      close (OUT_R2);

      #If the file was bzipped remove the temp file created
      if ($bzipped == 1){
	my $rm_cmd = "rm $file";
	system($rm_cmd);
      }

    }

  }else{
    print BLUE, "\nCreating concatenated raw prb file:\n\t$raw_prb_file\n\n", RESET;
    print LOG "\nCreating concatenated raw prb file:\n\t$raw_prb_file\n\n";

    foreach my $tile (sort {$a cmp $b} keys %prb_files){
      my $file = $prb_files{$tile}{file};
      my $file_name = $prb_files{$tile}{name};

      print YELLOW, "\n\tCat: $file", RESET;

      #If file is a _prb.txt.bz2 file, copy it to the working dir, bunzip2 it, cat it to the raw_read_file, and then delete it
      my $bzipped = 0;
      if ($file =~ /_prb\.txt\.bz2$/){
	$bzipped = 1;
	my $cp_cmd = "cp $file $raw_prb_dir";
	my $bunzip2_cmd = "bunzip2 "."$raw_prb_dir"."$file_name";
	system($cp_cmd);
	system($bunzip2_cmd);
	$file = "$raw_prb_dir"."$file_name";
	$file=~s/\.bz2$//; #the bunzip2'd file will have a new name
      }

      my $cmd = "cat $input_dir"."$file >> $raw_prb_file";
      system ($cmd);

      #If the file was bzipped remove the temp file created
      if ($bzipped == 1){
	my $rm_cmd = "rm $file";
	system($rm_cmd);
      }

    }
  }

  print "\n\n";

  return(@return_files);
}


#######################################################################################################################################
#Check file - If the file is already there, delete it
#######################################################################################################################################
sub checkFile{
  my %args = @_;
  my $file = $args{'-file'};

  if (-e $file && -s $file){

    print YELLOW, "\nFile ($file) already exists.  Overwrite (y/n)? ", RESET;
    my $answer = <>;
    chomp($answer);

      if ($answer =~ /^y$|^yes$/i){
	my $cmd = "rm -f $file";
	print YELLOW, "\nDeleting $file\n\n", RESET;
	system($cmd);
      }else{
	print RED, "\nNo point in continuing... exiting\n\n", RESET;
	exit();
      }
    }
  return();
}
