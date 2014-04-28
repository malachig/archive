#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);

#Load the ALEXA modules
BEGIN {
  use Cwd 'abs_path';
  if (abs_path($0) =~ /(.*)\/.*\.pl/){
    push (@INC, $1);
  }
}
use utilities::utility qw(:all);

my $fastq_dir = '';
my $working_dir = '';
my $maq_bin = '';
my $batch_file_dir = '';

GetOptions ('fastq_dir=s'=>\$fastq_dir, 'working_dir=s'=>\$working_dir, 'maq_bin=s'=>\$maq_bin, 'batch_file_dir=s'=>\$batch_file_dir);

#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\nExample: createMaqMapBatch.pl  --fastq_dir=/projects/malachig/solexa/maq_analysis/HS04391/  --working_dir=/projects/malachig/solexa/maq_analysis/HS04391/map_temp/  --maq_bin=/home/malachig/tools/MAQ/maq-0.6.6_i686-linux/maq  --batch_file_dir=/projects/malachig/solexa/batch_jobs/HS04391/MAQ_batchwise/  \n\n", RESET;

unless ($fastq_dir && $working_dir && $maq_bin && $batch_file_dir){
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  exit();
}

unless (-e $maq_bin){
  print RED, "\nSpecified MAQ binary doesnt seem valid\n\n", RESET;
  exit();
}

#Check working dir
$working_dir = &checkDir('-dir'=>$working_dir, '-clear'=>"no");
$fastq_dir = &checkDir('-dir'=>$fastq_dir, '-clear'=>"no");
$batch_file_dir = &checkDir('-dir'=>$batch_file_dir, '-clear'=>"no");

#Hard code paths to reference databases
my $reference_db_36mer = "/projects/rmorin/common/genomes/all_human_plus_junctions_36nt.bfa";
my $reference_db_42mer = "/projects/rmorin/common/genomes/all_human_plus_junctions_42nt.bfa";


#Open output batch files
#binary_fastq_1.sh  maq_map_2.sh  maq_merge_3.sh
my $binary_fq_file = "$batch_file_dir"."binary_fastq_1.sh";
my $maq_map_file = "$batch_file_dir"."maq_map_2.sh";
my $map_merge_file = "$batch_file_dir"."maq_merge_3.sh";

open (BFQ_OUT, ">$binary_fq_file") || die "\nCould not open $binary_fq_file\n\n";
open (MAP_OUT, ">$maq_map_file") || die "\nCould not open $maq_map_file\n\n";
open (MERGE_OUT, ">$map_merge_file") || die "\nCould not open $map_merge_file\n\n";

my %run_names;

#Look for matching R1 and R2 fastq files (assumed to be gzipped)
opendir(DIRHANDLE, "$fastq_dir") || die "\nCannot open directory: $fastq_dir\n\n";
my @temp = readdir(DIRHANDLE);
closedir(DIRHANDLE);

foreach my $file (@temp){
  if ($file =~ /(\w+_\w+)_R\d\.fastq\.gz/){
    $run_names{$1}{count}++;
  }
}

foreach my $run_name (sort {$a cmp $b} keys %run_names){

  print YELLOW, "\n\nPROCESSING: $run_name", RESET;

  opendir(DIRHANDLE, "$fastq_dir") || die "\nCannot open directory: $fastq_dir\n\n";
  my @temp = readdir(DIRHANDLE);
  closedir(DIRHANDLE);

  my $r1_name = "$run_name"."_R1.fastq.gz";
  my $r2_name = "$run_name"."_R2.fastq.gz";
  my $r1_file_path;
  my $r2_file_path;
  my $r1_file_path_new;
  my $r2_file_path_new;
  my $r1_file_found = 0;
  my $r2_file_found = 0;

  foreach my $file (@temp){

    if ($file =~ /$r1_name/){
      $r1_file_path = "$fastq_dir"."$r1_name";
      $r1_file_path_new = "$working_dir"."$r1_name";
      print YELLOW, "\nFound: $r1_file_path", RESET;
      $r1_file_found = 1;
    }
    if ($file =~ /$r2_name/){
      $r2_file_path = "$fastq_dir"."$r2_name";
      $r2_file_path_new = "$working_dir"."$r2_name";
      print YELLOW, "\nFound: $r2_file_path", RESET;
      $r2_file_found = 1;
    }
  }

  unless ($r1_file_found == 1 && $r2_file_found == 1){
    print RED, "\nCould not find expected files: $r1_name and $r2_name\n\n", RESET;
    exit();
  }

  #Clear the working directory of files for this run name
  my $cmd_clear = "rm -f $working_dir"."$run_name"."*";
  print BLUE, "\n\n$cmd_clear\n\n", RESET;
  system($cmd_clear);

  #Copy both files to the working dir
  my $cmd_r1_cp = "cp -f $r1_file_path $r1_file_path_new";
  print BLUE, "\n\n$cmd_r1_cp", RESET;
  system($cmd_r1_cp);

  my $cmd_r2_cp = "cp -f $r2_file_path $r2_file_path_new";
  print BLUE, "\n$cmd_r2_cp", RESET;
  system($cmd_r2_cp);

  #split these files into pieces
  chdir($working_dir);
  my $cmd_r1_split = "gunzip -c $r1_file_path_new | split -l 800000 -a 3 - $run_name"."_R1_";
  print BLUE, "\n\n$cmd_r1_split", RESET;
  system($cmd_r1_split);

  my $cmd_r2_split = "gunzip -c $r2_file_path_new | split -l 800000 -a 3 - $run_name"."_R2_";
  print BLUE, "$cmd_r2_split", RESET;
  system($cmd_r2_split);

  #Delete the copies of fastq files no longer needed
  my $cmd_r1_fastq_rm = "rm -f $r1_file_path_new";
  print BLUE, "\n\n$cmd_r1_fastq_rm", RESET;
  system($cmd_r1_fastq_rm);

  my $cmd_r2_fastq_rm = "rm -f $r2_file_path_new";
  print BLUE, "\n$cmd_r2_fastq_rm\n", RESET;
  system($cmd_r2_fastq_rm);

  #Get all the fastq file pieces
  opendir(DIRHANDLE, "$working_dir") || die "\nCannot open directory: $working_dir\n\n";
  @temp = readdir(DIRHANDLE);
  closedir(DIRHANDLE);

  my %r1_files;
  my %r2_files;

  $r1_name = "$run_name"."_R1_";
  $r2_name = "$run_name"."_R2_";

  my $r1_count = 0;
  my $r2_count = 0;

  foreach my $file (@temp){
    if ($file =~ /($r1_name)(\w+)/){
      $r1_count++;
      $r1_files{$2}{name} = "$1$2";
      $r1_files{$2}{path} = "$working_dir"."$1$2";
      $r1_files{$2}{bfq_path} = "$working_dir"."$1$2".".bfq";
    }
    if ($file =~ /($r2_name)(\w+)/){
      $r2_count++;
      $r2_files{$2}{name} = "$1$2";
      $r2_files{$2}{path} = "$working_dir"."$1$2";
      $r2_files{$2}{bfq_path} = "$working_dir"."$1$2".".bfq";
    }
  }

  #Determine the read length for each file
  foreach my $r1_file_code (sort {$a cmp $b} keys %r1_files){

    #Open a read1 file
    open (TEMP, "$r1_files{$r1_file_code}{path}") || die "\nCould not open file: $r1_files{$r1_file_code}{path}\n\n";
    my $line_count = 0;
    while(<TEMP>){
      $line_count++;
      if ($line_count == 2){
	chomp($_);
	my $seq = $_;
	my $r1_read_length = length($seq);
	$r1_files{$r1_file_code}{read_length} = $r1_read_length;
	last();
      }
    }
    close(TEMP);

    #Open the corresponding read2 file
    my $r2_file_code = $r1_file_code;
    open (TEMP, "$r2_files{$r2_file_code}{path}") || die "\nCould not open file: $r2_files{$r2_file_code}{path}\n\n";
    $line_count = 0;
    while(<TEMP>){
      $line_count++;
      if ($line_count == 2){
	chomp($_);
	my $seq = $_;
	my $r2_read_length = length($seq);
	$r2_files{$r2_file_code}{read_length} = $r2_read_length;
	last();
      }
    }
    close(TEMP);
  }

  #Now create batch commands for each of the files pieces
  #Create two batch jobs:
  #A.) Create binary fastq files - then delete the raw fastq files
  #B.) Perform MAQ map operation using the appropriate database based on read length - then delete the binary fasta files 
  #C.) Perform a MAQ merge operation for all map files of a lane - then delete individual map files when done

  my %mapfiles;

  foreach my $r1_file_code (sort {$a cmp $b} keys %r1_files){

    my $r2_file_code = $r1_file_code;

    #A.) Create binary fastq files - then delete the raw fastq files
    print BFQ_OUT "$maq_bin fastq2bfq $r1_files{$r1_file_code}{path} $r1_files{$r1_file_code}{bfq_path}; rm -f $r1_files{$r1_file_code}{path}\n";
    print BFQ_OUT "$maq_bin fastq2bfq $r2_files{$r2_file_code}{path} $r2_files{$r2_file_code}{bfq_path}; rm -f $r2_files{$r2_file_code}{path}\n";

    #B.) Perform MAQ map operation using the appropriate database based on read length - then delete the binary fasta files 
    my $map_file = "$working_dir"."$run_name"."_"."$r1_file_code".".map";

    $mapfiles{$r1_file_code}{map_file_path} = $map_file;

    my $target_ref_db;
    if ($r1_files{$r1_file_code}{read_length} == 36){
      $target_ref_db = $reference_db_36mer;
    }elsif($r1_files{$r1_file_code}{read_length} == 42){
      $target_ref_db = $reference_db_42mer;
    }else{
      print RED, "\nCould not find a matching reference database for length = $r1_files{$r1_file_code}{read_length}\n\n", RESET;
      exit();
    }

    print MAP_OUT "$maq_bin map $map_file $target_ref_db $r1_files{$r1_file_code}{bfq_path} $r2_files{$r2_file_code}{bfq_path}; rm -f $r1_files{$r1_file_code}{bfq_path}; rm -f $r2_files{$r2_file_code}{bfq_path}\n";
  }

  #C.) Perform a MAQ merge operation for all map files of a lane - then delete individual map files when done
  print MERGE_OUT "$maq_bin mapmerge $working_dir"."$run_name".".map ";
  foreach my $file_code (sort {$a cmp $b} keys %mapfiles){
    print MERGE_OUT "$mapfiles{$file_code}{map_file_path} ";
  }
  print MERGE_OUT "\n";
  foreach my $file_code (sort {$a cmp $b} keys %mapfiles){
    print MERGE_OUT "rm -f $mapfiles{$file_code}{map_file_path}\n";
  }

  #print Dumper %r1_files;
  #print Dumper %r2_files;


}


close (BFQ_OUT);
close (MAP_OUT);
close (MERGE_OUT);

exit();





