#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#Take an input fasta file, divide into pieces and create a map batch job

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;

#Load the ALEXA modules
BEGIN {
  use Cwd 'abs_path';
  if (abs_path($0) =~ /(.*)\/.*\/.*\.pl/){
    push (@INC, $1);
  }
}
use utilities::utility qw(:all);

my $input_fasta_file = '';
my $working_dir = '';
my $map_bin = '';           #Full path to the map executable to be used in the map command - blastall or bwa
my $map_database = '';      #Full path to map database created with the same version of BLAST or BWA specied as map_bin
my $word_size = '';         #Word search size (default for blastall is 11, not applicable for bwa)
my $min_bit_score = '';     #Alignment output will be piped to a filter script which will drop all hits that dont meet this minimum (use 0 for no filter)
my $map_filter_script = ''; #Location of map filter script
my $batch_file = '';        #Batch script file to create
my $map_results_dir = '';
my $job_size = '';
my $target_size = '';
my $min_overlap = '';
my $aligner = '';

GetOptions ('input_fasta_file=s'=>\$input_fasta_file, 'working_dir=s'=>\$working_dir,
	    'map_bin=s'=>\$map_bin, 'map_database=s'=>\$map_database, 'word_size=i'=>\$word_size,
	    'min_bit_score=f'=>\$min_bit_score, 'map_filter_script=s'=>\$map_filter_script,
	    'batch_file=s'=>\$batch_file, 'map_results_dir=s'=>\$map_results_dir, 'job_size=i'=>\$job_size,
            'target_size=i'=>\$target_size, 'min_overlap=i'=>\$min_overlap, 'aligner=s'=>\$aligner);

print GREEN, "\n\nExample usage:\n\ncreateMapBatch.pl  --input_fasta_file=/projects/malachig/solexa/fasta_seq_data/HS04391/HS04391_Lanes1-8_QualityFiltered_Unpaired.fa  --working_dir=/projects/malachig/solexa/fasta_seq_data/HS04391/fasta_blocks/  --map_bin=/home/pubseq/BioSw/BLAST2/blast2.2.18_x64/bin/blastall  --map_database=/projects/malachig/sequence_databases/ensembl_transcripts_hs_49_36k/alexa_transcripts  --word_size=11  --min_bit_score=40.0  --map_filter_script=/home/malachig/svn/solexa_analysis/filterBlastStream.pl  --batch_file=/projects/malachig/solexa/batch_jobs/HS04391/map_versus_EnsemblTranscripts_v49.sh  --map_results_dir=/projects/malachig/solexa/blast_results/HS04391/ensembl_transcripts_v49/  --job_size=250000  [--target_size=62  --min_overlap=6]  --aligner=blast\n\n", RESET;

unless ($input_fasta_file && $working_dir && $map_bin && $map_database && $word_size && ($min_bit_score =~ /^\d+/) && $map_filter_script && $batch_file && $map_results_dir && $job_size && $aligner){
  print RED, "\nRequired input parameter(s) missing!\n\n", RESET;
  exit();
}
unless ($aligner =~ /blast|bwa/){
  print RED, "\n--aligner option must be either 'bwa' or 'blast'\n\n", RESET;
  exit();
}
chomp($min_bit_score);
unless ($min_bit_score =~ m/^\d+.\d+$/ || $min_bit_score =~ m/^\d+$/){
  print "\n$map_filter_script was supplied a --min_bit_score value that does not appear to be a valid number!\n\n";
  exit();
}

#Check validity of input files and directories
unless (-e $map_filter_script){
  print RED, "\nInput map filter script does not appear to be valid!\n\n", RESET;
  exit();
}
unless (-e $input_fasta_file){
  print RED, "\nInput fasta file does not appear to be valid!\n\n", RESET;
  exit();
}
unless ($input_fasta_file =~ /(.*)\.gz$/){
  print RED, "\nFound an uncompressed file: $input_fasta_file\n\n\tMake sure all files are compressed before proceeding\n\n\t- A mix of compressed and uncompressed files may indicate a problem (i.e. you need to figure out which is complete and which might be partial!!)\n\n", RESET;
  exit();
}
unless (-e $map_bin){
  print RED, "\nMap binary ($map_bin) does not appear to be valid!\n\n", RESET;
  exit();
}


#See if fasta files are already present in the working directory
my $fasta_present = 0;
my $line_count = `zcat $input_fasta_file | wc -l`;
chomp($line_count);
my $expected_files = $line_count/$job_size;

my $floor_expected_files;
if ($expected_files =~ /^(\d+)\.\d+$/){
  $floor_expected_files = $1;
}elsif($expected_files =~ /^(\d+)$/){
  $floor_expected_files = $1;
}else{
  print RED, "\nNumber not understood: $expected_files\n\n", RESET;
  exit();
}

if ($floor_expected_files == $expected_files){
  $expected_files = $floor_expected_files;
}else{
  $expected_files = $floor_expected_files+1;
}

opendir(DIRHANDLE, "$working_dir") || die "\nCannot open directory: $working_dir\n\n";
my @files = readdir(DIRHANDLE);
closedir(DIRHANDLE);
foreach my $file(@files){
  if ($file =~ /fasta/){
    $fasta_present++;
  }
}
print BLUE, "\nFound $fasta_present fasta files already in directory - expect ($expected_files)\n\n", RESET;
#Split the input fasta file into pieces using unix split
if ($fasta_present == $expected_files){
  print YELLOW, "\nLeaving existing fasta and map files!!\n\n", RESET;
  $map_results_dir = &checkDir('-dir'=>$map_results_dir, '-clear'=>"no");
}else{
  print YELLOW, "\nGenerating new fasta files!!\n\n", RESET;
  $working_dir = &checkDir('-dir'=>$working_dir, '-clear'=>"yes");
  my $split_cmd = "zcat $input_fasta_file | /usr/bin/split -a 4 -d -l $job_size - fasta_";
  print BLUE, "\nExecuting:  $split_cmd\n", RESET;
  chdir($working_dir);
  system ($split_cmd);
  $map_results_dir = &checkDir('-dir'=>$map_results_dir, '-clear'=>"yes");
}

#Get a list of the files created
opendir(DIRHANDLE, "$working_dir") || die "\nCannot open directory: $working_dir\n\n";
@files = readdir(DIRHANDLE);
closedir(DIRHANDLE);

#Make sure the file was divided into too many pieces for the naming scheme
my $file_count = scalar(@files);

print BLUE, "\nFound $file_count files in the working directory\n", RESET;
if ($file_count >= 9999){
  print RED, "\nNumber of files suggests that split ran out of numbers when dividing the input file - use a larger job_size!!\n\n", RESET;
  exit();
}

#Get a list of fasta block files just created and organize them
my %files;
foreach my $file_name (@files){
  unless ($file_name =~ /^fasta_/){
    print YELLOW, "\n$file_name does not match expected pattern, skipping", RESET;
    next();
  }

  my $file_num;
  if ($file_name =~ /fasta_(\d+)/){
    $file_num = $1;
  }else{
    print RED, "\nFile name: $file_name not understood", RESET;
    exit();
  }

  $files{$file_num}{name} = $file_name;
  $files{$file_num}{path} = "$working_dir"."$file_name";

}

#Open the batch output file and print a line for each fasta file to be processed
print BLUE, "\n\nPrint map batch file: $batch_file\n", RESET;
open (BATCH, ">$batch_file") || die "\nCould not open output batch file: $batch_file\n\n";
foreach my $file_num (sort keys %files){

  my $file_path = $files{$file_num}{path};

  my $map_results_file = "$map_results_dir"."blast_"."$file_num";

  if ($aligner =~ /blast/){
    if ($target_size && $min_overlap){
      print BATCH "source ~/.bashrc; $map_bin -p blastn -d $map_database -i $file_path -m 8 -F F -W $word_size | $map_filter_script --min_bit_score=$min_bit_score --centre_span=$min_overlap --target_size=$target_size > $map_results_file; gzip -f $map_results_file\n";
    }else{
      print BATCH "source ~/.bashrc; $map_bin -p blastn -d $map_database -i $file_path -m 8 -F F -W $word_size | $map_filter_script --min_bit_score=$min_bit_score > $map_results_file; gzip -f $map_results_file\n";
    }
  }elsif($aligner =~ /bwa/){
    #IMPORTANT: The -n option supplied to samse should be greater than the max number of known transcripts for a single gene in the transcriptome being analyzed
    print BATCH "source ~/.bashrc; $map_bin aln $map_database $file_path 2>/dev/null | $map_bin samse -n 100 $map_database - $file_path  2>/dev/null | $map_filter_script --min_bit_score=$min_bit_score > $map_results_file; gzip -f $map_results_file\n";

  }

}
close (BATCH);

print BLUE, "\n\nTo submit job, log into cluster head node and submit the jobs:\n", RESET;
print BLUE, "\nTo run without a cluster use: bash $batch_file\n\n", RESET;

exit();

