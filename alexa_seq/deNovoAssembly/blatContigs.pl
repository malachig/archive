#!/usr/bin/perl -w
#Written by Kim Wong
#Modified heavily by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#This script is a wrapper for running the Velvet de novo assembler for short read sequences
#Example usage of Velvet:

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);

#Initialize command line options
my $library = '';
my $working_dir = '';
my $hash_lengths = '';
my $blat_dir = '';
my $reference_dir = '';
my $logfile = '';

GetOptions ('library=s'=>\$library, 'working_dir=s'=>\$working_dir, 'hash_lengths=s'=>\$hash_lengths,
	    'blat_dir=s'=>\$blat_dir, 'reference_dir=s'=>\$reference_dir, 'logfile=s'=>\$logfile);


#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify a working directory for output files and sub-directories using: --working_dir", RESET;
print GREEN, "\n\tSpecify the hash length(s) that were used.  e.g --hash_lengths='19 21 23 25 27 29 31'", RESET;
print GREEN, "\n\tSpecify the full path to a directory containing a BLAT binary using: --blat_dir", RESET;
print GREEN, "\n\tSpecify the full path to a directory containing reference .nib files (one per chromosome) for BLAT alignment using: --reference_dir", RESET;
print GREEN, "\n\tSpecify the full path to a logfile which will save your run parameters for reference using:  --logfile", RESET;
print GREEN, "\n\nUnpaired mode example:", RESET;
print GREEN, "\n\tblatContigs.pl  --library=HS0439  --working_dir=/projects/malachig/solexa/abyss_analysis/HS04391_16Lanes/  --hash_lengths='22 24 26 28 30 32 34 36 38 40 42'  --blat_dir=/home/malachig/tools/blat/blat34/  --reference_dir=/projects/malachig/sequence_databases/hg18_genome_blatdb/  --logfile=/projects/malachig/solexa/abyss_analysis/HS04391_16Lanes/blatContigs_LOG.txt", RESET;

unless ($library && $working_dir && $hash_lengths && $blat_dir && $reference_dir && $logfile){
  print RED, "\n\nRequired parameter missing on not understood!\n\n", RESET;
  exit();
}

#Add trailing '/' to directory paths if they were forgotten
unless ($working_dir =~ /.*\/$/){
  $working_dir = "$working_dir"."/";
}
unless ($blat_dir =~ /.*\/$/){
  $blat_dir = "$blat_dir"."/";
}
unless ($reference_dir =~ /.*\/$/){
  $reference_dir = "$reference_dir"."/";
}

#Check the validity of required directories and binaries specified by the user
unless (-d $working_dir && -e $working_dir){
  mkdir($working_dir);
}
unless (-d $blat_dir && -e $blat_dir){
  print RED, "\nBLAT binary directory does not appear to be valid:\n\tblat_dir = $blat_dir\n\n", RESET;
  exit();
}
unless (-d $reference_dir && -e $reference_dir){
  print RED, "\nReference directory does not appear to be valid:\n\treference_dir = $reference_dir\n\n", RESET;
  exit();
}

#Store run parameters in a log file
open(LOG, ">$logfile") || die "\n\nCould not open logfile: $logfile\n\n";
print LOG "The following parameters were specified to blatContigs.pl:\n\n\tworking_dir = $working_dir\n\thash_lengths = $hash_lengths\n\tblat_dir = $blat_dir\n\treference_dir = $reference_dir\n\tlogfile = $logfile\n\n";

my @hash_lengths = split(" ", $hash_lengths);

#Now BLAT the resulting contigs against the specified reference sequence (i.e. complete reference human genome)
print BLUE, "\n\nUsing BLAT to align contigs to the specified reference sequences in: $reference_dir...\n\n", RESET;
print BLUE, "\nFirst making a list of all contigs for each k value used", RESET;
print LOG "\nUsing BLAT to align contigs to the specified reference sequences in: $reference_dir...\n\n";
print LOG "\nFirst making a list of all contigs for each k value used";

#Create a list of all contigs generated
foreach my $k (@hash_lengths){
  chdir($working_dir);

  my $makelist_cmd = "cat k$k/$library-contigs.fa | grep \">\" | cut -f 2 -d \">\" > $working_dir"."k$k/contigslist.txt";

  print BLUE, "\nExecuting: $makelist_cmd\n", RESET;
  print LOG "\nExecuting: $makelist_cmd\n";
  system ("$makelist_cmd");
}

#Run blat - Standalone blat has memory issues and will need to be run chromosome by chromosome
# Note that the '-fastMap' may result in small contigs failing to align
# - '-fastMap' will also do ungapped alignments i believe, so is not ideal for transcript sequences.  May still be suitable for crudely assessing assembly success??


print BLUE, "\nNow running the actual BLAT job for each k value used", RESET;
print LOG "\nNow running the actual BLAT job for each k value used";

foreach my $k (@hash_lengths){

  chdir($working_dir);

  #Make a temp dir for BLAT output files (one psl file per chromosome)
  my $temp_dir = "$working_dir"."k$k"."/blat_temp/";
  mkdir($temp_dir);

  #Get all .nib files from the specified reference directory

  opendir(DIRHANDLE, "$reference_dir") || die "\nCannot open reference directory: $reference_dir\n\n";
  my @files = readdir(DIRHANDLE);
  my @ref_files;

  foreach my $file (@files){

    #Skip all files except .nib files
    unless ($file =~ /\.nib$/){
      next();
    }
    my $file_path = "$reference_dir"."$file";
    push (@ref_files, $file_path);
  }

  my $ref_file_count = scalar(@ref_files);
  print BLUE, "\nFound $ref_file_count reference .nib files\n\n", RESET;
  print LOG "\nFound $ref_file_count reference .nib files\n\n";
  closedir(DIRHANDLE);


  #Run blat once for each .nib file (chromosome) in the source reference directory
  my $blat_log_file = "$working_dir"."k$k/k$k"."_blat_logile.txt";
  foreach my $ref_file (sort @ref_files){

    my $chr;
    if ($ref_file =~ /chr(.*)\.nib/){
      $chr = $1;
    }else{
      print RED, "\nReference file name not understood: $ref_file\n\n", RESET;
      exit();
    }

    my $contig_fasta_file = "$working_dir"."k$k/$library-contigs.fa";
    my $psl_out_file = "$temp_dir"."contigs_fastMap.chr$chr.psl";

    #With '-fastMap' option
    my $blat_cmd =  "$blat_dir"."blat -q=dna -fastMap -noHead $ref_file $contig_fasta_file $psl_out_file >> $blat_log_file";

    #Without '-fastMap' option
    #my $blat_cmd =  "$blat_dir"."blat -q=dna -noHead $ref_file $contig_fasta_file $psl_out_file >> $blat_log_file";


    print YELLOW, "\n\tExecuting: $blat_cmd\n", RESET;
    print LOG "\n\tExecuting: $blat_cmd\n";

    system($blat_cmd);
  }

  #Grab and concatenate all of the individual BLAT files to one grand file and store in the main directory
  my $grand_psl_out_file = "$working_dir"."k$k/"."k$k"."_psl_outfile.txt";
  print BLUE, "\nCreating concatenated PSL output file:\n\t$grand_psl_out_file\n\n", RESET;
  print LOG "\nCreating concatenated PSL output file:\n\t$grand_psl_out_file\n\n";

  opendir(DIRHANDLE, "$temp_dir") || die "\nCannot open reference directory: $temp_dir\n\n";
  @files = readdir(DIRHANDLE);
  my @psl_files;

  foreach my $psl_file (sort @files){

    #Skip all files except .nib files
    unless ($psl_file =~ /\.psl$/){
      next();
    }
    my $file_path = "$temp_dir"."$psl_file";
    push (@psl_files, $file_path);
  }

  my $psl_file_count = scalar(@psl_files);
  print BLUE, "\nFound $psl_file_count reference .psl files\n\n", RESET;
  print LOG "\nFound $psl_file_count reference .psl files\n\n";
  closedir(DIRHANDLE);


  foreach my $psl_file (sort @psl_files){
    my $cmd = "cat $psl_file >> $grand_psl_out_file";
    system ($cmd);
  }

  #Make sure one .psl results file was generated from every .nib file
  unless ($ref_file_count == $psl_file_count){
    print RED, "\nNumber of resulting .psl files does not match number of .nib reference files!!\n\n", RESET;
    exit();
  }

  #Delete the temp directory containing all the individual BLAT files
  my $rm_cmd = "rm -fr $temp_dir";
  system($rm_cmd);
}

print BLUE, "\nNow parsing the PSL output of the BLAT job and summarizing using Kim's PSL parser\n\n", RESET;
print LOG "\nNow parsing the PSL output of the BLAT job and summarizing using Kim's PSL parser\n\n";
foreach my $k (@hash_lengths){
  my $grand_psl_out_file = "$working_dir"."k$k/"."k$k"."_psl_outfile.txt";

  chdir($working_dir);
  my $stats_file = "$working_dir"."k$k"."_contigs.BlatFastMap.stats";
  my $pslparse_cmd = "/home/malachig/svn/solexa_analysis/parsepsl.pl k$k/contigslist.txt $grand_psl_out_file > $stats_file";
  system ($pslparse_cmd);
  print BLUE, "\n\tStats were stored in: $stats_file", RESET;
  print LOG "\n\tStats were stored in: $stats_file";
}

print "\n\n";
print LOG "\n\n";

#Compress the files in the working directory to save disk space
print BLUE, "\nNow compressing files in the working directories created\n\n", RESET;
print LOG "\nNow compressing files in the working directories created\n\n";
foreach my $k (@hash_lengths){
  chdir($working_dir);

  my $gzip_cmd = "gzip -f $working_dir"."k$k"."/*";
  system($gzip_cmd);
}


close(LOG);

exit();





