#!/usr/bin/perl -w
#Written by Kim Wong
#Modified heavily by Malachi Griffith
#Copyright 2009 Malachi Griffith and Kim Wong
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
my $working_dir = '';
my $velvet_dir = '';
my $read_type = '';
my $hash_lengths = '';
my $input_file = '';
my $insert_size = '';
my $coverage_cutoff = '';
my $min_contig_length = '';
my $blat_dir = '';
my $reference_dir = '';
my $logfile = '';

GetOptions ('working_dir=s'=>\$working_dir, 'velvet_dir=s'=>\$velvet_dir, 'read_type=s'=>\$read_type, 'hash_lengths=s'=>\$hash_lengths,
	    'input_file=s'=>\$input_file, 'insert_size=i'=>\$insert_size, 'coverage_cutoff=f'=>\$coverage_cutoff, 'min_contig_length=i'=>\$min_contig_length,
	    'blat_dir=s'=>\$blat_dir, 'reference_dir=s'=>\$reference_dir, 'logfile=s'=>\$logfile);


#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify a working directory for output files and sub-directories using: --working_dir", RESET;
print GREEN, "\n\tSpecify the complete path to the velvet binary directory using: --velvet_dir", RESET;
print GREEN, "\n\tSpecify the read type using: --read_type=short or --read_type=shortPaired", RESET;
print GREEN, "\n\t\tNote that if using paired reads, each read1 in the input fasta file must be DIRECTLY followed by its corresponding read2", RESET;
print GREEN, "\n\tSpecify the hash length(s) you wish to use as a space seperated list.  e.g --hash_lengths='19 21 23 25 27 29 31'", RESET;
print GREEN, "\n\t\tNote that these values should be odd numbers and must be 31 or less", RESET;
print GREEN, "\n\tSpecify the input file containing your read sequences using: --input_file", RESET;
print GREEN, "\n\tIf you are using paired reads mode, you must specify a maximum expected insert size using: --insert_size", RESET;
print GREEN, "\n\tSpecify the desired coverage cutoff for contigs using: --coverage_cutoff (say 5)", RESET;
print GREEN, "\n\tIf you wish to limit output to only those contigs over some minimum size N, use: --min_contig_length=N (say 50)", RESET;
print GREEN, "\n\tSpecify the full path to a directory containing a BLAT binary using: --blat_dir", RESET;
print GREEN, "\n\tSpecify the full path to a directory containing reference .nib files (one per chromosome) for BLAT alignment using: --reference_dir", RESET;
print GREEN, "\n\tSpecify the full path to a logfile which will save your run parameters for reference using:  --logfile", RESET;
print GREEN, "\n\nUnpaired mode example:", RESET;
print GREEN, "\n\trunVelvet.pl  --working_dir=/projects/malachig/solexa/velvet_analysis/HS04391/unpaired  --velvet_dir=/home/malachig/tools/Velvet/x64/velvet_0.6.01/  --read_type=short  --hash_lengths='19 21 23 25 27 29 31'  --input_file=/projects/malachig/solexa/fasta_seq_data/HS04391/20821AAXX_Lane1_QualityFiltered.fa  --coverage_cutoff=5  --min_contig_length=50  --blat_dir=/home/malachig/tools/blat/blat34/  --reference_dir=/projects/malachig/sequence_databases/hg18_genome_blatdb/  --logfile=/projects/malachig/solexa/velvet_analysis/HS04391/runVelvet_unpaired_LOG.txt", RESET;
print GREEN, "\n\nPaired mode example:", RESET;
print GREEN, "\n\trunVelvet.pl  --working_dir=/projects/malachig/solexa/velvet_analysis/HS04391/paired  --velvet_dir=/home/malachig/tools/Velvet/x64/velvet_0.6.01/  --read_type=shortPaired  --hash_lengths='19 21 23 25 27 29 31'  --input_file=/projects/malachig/solexa/fasta_seq_data/HS04391/20821AAXX_Lane1_QualityFiltered.fa  --insert_size=500  --coverage_cutoff=5  --min_contig_length=50  --blat_dir=/home/malachig/tools/blat/blat34/  --reference_dir=/projects/malachig/sequence_databases/hg18_genome_blatdb/  --logfile=/projects/malachig/solexa/velvet_analysis/HS04391/runVelvet_unpaired_LOG.txt\n\n", RESET;

unless ($working_dir && $velvet_dir && ($read_type =~ /short|shortPaired/) && $hash_lengths && $input_file && $coverage_cutoff && ($min_contig_length =~ /\d+/) && $blat_dir && $reference_dir && $logfile){
  print RED, "\nRequired parameter missing on not understood!\n\n", RESET;
  exit();
}
chomp($read_type);
if ($read_type =~ /^shortPaired$/){
  unless ($insert_size =~ /\d+/){
    print RED, "\nIf the paired read type is specified you must also specify an insert size with: --insert_size\n\n", RESET;
  }
}

#Add trailing '/' to directory paths if they were forgotten
unless ($working_dir =~ /.*\/$/){
  $working_dir = "$working_dir"."/";
}
unless ($velvet_dir =~ /.*\/$/){
  $velvet_dir = "$velvet_dir"."/";
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
unless (-d $velvet_dir && -e $velvet_dir){
  print RED, "\nVelvet binary directory does not appear to be valid:\n\tvelvet_dir = $velvet_dir\n\n", RESET;
  exit();
}
unless (-d $blat_dir && -e $blat_dir){
  print RED, "\nBLAT binary directory does not appear to be valid:\n\tblat_dir = $blat_dir\n\n", RESET;
  exit();
}
unless (-e $input_file){
  print RED, "\nInput fasta file does not appear to be valid:\n\tinput_file = $input_file\n\n", RESET;
  exit();
}
unless (-d $reference_dir && -e $reference_dir){
  print RED, "\nReference directory does not appear to be valid:\n\treference_dir = $reference_dir\n\n", RESET;
  exit();
}

#Store run parameters in a log file
open(LOG, ">$logfile") || die "\n\nCould not open logfile: $logfile\n\n";
print LOG "The following parameters were specified to runVelvet.pl:\n\n\tworking_dir = $working_dir\n\tvelvet_dir = $velvet_dir\n\tread_type = $read_type\n\thash_lengths = $hash_lengths\n\tinput_file = $input_file\n\tinsert_size = $insert_size\n\tcoverage_cutoff = $coverage_cutoff\n\tmin_contig_length = $min_contig_length\n\tblat_dir = $blat_dir\n\treference_dir = $reference_dir\n\tlogfile = $logfile\n\n";

my $velveth_bin = "$velvet_dir"."velveth";
my $velvetg_bin = "$velvet_dir"."velvetg";

my @hash_lengths = split(" ", $hash_lengths);

#First run velveth at each of the specified k values (hash lengths)
print BLUE, "\nRunning velveth with k values $hash_lengths...\n\n", RESET;
print LOG "\nRunning velveth with k values $hash_lengths...\n\n";

foreach my $k (@hash_lengths){
  chdir($working_dir);
  mkdir("k$k");

  my $file_type;
  if ($input_file =~ /\.fa$/){
    $file_type = "-fasta";
  }elsif($input_file =~ /\.fa\.gz$/){
    $file_type = "-fasta.gz";
  }else{
    print RED, "\nInput file extension not understood: $input_file\n\n", RESET;
    exit();
  }

  if ($read_type =~ /^short$/){
    #Unpaired mode
    my $cmd = "$velveth_bin k$k $k $file_type -short $input_file > $working_dir"."k$k/k$k"."_velveth_logile.txt";
    print BLUE, "\nExecuting: $cmd", RESET;
    print LOG "\nExecuting: $cmd";
    system($cmd);

  }else{
    #Paired mode
    my $cmd = "$velveth_bin k$k $k $file_type -shortPaired $input_file > $working_dir"."k$k/k$k"."_velveth_logile.txt";
    print BLUE, "\nExecuting: $cmd", RESET;
    print LOG "\nExecuting: $cmd";
    system($cmd);
  }
}

#Now run velvetg at each of the specified k values (hash lengths)
print BLUE, "\nRunning velvetg with k values $hash_lengths...\n\n", RESET;
print LOG "\nRunning velvetg with k values $hash_lengths...\n\n";

foreach my $k (@hash_lengths){
  chdir($working_dir);

  if ($read_type =~ /^short$/){
    #Unpaired mode
    my $cmd = "$velvetg_bin k$k -cov_cutoff $coverage_cutoff -min_contig_lgth $min_contig_length > $working_dir"."k$k/k$k"."_velvetg_logile.txt";
    print BLUE, "\nExecuting: $cmd", RESET;
    print LOG "\nExecuting: $cmd";
    system($cmd);

    # count number of contigs
    my $contigs = `grep ">" k$k/contigs.fa | wc -l`;
    print BLUE, "\n\tFound the following number of contigs: $contigs", RESET;
    print LOG "\n\tFound the following number of contigs: $contigs";

  }else{
    #Paired mode
    my $cmd = "$velvetg_bin k$k -cov_cutoff $coverage_cutoff -ins_length $insert_size -min_contig_lgth $min_contig_length > $working_dir"."k$k/k$k"."_velvetg_logile.txt";
    print BLUE, "\nExecuting: $cmd", RESET;
    print LOG "\nExecuting: $cmd";
    system($cmd);

    # count number of contigs
    my $contigs = `grep ">" k$k/contigs.fa | wc -l`;
    print BLUE, "\n\tFound the following number of contigs: $contigs", RESET;
    print LOG "\n\tFound the following number of contigs: $contigs";
  }
}

#Now BLAT the resulting contigs against the specified reference sequence (i.e. complete reference human genome)
print BLUE, "\nUsing BLAT to align contigs to the specified reference sequences in: $reference_dir...\n\n", RESET;
print BLUE, "\nFirst making a list of all contigs for each k value used", RESET;
print LOG "\nUsing BLAT to align contigs to the specified reference sequences in: $reference_dir...\n\n";
print LOG "\nFirst making a list of all contigs for each k value used";

#Create a list of all contigs generated
foreach my $k (@hash_lengths){
  chdir($working_dir);

  my $makelist_cmd = "grep \">\" k$k/contigs.fa | cut -f 2 -d \">\" > $working_dir"."k$k/contigslist.txt";

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

    my $contig_fasta_file = "$working_dir"."k$k/contigs.fa";
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

  my $gzip_cmd = "gzip $working_dir"."k$k"."/*";
  system($gzip_cmd);
}



#EXONERATE PROCESSING - Do not attempt this for now
#print STDERR "Running Exonerate gapped alignment. This might take a while...\n";

# run exonerate
#my $exoncommand = qq(/home/pubseq/BioSw/exonerate/exonerate-2.0.0_x86_64/bin/exonerate  --model affine:local --bestn 2 --showvulgar 0 --ryo \"%S %r %ql %tl %em %V\\n\" k\$k/contigs.fa $ref > k\$k/contigs.exon);
#print STDERR   qq(/home/pubseq/BioSw/exonerate/exonerate-2.0.0_x86_64/bin/exonerate  --model affine:local --bestn 1 --showvulgar 0 --ryo \"%S %r %ql %tl %em %V\\n\" k\$k/contigs.fa $ref );
#my $runthis = qq(
#for k in $ks; 
#do $exoncommand;
#echo '$exoncommand'  >> k\$k/logfile;
#done;
#);
#print STDERR "$runthis\n";
#`$runthis`;


# run exonerate parser

#my $exonparsecommand = "/home/kwong/svn/assemblies/parseExonerate2.pl k\$k/contigslist k\$k/contigs.exon > k\$k/contigs.exon.stats";
#system qq(
#for k in $ks; 
#do $exonparsecommand;
#echo '$exonparsecommand'  >> k\$k/logfile;
#done;
#);
#print STDERR "\nDONE!\n\n";

close(LOG);

exit();





