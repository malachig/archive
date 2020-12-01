#!/usr/bin/perl -w
# Malachi Griffith

#Wrapper that sets up splign jobs for the cluster
#Check all current splign results and only select those junctions that have not already been analyzed
#Use as a starting point, a master merge of all junctions across all projects... or a for a single project


use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;
use IO::File;
use Tie::File;
use Fcntl 'O_RDONLY';
use Benchmark;

my $junction_list = '';
my $splign_bin_dir = '';
my $blast_bin_dir = '';
my $genome_fasta = '';
my $splign_results_dir = '';
my $job_size = '';
my $job_name = '';

GetOptions ('junction_list=s'=>\$junction_list, 'splign_bin_dir=s'=>\$splign_bin_dir, 'blast_bin_dir=s'=>\$blast_bin_dir, 'genome_fasta=s'=>\$genome_fasta, 
            'splign_results_dir=s'=>\$splign_results_dir, 'job_size=i'=>\$job_size, 'job_name=s'=>\$job_name);

if ($junction_list && $splign_bin_dir && $blast_bin_dir && $genome_fasta && $splign_results_dir && $job_size && $job_name){
  #Parameters ok
}else{
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  print GREEN, "\nSpecify the junction list using: --junction_list", RESET;
  print GREEN, "\nSpecify the path to the splign binary dir using: --splign_bin_dir", RESET;
  print GREEN, "\nSpecify the path to the splign binary dir using: --blast_bin_dir", RESET;
  print GREEN, "\nSpecify the path to a junction fasta file using: --junction_fasta", RESET;
  print GREEN, "\nSpecify the path to a genome fasta file using: --genome_fasta", RESET;
  print GREEN, "\nSpecify the path to a splign results dir using: --splign_results_dir", RESET;
  print GREEN, "\nSpecify the number of reads to process as a block using: --job_size", RESET;
  print GREEN, "\nSpecify a base name for jobs using: --job_name", RESET;
  print GREEN, "\n\nExample: createSplignJobs.pl  --junction_list=/projects/alexa2/hmmSplicer/MASTER.junction.list.txt  --splign_bin_dir=/home/malachig/tools/splign/  --blast_bin_dir=/home/malachig/tools/BLAST2/blast2.2.18_x64/bin/  --genome_fasta=/projects/malachig/sequence_databases/hg18_genome/all_human.fa  --splign_results_dir=/projects/alexa2/hmmSplicer/SplignValidation/hg18/  --job_name=Nov25  --job_size=10000\n\n", RESET;
  exit();
}

#Check dirs and files
unless (-e $junction_list){
  print RED, "\n\nJunction list: $junction_list not found\n\n", RESET;
  exit();
}
unless ($splign_bin_dir =~ /\/$/){
  $splign_bin_dir .= "/";
}
unless (-e $splign_bin_dir && -d $splign_bin_dir){
  print RED, "\n\nSplign binary dir: $splign_bin_dir not found\n\n", RESET;
  exit();
}
unless ($blast_bin_dir =~ /\/$/){
  $blast_bin_dir .= "/";
}
unless (-e $blast_bin_dir && -d $blast_bin_dir){
  print RED, "\n\nBLAST binary dir: $blast_bin_dir not found\n\n", RESET;
  exit();
}
unless (-e $genome_fasta){
  print RED, "\n\nGenome fasta: $genome_fasta not found\n\n", RESET;
  exit();
}
unless ($splign_results_dir =~ /\/$/){
  $splign_results_dir .= "/";
}
unless (-e $splign_results_dir && -d $splign_results_dir){
  print RED, "\n\nSplign results dir: $splign_results_dir not found\n\n", RESET;
  exit();
}

#Get a list of all junction IDs already processed
my $splign_master_file = "$splign_results_dir"."splignMasterResults.txt";
my %junc_proc;
open (SPLIGN, "$splign_master_file") || die "\n\nCould not open splign file: $splign_master_file\n\n";
my $header = 1;
while(<SPLIGN>){
  chomp($_);
  my @line = split("\t", $_);
  if ($header == 1){
    $header = 0;
    next();
  }
  $junc_proc{$line[0]}=1;
}
close(SPLIGN);

#Get a list of all junction IDs in the input list - skip those that have already been processed - create fasta files of size $job_size for the remainder
my %junc_new;
my %fasta_files;
open (JUNC, $junction_list) || die "\n\nJunction list could not be opened: $junction_list\n\n";
my $c = 0;
my $block = 1;
$header = 1;
my %columns;
while(<JUNC>){
  chomp($_);
  my @line = split("\t", $_);
  if ($header){
    my $pos = 0;
    foreach my $head (@line){
      $columns{$head}{pos} = $pos;
      $pos++;
    }
    $header = 0;
    next();
  }
  my $jid = $line[$columns{'JID'}{pos}];
  my $seq = $line[$columns{'Sequence'}{pos}];
  if ($junc_proc{$jid}){
    next();
  }
  $c++;
  if($c == 1){
    my $job_id = "$job_name"."_"."$block";
    my $outfile = "$splign_results_dir"."fasta/$job_id".".fa";
    open (OUT, ">$outfile") || die "\n\nCould not open outfile: $outfile\n\n";
    $fasta_files{$block}{outfile} = $outfile;
    $fasta_files{$block}{job_id} = $job_id;
  }

  print OUT ">$jid\n$seq\n";

  if ($c == $job_size){
    $c = 0;
    $block++;
    close(OUT);
  }
}
close(JUNC);


#For each fasta file generated, create a splignRun.pl command for each fasta file and write this command to the output bash file
my $bash_file = "$splign_results_dir"."splignRun_"."$job_name".".sh";
open (BASH, ">$bash_file") || die "\n\nCould not open bash file: $bash_file\n\n";
foreach my $block (sort {$a <=> $b} keys %fasta_files){
  my $fasta_file = $fasta_files{$block}{outfile};
  my $job_id = $fasta_files{$block}{job_id};
  print BASH "/home/malachig/svn/hmmSplicer/splignRun.pl  --splign_bin_dir=$splign_bin_dir  --blast_bin_dir=$blast_bin_dir  --junction_fasta=$fasta_file  --genome_fasta=$genome_fasta  --splign_results_dir=$splign_results_dir  --job_id=$job_id\n";
}
close(BASH);


