#!/usr/bin/perl -w
# Malachi Griffith

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;
use IO::File;
use Tie::File;
use Fcntl 'O_RDONLY';
use Benchmark;

my $splign_bin_dir = '';
my $blast_bin_dir = '';
my $junction_fasta = '';
my $genome_fasta = '';
my $splign_results_dir = '';
my $job_id = '';

GetOptions ('splign_bin_dir=s'=>\$splign_bin_dir, 'blast_bin_dir=s'=>\$blast_bin_dir, 'junction_fasta=s'=>\$junction_fasta, 'genome_fasta=s'=>\$genome_fasta, 
            'splign_results_dir=s'=>\$splign_results_dir, 'job_id=s'=>\$job_id);

if ($splign_bin_dir && $blast_bin_dir && $junction_fasta && $genome_fasta && $splign_results_dir && $job_id){
  #Parameters ok
}else{
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  print GREEN, "\nSpecify the path to the splign binary dir using: --splign_bin_dir", RESET;
  print GREEN, "\nSpecify the path to the splign binary dir using: --blast_bin_dir", RESET;
  print GREEN, "\nSpecify the path to a junction fasta file using: --junction_fasta", RESET;
  print GREEN, "\nSpecify the path to a genome fasta file using: --genome_fasta", RESET;
  print GREEN, "\nSpecify the path to a splign results dir using: --splign_results_dir", RESET;
  print GREEN, "\nSpecify a job id using: --job_id", RESET;
  print GREEN, "\n\nExample: splignRun.pl  --splign_bin_dir=/home/malachig/tools/splign/  --blast_bin_dir=/home/malachig/tools/BLAST2/blast2.2.18_x64/bin/  --junction_fasta=junctions.fa  --genome_fasta=/projects/malachig/sequence_databases/hg18_genome/all_human.fa  --splign_results_dir=/projects/alexa2/hmmSplicer/SplignValidation/hg18/  --job_id=junction_test\n\n", RESET;
  exit();
}

#Check dirs and files
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
unless (-e $junction_fasta){
  print RED, "\n\nJunction fasta: $junction_fasta not found\n\n", RESET;
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
my $temp_dir = "$splign_results_dir"."working/$job_id/";
my $splign_result_file = "$temp_dir"."splign.result.txt";
my $splign_summary_file = "$splign_results_dir"."results/$job_id"."_summary.txt";

#Global time variable
my $t1 = new Benchmark;

#Get the list of junction IDs being considered from the fasta file
my %hmm_junctions;
open (JID, $junction_fasta) || die "\n\nCould not open junction fasta: $junction_fasta\n\n";
while(<JID>){
  chomp($_);
  if ($_ =~ /^\>(.*)/){
    $hmm_junctions{$1}{splign_junctions} = 0;
    $hmm_junctions{$1}{splign_match} = 0;
  }
}
close(JID);

#my $debug = 1;
#unless ($debug){

#Set the temp dir using the splign results dir and the job id - make sure it does not already exist - then make it
if (-e $temp_dir){
  print RED, "\n\nTemp dir already exists: $temp_dir - aborting", RESET;
  exit();
}
mkdir($temp_dir);


#Symlink the fasta files to this dir
print BLUE, "\n\nSymlinking fasta files to working dir", RESET;
my $junction_fasta_name = '';
if ($junction_fasta =~ /fasta\/(\w+\.fa)/){
  $junction_fasta_name = $1;
}else{
  print RED, "\n\nCould not determine junction fasta name from: $junction_fasta", RESET;
  exit();
}
my $ln_cmd1 = "ln -s $junction_fasta $temp_dir$junction_fasta_name";
print BLUE, "\n\t$ln_cmd1", RESET;
system($ln_cmd1);

my $genome_fasta_base = '';
my $genome_fasta_name = '';
my $suffix = '';
if ($genome_fasta =~ /(.*\/)(\w+)(\.fa)/){
  $genome_fasta_base = $1;
  $genome_fasta_name = $2;
  $suffix = $3;
}else{
  print RED, "\n\nCould not determine genome fasta name from: $genome_fasta", RESET;
  exit();
}
my $ln_cmd2 = "ln -s $genome_fasta_base$genome_fasta_name$suffix $temp_dir$genome_fasta_name$suffix";
print BLUE, "\n\t$ln_cmd2", RESET;
system($ln_cmd2);


#Create LDS index in this dir
print BLUE, "\n\nCreating LDS index for splign", RESET;
my $mklds_cmd = "$splign_bin_dir"."splign -mklds $temp_dir";
print BLUE, "\n\t$mklds_cmd", RESET;
system($mklds_cmd);


#Symlink genome BLAST DB to this dir
print BLUE, "\n\nSymlinking genome BLAST DB to working dir", RESET;
my $ln_cmd3 = "ln -s $genome_fasta_base$genome_fasta_name$suffix".".nhr $temp_dir$genome_fasta_name$suffix".".nhr";
my $ln_cmd4 = "ln -s $genome_fasta_base$genome_fasta_name$suffix".".nin $temp_dir$genome_fasta_name$suffix".".nin";
my $ln_cmd5 = "ln -s $genome_fasta_base$genome_fasta_name$suffix".".nsd $temp_dir$genome_fasta_name$suffix".".nsd";
my $ln_cmd6 = "ln -s $genome_fasta_base$genome_fasta_name$suffix".".nsi $temp_dir$genome_fasta_name$suffix".".nsi";
my $ln_cmd7 = "ln -s $genome_fasta_base$genome_fasta_name$suffix".".nsq $temp_dir$genome_fasta_name$suffix".".nsq";
print BLUE, "\n\t$ln_cmd3", RESET; system($ln_cmd3);
print BLUE, "\n\t$ln_cmd4", RESET; system($ln_cmd4);
print BLUE, "\n\t$ln_cmd5", RESET; system($ln_cmd5);
print BLUE, "\n\t$ln_cmd6", RESET; system($ln_cmd6);
print BLUE, "\n\t$ln_cmd7", RESET; system($ln_cmd7);


#Create junction BLAST DB in this dir
print BLUE, "\n\nCreating BLAST DB for junction fasta", RESET;
chdir($temp_dir);
my $format_db_cmd = "$blast_bin_dir"."formatdb -pF -oT -i $temp_dir$junction_fasta_name";
print BLUE, "\n\t$format_db_cmd\n", RESET;
system($format_db_cmd);


#Run initial mapping ('compart') step
print BLUE, "\n\nRunning 'compart' step on junction and genome fastas", RESET;
my $compart_cmd = "$splign_bin_dir"."compart -qdb $temp_dir$junction_fasta_name -sdb $genome_fasta > $temp_dir"."cdna.compartments";
print BLUE, "\n\t$compart_cmd\n", RESET;
system($compart_cmd);


#Run actual 'splign' step
print BLUE, "\n\nRunning 'splign' step on comparts generated in previous step", RESET;
my $splign_cmd = "$splign_bin_dir"."splign -type est -ldsdir $temp_dir -comps $temp_dir"."cdna.compartments > $splign_result_file";
print BLUE, "\n\t$splign_cmd\n", RESET;
system($splign_cmd);

#}

#Parse splign results file - store one result per junction in the original junction fasta file
my $splign_ref = &parseSplignResult('-infile'=>$splign_result_file);


#Write summary for each junction to a results file in the master results dir
open (OUT, ">$splign_summary_file") || die "\n\nCould not open splign summary file for writing: $splign_summary_file\n\n";
print OUT "JID\tSplign_Match\tSplign_Junction_Count\n";
foreach my $jid (sort keys %hmm_junctions){
  print OUT "$jid\t$hmm_junctions{$jid}{splign_match}\t$hmm_junctions{$jid}{splign_junctions}\n";
}

#Remove the temp dir
my $rm_cmd = "rm -fr $temp_dir";
print BLUE, "\n\nCleaning up:\n$rm_cmd", RESET;
system($rm_cmd);


#print Dumper $splign_ref;
#print Dumper %hmm_junctions;

#Get time since last check
my $t2 = new Benchmark;
my $td = timediff($t2, $t1);
my $time_string = timestr($td);
my $seconds = "N";
my $minutes = "N";
my $hours = "N";
if ($time_string =~ /(\d+)\s+wallclock/){
  $seconds = $1;
  $minutes = sprintf("%.2f", ($seconds/60));
  $hours = sprintf("%.2f", ($minutes/60));
}
print BLUE, "\n\nProcessing time = $seconds seconds | $minutes minutes | $hours hours", RESET;

print "\n\n";

exit();


########################################################################################################################################
#Parse a splign result file ...                                                                                                        #
########################################################################################################################################
sub parseSplignResult{
  my %args = @_;
  my $infile = $args{'-infile'};
  my %splign;

  print BLUE, "\n\nParsing splign results file...", RESET;

  #Open the splign results file and dump it into an array... (one line per record)
  my @file_array;
  tie @file_array, 'Tie::File', $infile, mode => O_RDONLY;
  my $line_count = scalar(@file_array)-1;
  print BLUE, "\n\tFound $line_count lines in this file", RESET;

  #process two lines at a time.
  for (my $i = 0; $i < ($line_count-1); $i++){
    my $l1 = $file_array[$i];
    my $l2 = $file_array[$i+1];
    my @line1 = split("\t", $l1);
    my @line2 = split("\t", $l2);

    #Make sure both lines are part of the same splign block
    unless ($line1[0] == $line2[0]){
      next();
    }

    #Skip non-exon lines (these are missing basic info like coords)
    unless (($l1 =~ /exon/) && ($l2 =~ /exon/)){
      next();
    }

    #Get the strand of the splign alignment
    my $splign_id = $line1[0];
    my $splign_strand;
    if ($splign_id =~ /([\+\-])\d+/){
      $splign_strand = $1;
    }else{
      print RED, "\n\nCould not interpret strand from splign id: $splign_id\n\n", RESET;
      exit();
    }


    #If it is a negative strand record, swap the order of the two lines as well as the chr coords
    if ($splign_strand eq "-"){
      #print YELLOW, "\nSwapping", RESET;
      my @tmp1 = @line1;
      $line1[7] = $tmp1[8];
      $line1[8] = $tmp1[7];
      my @tmp2 = @line2;
      $line2[7] = $tmp2[8];
      $line2[8] = $tmp2[7];
      my @tmp3 = @line1;
      @line1 = @line2;
      @line2 = @tmp3;
    }
    my $jid = $line1[1];
    my $chr = $line1[2];
    my $pid1 = $line1[3];
    my $pid2 = $line2[3];
    my $al1 = $line1[4];
    my $al2 = $line2[4];
    my $left1 = $line1[7];
    my $left2 = $line2[7];
    my $right1 = $line1[8];
    my $right2 = $line2[8];

    my $new_j_id = "$chr:$right1-$left2($splign_strand)";
    #print YELLOW, "\n$new_j_id", RESET;

    #Store splign junctions
    $splign{$jid}{splign_junctions}++;
    if ($jid eq $new_j_id){
      $splign{$jid}{match} = 1;
    }

    #Add counts to main list of junctions
    $hmm_junctions{$jid}{splign_junctions}++;
    if ($jid eq $new_j_id){
      $hmm_junctions{$jid}{splign_match} = 1;
    }
  }
  untie @file_array;

  return(\%splign);
}


