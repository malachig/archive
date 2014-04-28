#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#The purpose of this script is to create a database random sequence (by sampling without replacement from the exonic content of the human genome)

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;
use Benchmark;
use BerkeleyDB;

#Load the ALEXA modules
BEGIN {
  use Cwd 'abs_path';
  if (abs_path($0) =~ /(.*)\/.*\/.*\.pl/){
    push (@INC, $1);
  }
}
use utilities::utility qw(:all);
use utilities::ALEXA_DB qw(:all);
use utilities::Shuffle qw(shuffle);

my $database = '';
my $server = '';
my $user = '';
my $password = '';
my $outdir = '';
my $logfile = '';

GetOptions ('database=s'=>\$database, 'server=s'=>\$server, 'user=s'=>\$user, 'password=s'=>\$password, 
            'outdir=s'=>\$outdir, 'logfile=s'=>\$logfile);


#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the TARGET Database and Server to query using: --database and --server", RESET;
print GREEN, "\n\tSpecify the User and Password for access using: --user and --password\n", RESET;
print GREEN, "\n\tSpecify the path to be used for output files using: --outdir", RESET;
print GREEN, "\n\tAlso specify the name of a log file to be used to store a record of the design process using: --logfile", RESET;

print GREEN, "\n\nExample: createRandomSequenceDatabase.pl  --database=ALEXA_hs_49_36k  --server=jango.bcgsc.ca  --user=viewer  --password=viewer  --outdir=/projects/malachig/sequence_databases/hs_49_36k/ensembl_random_exonic_hs_49_36k/temp/  --logfile=/projects/malachig/solexa/logs/alternativeExpressionDatabase/createRandomSequenceDatabase/createRandomSequenceDatabase_LOG.txt\n\n", RESET;

#Make sure all options were specified
unless ($database && $server && $user && $password && $outdir && $logfile){
  print RED, "\nRequired input parameter(s) missing!\n\n", RESET;
  exit();
}

#Global time and process variables
my $t1 = new Benchmark;
my $pid = $$;

$outdir = &checkDir('-dir'=>$outdir, '-clear'=>"no");

open (LOG, ">$logfile") || die "\nCould not open logfile: $logfile\n\n";

#Set cache size for Berkley DB objects
my $cachesize = 256000000;

#Hash to store actual complete exonic content of human genome
my $seq_file_name = "$outdir"."ExonicSeq.btree";
my $seq_rm_cmd = "rm -f $seq_file_name";
system($seq_rm_cmd);
my %seq;
print YELLOW, "\nCreating binary tree file: $seq_file_name\n", RESET;
tie(%seq, 'BerkeleyDB::Btree', -Cachesize => $cachesize, -Filename=> $seq_file_name , -Flags => DB_CREATE) or die "can't open file $seq_file_name: $! $BerkeleyDB::Error\n";
my $seq_ref = \%seq;

#1.) Get gene models from ALEXA for all genes - get completed masked sequence for each gene locus

#Establish connection with the Alternative Splicing Expression database
my $alexa_dbh = &connectDB('-database'=>$database, '-server'=>$server, '-user'=>$user, '-password'=>$password);
my $genes_ref;

my @gene_ids = @{&getAllGenes ('-dbh'=>$alexa_dbh, '-gene_type'=>'All')};

$| = 1; print BLUE, "\n\n1.) Getting basic gene info", RESET; $| = 0;
print LOG "\n\n1.) Getting basic gene info";
&getBasicGeneInfo ('-gene_ids'=> \@gene_ids);

#Close database connection
$alexa_dbh->disconnect();

#How many bases in the seq object generated
my $base_total = keys %{$seq_ref};
print BLUE, "\n\nStored a total of $base_total exonic bases\n", RESET;
print LOG "\n\nStored a total of $base_total exonic bases\n";

#Now go through the position-base hash, retrieve a random key, get the base for this key, add the base to a new random seq hash, and delete the corresponding key from the source hash
my %seq_base_counts;
my %rseq_base_counts;

#Print out random seq to a fasta file
my $fasta_file = "$outdir"."randomExonic.fa"; 
print BLUE, "\n\nPrinting random exonic sequence to: $fasta_file", RESET;
print LOG "\n\nPrinting random exonic sequence to: $fasta_file";

open (FASTA, ">$fasta_file") || die "\nCould not open random seq fasta file: $fasta_file\n\n";
print FASTA ">RandomExonicSeq\n";

my $base_countdown = $base_total;
my $base_counter = 0;
foreach my $i (shuffle keys %{$seq_ref}){

  $base_counter++;
  $base_countdown--;

  if ($base_counter == 100000){
    $base_counter = 0;
    print BLUE, "\nBases remaining: $base_countdown", RESET;
    &processUpdate('-display_option'=>1);
  }

  $rseq_base_counts{$seq_ref->{$i}}++;
  print FASTA "$seq_ref->{$i}";
  #print "$seq_ref->{$i}";
}
close(FASTA);

#Summarize the base composition of the input exonic base content and the random output (should be the same!!)
print BLUE, "\n\nBase composition for input exonic content of human genome ($base_total bp)\n", RESET;
foreach my $base (sort {$a cmp $b} keys %seq_base_counts){
  my $count = $seq_base_counts{$base};
  my $percent = sprintf("%.1f",($count/$base_total)*100);
  print BLUE, "$base: $count ($percent%)\t", RESET;
  print LOG "$base: $count ($percent%)\t";
}

print BLUE, "\n\nBase composition for random exonic content generated from human genome ($base_total bp)\n", RESET;
foreach my $base (sort {$a cmp $b} keys %rseq_base_counts){
  my $count = $rseq_base_counts{$base};
  my $percent = sprintf("%.1f",($count/$base_total)*100);
  print BLUE, "$base: $count ($percent%)\t", RESET;
  print LOG "$base: $count ($percent%)\t";
}



#Clean-up temp files
untie(%seq);
system($seq_rm_cmd);
close(LOG);

exit();


################################################################################################
#Get gene models, masked gene sequence, intron content, etc.                                   #
################################################################################################
sub getBasicGeneInfo{
  my %args = @_;
  my @gene_ids = @{$args{'-gene_ids'}};

  #Get basic info about each gene
  $| = 1; print BLUE, "\n\nGetting gene models", RESET; $| = 0; 
  print LOG "\nGetting gene models";

  my $storable_name = "$database"."_AllGenes_GeneInfo_WithSeq.storable";
  my $genes_ref = &getGeneInfo ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids, '-sequence'=>"yes", '-storable'=>$storable_name);
  
  #Get exon content for each gene
  $| = 1; print BLUE, "\n\nGetting exon content", RESET; $| = 0; 
  print LOG "\n\nGetting exon content";
  $storable_name = "$database"."_AllGenes_ExonContent.storable";
  my $gene_exon_content_ref = &getExonContent ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids, '-storable'=>$storable_name);

  print BLUE, "\n\nStoring a hash of all exonic content for the entire genome", RESET;
  my $grand_base_count = 0;

  my $gene_countdown = keys %{$genes_ref};
  my $gene_countup = 0;
  my $gene_counter = 0;
  GENE:foreach my $gene_id (keys %{$genes_ref}){
    $gene_countdown--;
    $gene_counter++;
    $gene_countup++;

    if ($gene_counter == 100){
      $gene_counter = 0;
      print BLUE, "\nGenes remaining: $gene_countdown\tBases stored: $grand_base_count", RESET;
      &processUpdate('-display_option'=>1);
    }

    #DEBUG
    #if ($gene_countup == 1000){
    #  last GENE;
    #}

    my $exon_content_ref = $gene_exon_content_ref->{$gene_id}->{exon_content};
    my $gene_seq = $genes_ref->{$gene_id}->{sequence};

    foreach my $ec (keys %{$exon_content_ref}){

      my $start = $exon_content_ref->{$ec}->{start};
      my $end = $exon_content_ref->{$ec}->{end};
      my $size = ($end-$start)+1;

      #Get the sequence for this exon content block
      my $ec_seq = substr($gene_seq, $start, $size);

      #Split this sequence into an array (one base per array element)
      my @seq = split("", $ec_seq);

      #print "\n$ec_seq";

      #Transfer sequence to one massive indexed object of the form: %hash{$pos} = $base 
      foreach my $base (@seq){
        $grand_base_count++;
        $seq_base_counts{$base}++;
        $seq_ref->{$grand_base_count} = $base;
      }
    }
  }

  return();
}


##################################################################################################
#Display compute time and memory usage                                                           #
##################################################################################################
sub processUpdate{
  my %args = @_;
  my $display_option = $args{'-display_option'};

  unless($display_option){
    $display_option = 0;
  }

  #Get memory usage in k and percent 
  my $ps_cmd = "ps -p $pid -o \"rss pmem\"";
  my $result = `$ps_cmd`;
  #print Dumper $result;
  my $memory_k = "Unknown";
  my $memory_m = "Unknown";
  my $memory_g = "Unknown";
  my $pmem = "Unknown";
  if ($result =~ /(\d+)\s+(\d+\.\d+)/m){
    $memory_k = $1;
    $memory_m = $memory_k/1024;
    $memory_g = sprintf("%.2f", ($memory_m/1024));
    $pmem = $2;
  } 

  #Get time since last check
  my $t2 = new Benchmark;
  my $td = timediff($t2, $t1);
  my $time_string = timestr($td);
  my $seconds = "Unknown";
  if ($time_string =~ /(\d+)\s+wallclock/){
    $seconds = $1;
  }

  $| = 1;
  if ($display_option == 1){
    print CYAN, "\n\tPROCESS: Time since last update: $seconds seconds\tCurrent memory usage: $memory_g Gb ($pmem%)\n", RESET;
  }elsif($display_option == 2){
    print CYAN, "\tCurrent memory usage: $memory_g Gb ($pmem%)", RESET;
  }else{
    print CYAN, "\nPROCESS: Time since last update: $seconds seconds\tCurrent memory usage: $memory_g Gb ($pmem%)\n", RESET;
  }
  $| = 0;
  $t1 = $t2;
  return();
}




