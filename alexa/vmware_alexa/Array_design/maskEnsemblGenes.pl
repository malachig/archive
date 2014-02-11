#!/usr/bin/perl -w

#Written by Malachi Griffith
#Grab all complete gene sequences from ALEXA
#Write these gene sequences in fasta format files, such that each file contains a limited number of sequences
#Create commands to use RepeatMasker on each of these files on the cluster or simply run the batch script on a single computer

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;

#ALEXA libraries
#When a script is initiated, use the full path of the script location at execution to add the perl module libraries to @INC
#This should allow this scripts to work regardless of the current working directory or the script location (where it was unpacked).
#The /utilities directory must remain in the same directory as this script but the entire code directory can be moved around
BEGIN {
  my $script_dir = &File::Basename::dirname($0);
  push (@INC, $script_dir);
}
use utilities::utility qw(:all);
use utilities::ALEXA_DB qw(:all);

#Initialize command line options
my $database = '';
my $server = '';
my $user = '';
my $password = '';
my $species = '';
my $repeat_masker = '';
my $temp_dir = '';
my $batch_file = '';
my $logfile = '';

GetOptions ('database=s'=>\$database, 'server=s'=>\$server, 'user=s'=>\$user, 'password=s'=>\$password,
	    'species=s'=>\$species, 'repeat_masker=s'=>\$repeat_masker, 'temp_dir=s'=>\$temp_dir,
	    'batch_file=s'=>\$batch_file, 'logfile=s'=>\$logfile);

#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the database and server to query using: --database and --server", RESET;
print GREEN, "\n\tSpecify the user and password for access using: --user and --password", RESET;
print GREEN, "\n\tSpecify a valid species for RepeatMasking using: --species", RESET;
print GREEN, "\n\t\tValid species names include: human, mouse, rattus, chicken, elegans, danio, etc.", RESET;
print GREEN, "\n\tSpecify the full path to the RepeatMasker binary you wish to use with: --repeat_masker", RESET;
print GREEN, "\n\tSpecify the full path to the output directory for gene fasta files and masked files using: --temp_dir", RESET;
print GREEN, "\n\tSpecify the name of the output batch file for cluster jobs using: --batch_file", RESET;
print GREEN, "\n\tAlso specify the name of a log file to be used to store a record of this script: --logfile", RESET;
print GREEN, "\n\nExample: maskEnsemblGenes.pl  --database=ALEXA_hs_35_35h  --server=server_name  --user=username  --password=pwd  --species='homo sapiens'  --repeat_masker=/usr/local/RepeatMasker/RepeatMasker3.1.6/RepeatMasker  --temp_dir=/home/user/alexa/ALEXA_version/repeat_masking/  --batch_file=/home/user/alexa/ALEXA_version/batch_scripts/repeat_masking/RepeatMaskAllGenes.sh  --logfile=/home/user/alexa/ALEXA_version/logs/maskEnsemblGenes_LOG.txt\n\n", RESET;

#Make sure all options were specified
unless ($database && $server && $user && $password && $species && $repeat_masker && $temp_dir && $batch_file && $logfile){
  print RED, "\nRequired input parameter(s) missing!\n\n", RESET;
  exit();
}

#Create a $fasta_dir called 'ensembl_genes_fasta' and a $masked_dir called 'ensembl_genes_masked'
my $fasta_dir = &createNewDir('-path'=>$temp_dir, '-new_dir_name'=>"ensembl_genes_fasta");
my $masked_dir = &createNewDir('-path'=>$temp_dir, '-new_dir_name'=>"ensembl_genes_masked");

#Establish connection with the Alternative Splicing Expression database
my $alexa_dbh = &connectDB('-database'=>$database, '-server'=>$server, '-user'=>$user, '-password'=>$password);

open (LOG, ">$logfile") || die "\nCould not open logfile: $logfile\n\n";

#Print out the parameters supplied by the user to the logfile for future reference
print LOG "\nUser Specified the following options:\ndatabase = $database\nspecies = $species\nrepeat_masker = $repeat_masker\ntemp_dir = $temp_dir\nbatch_file = $batch_file\nlogfile = $logfile\n\n";

my %gene;
my $gene_count = 0;

my @gene_ids = @{&getAllGenes ('-dbh'=>$alexa_dbh, '-gene_type'=>'All')};

my $unique_genes = @gene_ids;

print BLUE, "\nFound $unique_genes genes in $database\n\n", RESET;
print LOG "\nFound $unique_genes genes in $database\n\n";

my $block_size = 5000000;  #Cutoff for number of bases of sequence to put in a single file
my $block_count = 0;
my $current_block = 1;

my $genes_ref = &getGeneInfo ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids,'-sequence'=>"yes");

my $gene_info_count = keys %{$genes_ref};
print BLUE, "\n\nFound gene info for $gene_info_count genes\n\n", RESET;
print LOG "\n\nFound gene info for $gene_info_count genes\n\n";

#open the batch command output file
open (BATCH, ">$batch_file") || die "\nCould not open $batch_file\n\n";

#Open the first output file
my $out_file = "$fasta_dir"."Gene_Block$current_block";
open (OUTFILE, ">$out_file") || die "\nCould not open $out_file\n\n";

foreach my $gene_id (sort {$a <=> $b} keys %{$genes_ref}){
  $gene_count++;

  my $length = length($genes_ref->{$gene_id}->{sequence});
  $block_count += $length;

  print OUTFILE ">$gene_id\n$genes_ref->{$gene_id}->{sequence}\n";

  #Once the block size has been exceeded, close that file and start on the next one
  if ($block_count >= $block_size){
    close OUTFILE;

    print BATCH "$repeat_masker -species '$species' -dir $masked_dir $out_file\n";
    print BLUE, "$repeat_masker -species '$species' -dir $masked_dir $out_file\n", RESET;

    $current_block++;

    $out_file = "$fasta_dir"."Gene_Block$current_block";
    open (OUTFILE, ">$out_file") || die "\nCould not open $out_file";

    $block_count = 0;
  }

}
close (OUTFILE);

print BATCH "$repeat_masker -species '$species' -dir $masked_dir $out_file\n";
print BLUE, "$repeat_masker -species '$species' -dir $masked_dir $out_file\n", RESET;

print BLUE, "\n\nWrote $gene_count gene sequences into $current_block files in $fasta_dir\n\n", RESET;
print LOG "\n\nWrote $gene_count gene sequences into $current_block files in $fasta_dir\n\n";

close (BATCH);
close (LOG);

#Close database connection
$alexa_dbh->disconnect();

exit();
