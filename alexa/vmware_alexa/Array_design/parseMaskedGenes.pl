#!/usr/bin/perl -w

#Written by Malachi Griffith
#This script parses through all file in a directory where each file contains gene sequences that were masked by RepeatMasker
#These masked gene sequences will then be stored in ALEXA in a seperate table than the Gene table (for performance reasons)

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use DBI;
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
my $masked_dir = '';
my $expected_genes = '';
my $update_database = '';
my $logfile = '';

GetOptions ('database=s'=>\$database, 'server=s'=>\$server, 'user=s'=>\$user, 'password=s'=>\$password,
	    'masked_dir=s'=>\$masked_dir, 'expected_genes=i'=>\$expected_genes, 'update_database=s'=>\$update_database, 'logfile=s'=>\$logfile);

#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the database and server to query using: --database and --server", RESET;
print GREEN, "\n\tSpecify the user and password for access using: --user and --password", RESET;
print GREEN, "\n\tSpecify the path of the directory containing repeat masked fasta files using: --masked_dir", RESET;
print GREEN, "\n\tSpecify the number of gene records that are expected to be found in this directory using: --expected_genes=1234", RESET;
print GREEN, "\n\tOnce this script has been tested, use the --update_database option to allow a database update (ex. --update_database=yes)", RESET;
print GREEN, "\n\tAlso specify the name of a log file to be used to store a record of the database population using: --logfile", RESET;
print GREEN, "\n\nExample: parseMaskedGenes.pl  --database=ALEXA_hs_35_35h  --server=server_name  --user=user_name  --password=pwd  --masked_dir=/home/user/alexa/ALEXA_version/repeat_masking/ensembl_genes_masked/  --expected_genes=1234  --update_database=no  --logfile=/home/user/alexa/ALEXA_version/logs/database_population/parseMaskedGenes_LOG.txt\n\n", RESET;

#Make sure all options were specified
unless ($database && $server && $user && $password && $masked_dir && $expected_genes && $update_database && $logfile){
  print RED, "\nRequired input parameter(s) missing!\n\n", RESET;
  exit();
}

#Establish connection with the Alternative Splicing Expression database
my $alexa_dbh = &connectDB('-database'=>$database, '-server'=>$server, '-user'=>$user, '-password'=>$password);

unless ($update_database eq "yes"){
  print YELLOW, "\nOnce you are ready, use the --update_database=yes option to update the exon entries in the database\n\n", RESET;
}

open (LOG, ">$logfile") || die "\nCould not open logfile: $logfile\n\n";

#Print out the parameters supplied by the user to the logfile for future reference
print LOG "\nUser Specified the following options:\ndatabase = $database\nmasked_dir = $masked_dir\nexpected_genes = $expected_genes\nupdate_database = $update_database\nlogfile = $logfile\n\n";

#First get all sequences from the masked repeat file 
my $masked_genes_ref = &parseMaskedSequences('-dir'=>$masked_dir);
my $genes_found = keys %{$masked_genes_ref};

unless ($expected_genes == $genes_found){
  print RED, "\nDid not find the expected number of genes! (Expected: $expected_genes and Found: $genes_found)\n\n", RESET;
  close (LOG);
  exit();
}
print BLUE, "\nFound $genes_found masked gene sequences in the input files\n\n", RESET;
print LOG "\nFound $genes_found masked gene sequences in the input files\n\n";

#Make sure the max allowed packet size is large enough to allow large genes
&checkMaxPacketSize('-dbh'=>$alexa_dbh);

#Update masked gene entries in the database
if ($update_database eq "yes"){
  &updateDatabase('-dbh'=>$alexa_dbh);
}else{
  print YELLOW, "\n\nMasked gene sequences not parsed to database!! - Run script again using -update_database=yes option\n\n", RESET;
}

#Go through the masked sequences and count the number of masked bases
print BLUE, "\n\nSummarizing the number of masked and unmasked bases\n", RESET;
print LOG "\n\nSummarizing the number of masked and unmasked bases\n";
my $total_bases = 0;
my $masked_bases = 0;
foreach my $gene_id (keys %{$masked_genes_ref}){
  my $seq = $masked_genes_ref->{$gene_id}->{MASKED_SEQ};
  $total_bases += length($seq);
  $masked_bases += ($seq =~ tr/N/N/);
}
my $percent_masked = sprintf("%.3f",(($masked_bases/$total_bases)*100));
print BLUE, "\nFound a total of $masked_bases masked bases of $total_bases total bases ($percent_masked%)", RESET;
print BLUE, "\nIf this does not seem reasonable, perhaps the incorrect repeat library specied was specified\n\n", RESET;
print LOG "\nFound a total of $masked_bases masked bases of $total_bases total bases ($percent_masked%)";
print LOG "\nIf this does not seem reasonable, perhaps the incorrect repeat library specied was specified\n\n";

#Close database connection
$alexa_dbh->disconnect();

close (LOG);

exit();


#################################################################################################################################
#parseMaskedSequence - Get masked meta_transcript sequences from RepeatMasker masked sequence file                              #
#################################################################################################################################
sub parseMaskedSequences{
  my %args = @_;
  my $repeat_dir = $args{'-dir'};

  my %seqs;

  unless ($repeat_dir =~ /.*\/$/){
    $repeat_dir = "$repeat_dir"."/";
  }

  #First make sure the specified base path exists and is a directory
  unless (-e $repeat_dir && -d $repeat_dir){
    print RED, "\nSpecified directory: $repeat_dir does not appear valid!\n\n", RESET;
    close (LOG);
    exit();
  }

  #Get all the input files in the repeat masked result directory
  print BLUE, "\nSearching $repeat_dir for files\n", RESET;
  print LOG "\nSearching $repeat_dir for files\n";

  opendir(DIRHANDLE, "$repeat_dir") || die "\nCannot open directory: $repeat_dir\n\n";
  my @repeat_files = readdir(DIRHANDLE);

  #Order the files
  my %files;
  foreach my $file (@repeat_files){
    chomp($file);

    my $file_path = "$repeat_dir"."$file";

    #Skip directories within the specified directory
    if (-e $file_path && -d $file_path){
      print YELLOW, "\n\t$file_path  is a directory - skipping", RESET;
      next();
    }

    #Skip all files except masked files
    unless ($file_path =~ /\.masked/){
      next();
    }

    if ($file_path =~ /.*Gene_Block(\d+)\.masked/){
      my $file_number = $1;
      $files{$file_number}{file_path} = $file_path;
    }else{
      print RED, "\nCould not determine file number from filename:\n$file_path\n\n", RESET;
      close (LOG);
      $alexa_dbh->disconnect();
      exit();
    }
  }

  print BLUE, "\n\n", RESET;

  my $file_count = 0;
  foreach my $file (sort {$a <=> $b} keys %files){
    $file_count++;

    my $infile = $files{$file}{file_path};
    print BLUE, "\nProcessing: $infile", RESET;
    print LOG "\nProcessing: $infile";

    open (INFILE, "$infile") || die "\nCould not open input file of masked sequence $infile\n\n";

    #set input seperator to ">".
    my $save_input_seperator = $/;
    $/ = ">";

    my $first_line = 1;

    while (<INFILE>){

      #Skip first junk entry
      if ($first_line == 1){
	$first_line = 0;
	next();
      }

      my ($gene_id, $seq);

      #Seperate the header from the sequence 
      if ($_ =~ /^(\d+)(.*)/s){
	$gene_id = $1;
	$seq = $2;

	#clean-up sequence data by removing spaces, >,  and making uppercase
	$seq =~ s/[\s\>]//g;
	my $upper_seq = uc($seq);

	$seqs{$gene_id}{MASKED_SEQ} = $upper_seq;

      }else{
	print RED, "\nCould not interpret entry: $_\n\n", RESET;
	close (LOG);
	$alexa_dbh->disconnect();
	exit();
      }
    }

    #Reset input seperator
    $/ = $save_input_seperator;

    close (INFILE);
  }
  return(\%seqs);
}


#################################################################################################################################
#update exon entries in ALEXA to included masked exon sequence                                                                  #
#################################################################################################################################
sub updateDatabase{
  my %args = @_;
  my $dbh = $args{'-dbh'};

  my $gene_count = 0;
  my $gene_insert_count = 0;

  my @gene_ids = keys %{$masked_genes_ref};
  my $total_genes = @gene_ids;

  #First make sure the database is empty (if it has already been populate the user may have selected the wrong database to populate)
  print BLUE, "\n\nChecking status of MaskedGene table in $database before proceeding with population\n\n", RESET;
  print LOG "\n\nChecking status of MaskedGene table in $database before proceeding with population\n\n";

  my $sql_test = "select count(id) from MaskedGene;";
  my $sth_test = $dbh->prepare("$sql_test");
  $sth_test->execute();
  my $current_record_count = $sth_test->fetchrow_array();
  $sth_test->finish();

  if ($current_record_count > 0){
    print BLUE, "\nDatabase appears to have been populated already ($current_record_count records found in MaskedGene table) - Aborting\n\n", RESET;
    print LOG "\nDatabase appears to have been populated already ($current_record_count records found in MaskedGene table) - Aborting\n\n";
    close (LOG);
    $alexa_dbh->disconnect();
    exit();
  }else{
    print BLUE, "\nMaskedGene table appears empty - Proceeding with population\n\n", RESET;
    print LOG "\nMaskedGene table appears empty - Proceeding with population\n\n";
  }

  #Now get all the non-masked gene sequences to allow a basic length check
  print BLUE, "\nGet non-masked sequence from database to allow basic sanity check of gene length before/after masking\n\n", RESET;
  print LOG "\nGet non-masked sequence from database to allow basic sanity check of gene length before/after masking\n\n";

  my $gene_info_ref = &getGeneInfo('-dbh'=>$dbh, '-gene_ids'=>\@gene_ids, '-sequence'=>"no");

  print BLUE, "\n\nBegin database update\n", RESET;
  print LOG "\n\nBegin database update\n";

  foreach my $gene_id (sort {$a <=> $b} keys %{$masked_genes_ref}){
    $gene_count++;

    print BLUE, "\nProcessing $gene_count : $gene_id", RESET;
    print LOG "\nProcessing $gene_count : $gene_id";

    my $masked_seq = $masked_genes_ref->{$gene_id}->{MASKED_SEQ};

    unless (length($masked_seq) == $gene_info_ref->{$gene_id}->{seq_length}){
      print RED, "\nMasked and Unmasked sequence lengths do not agree", RESET;
      close (LOG);
      $alexa_dbh->disconnect();
      exit();
    }

    my $sql = "INSERT INTO MaskedGene (fk_Gene__id, sequence) VALUES ($gene_id,\'$masked_seq\');";

    #Actually insert the record
    my $sth = $dbh->prepare("$sql");
    $sth->execute();
    $sth->finish();
    $gene_insert_count++;
    print BLUE, "\n\t->New MaskedGene Record for gene: $gene_id", RESET;
    print LOG "\n\t->New MaskedGene Record for gene: $gene_id";
  }
  return();
}

