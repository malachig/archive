#!/usr/bin/perl -w
#Written by Malachi Griffith
#The purpose of this script is to get data for a chromosome segment specified by coordinates

#Make sure you select a version of the ensembl API that matches the database and genome build you wish to use!!
#For example to use homo_sapiens_core_35_35h  (Ensembl version 35, genome version 35h which equals hg17)

#NOTE: Before running this script you must have access to a local copy of an EnsEMBL database
#To see what EnsEMBL databases are available, log into the local ensembl server with mysql as follows:
#   mysql -h ensembl01.bcgsc.ca -u ensembl -pensembl
#Then use the command 'show databases'

use DBI;
use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);

#Initialize command line options
my $ensembl_api_version = ''; #Version of EnsEMBL to use
my $ensembl_database = '';
my $ensembl_server = '';
my $ensembl_user = '';
my $ensembl_password = '';
my $segment_file = '';
my $connect_type = '';
my $species = '';
my $fasta_file = '';
my $gene_file = '';

GetOptions ('ensembl_api_version=s'=>\$ensembl_api_version,
	    'ensembl_database=s'=>\$ensembl_database, 'ensembl_server=s'=>\$ensembl_server,
	    'ensembl_user=s'=>\$ensembl_user, 'ensembl_password=s'=>\$ensembl_password,
	    'connect_type=s'=>\$connect_type, 'species=s'=>\$species,
	    'segment_file=s'=>\$segment_file, 'fasta_file=s'=>\$fasta_file, 'gene_file=s'=>\$gene_file);

#Provide instruction to the user
print BLUE, "\n\nNOTE: Before using this script, make sure the correct API version is hard coded!!\n\n", RESET;
print GREEN, "\n\tThis script takes a list of coordinates as input, gets the chromosomal sequences and writes them to a fasta file", RESET;
print GREEN, "\n\tA list of all genes overlaping these coordinates is also written", RESET;
print BLUE, "\n\nBasic required parameters:", RESET;
print GREEN, "\n\tSpecify the correct EnsEMBL API version using: --ensembl_api_version (41, 42, etc.)", RESET;
print GREEN, "\n\tSpecify the desired connection type using: --connect_type (see below for examples)", RESET;
print GREEN, "\n\tSpecify the input file containing chromosome regions of interest", RESET;
print GREEN, "\n\t\tThis file should be tab-delimited (no header line) with each line containing: chromosome  start_position  end_position", RESET;
print GREEN, "\n\tSpecify output fasta file using: --fasta_file", RESET;
print GREEN, "\n\tSpecify output gene list file using: --gene_file", RESET;

print BLUE, "\n\nOption 1: Local database connection:", RESET;
print GREEN, "\n\tTo use a local EnsEMBL database specify the following:", RESET;
print GREEN, "\n\tUse: --connect_type=local", RESET;
print GREEN, "\n\tThe Ensembl Server, User and Password for access using: --ensembl_server, --ensembl_user and --ensembl_password", RESET;
print GREEN, "\n\tThe species using: --species (e.g. --species=Human or --species='Homo sapiens')", RESET;
print GREEN, "\n\t\tThe 'Core' database for this species for the EnsEMBL version you specified above will be used automatically", RESET;
print GREEN, "\n\nLocal example: getEnsemblChromosomeSegmentData.pl  --ensembl_api_version=35  --connect_type=local  --ensembl_server=ensembl01.bcgsc.ca  --ensembl_user=ensembl  --ensembl_password=ensembl  --species=Human  --segment_file=ExampleSegments_hg17.txt  --fasta_file=segments.fa  --gene_file=gene_list.txt\n", RESET;

print BLUE, "\n\nOption 2: Remote database connection over the web:", RESET;
print GREEN, "\n\tTo connect to EnsEMBL over the web specify the following:", RESET;
print GREEN, "\n\tUse: --connect_type=remote", RESET;
print GREEN, "\n\tThe species using: --species (e.g. --species=Human or --species='Homo sapiens')", RESET;
print GREEN, "\n\t\tThe 'Core' database for this species for the EnsEMBL version you specified above will be used automatically", RESET;
print GREEN, "\n\nWeb example: getEnsemblChromosomeSegmentData.pl  --ensembl_api_version=35  --connect_type=remote  --species=Human  --segment_file=ExampleSegments_hg17.txt  --fasta_file=segments.fa  --gene_file=gene_list.txt\n", RESET;

print BLUE, "\n\nOption 3: Legacy connection (required if EnsEMBL API version is X or earlier):", RESET;
print GREEN, "\n\tUse: --connect_type=legacy", RESET;

print GREEN, "\n\nlocal example: getEnsemblChromosomeSegmentData.pl  --ensembl_api_version=31  --connect_type=legacy  --ensembl_database=homo_sapiens_core_31_35d  --ensembl_server=ensembl01.bcgsc.ca  --ensembl_user=ensembl  --ensembl_password=ensembl  --segment_file=ExampleSegments_hg17.txt", RESET;

print GREEN, "\n\nWeb example: getEnsemblChromosomeSegmentData.pl  --ensembl_api_version=31  --connect_type=legacy  --ensembl_database=homo_sapiens_core_31_35d  --ensembl_server=ensembldb.ensembl.org  --ensembl_user=anonymous  --ensembl_password=anything  --segment_file=ExampleSegments_hg17.txt  --fasta_file=segments.fa  --gene_file=gene_list.txt\n\n", RESET;


unless (($ensembl_api_version =~ /^\d+/) && $segment_file && $connect_type && $fasta_file && $gene_file){
  print RED, "\nBasic option(s) missing or incorrect format\n\n", RESET;
  exit();
}
chomp($connect_type);

#Depending on the connection type specified by the user, connect to the database
#This code is required because the style of connection differs for local vs. remote and has changed with newer version of the API (thus the 'legacy' option)
if ($connect_type =~ /^local$/i){
  unless ($ensembl_server && $ensembl_user && $ensembl_password && $species){
    print RED, "'LOCAL' connections require the following parameters to be specified: --ensembl_server, --ensembl_user, --ensembl_password, and --species\n\n", RESET;
    exit();
  }
  print BLUE, "\nAttempting to connect to the specified LOCAL EnsEMBL database (Version $ensembl_api_version)\n\n", RESET;
}elsif ($connect_type =~ /^remote$/i){
  unless ($species){
    print RED, "'REMOTE' connections require the following parameters to be specified: --species\n\n", RESET;
    exit();
  }
  print BLUE, "\nAttempting to connect to a REMOTE EnsEMBL database (Version $ensembl_api_version) over the web - this will be much slower\n\n", RESET;
}elsif ($connect_type =~ /^legacy$/i){
  unless ($ensembl_database && $ensembl_server && $ensembl_user && $ensembl_password){
    print RED, "'LEGACY' connections require the following parameters to be specified: --ensembl_database, --ensembl_server, --ensembl_user, and --ensembl_password\n\n", RESET;
    exit();
  }
  print BLUE, "\nAttempting a legacy connection to EnsEMBL (Version $ensembl_api_version)\n\n", RESET;
}else{
  print RED, "\nSpecified '--connect_type' was not understood! - see examples above!\n\n", RESET;
  exit();
}


#**********************************************************************************************************
#IMPORTANT NOTE: You must have the correct Ensembl API installed locally AND bioperl 1.2 or greater!!
if ($ensembl_api_version =~ /^\d+/){

  if ($ensembl_api_version eq "27"){
    unshift(@INC, "/home/malachig/perl/ensembl_27_perl_API/ensembl/modules");
  }elsif ($ensembl_api_version eq "28"){
    unshift(@INC, "/home/malachig/perl/ensembl_28_perl_API/ensembl/modules");
  }elsif($ensembl_api_version eq "29"){
    unshift(@INC, "/home/malachig/perl/ensembl_29_perl_API/ensembl/modules");
  }elsif($ensembl_api_version eq "30"){
    unshift(@INC, "/home/malachig/perl/ensembl_30_perl_API/ensembl/modules");
  }elsif($ensembl_api_version eq "31"){
    unshift(@INC, "/home/malachig/perl/ensembl_31_perl_API/ensembl/modules");
  }elsif($ensembl_api_version eq "32"){
    unshift(@INC, "/home/malachig/perl/ensembl_32_perl_API/ensembl/modules");
  }elsif($ensembl_api_version eq "33"){
    unshift(@INC, "/home/malachig/perl/ensembl_33_perl_API/ensembl/modules");
  }elsif($ensembl_api_version eq "34"){
    unshift(@INC, "/home/malachig/perl/ensembl_34_perl_API/ensembl/modules");
  }elsif($ensembl_api_version eq "35"){
    unshift(@INC, "/home/malachig/perl/ensembl_35_perl_API/ensembl/modules");
  }elsif($ensembl_api_version eq "36"){
    unshift(@INC, "/home/malachig/perl/ensembl_36_perl_API/ensembl/modules");
  }elsif($ensembl_api_version eq "37"){
    unshift(@INC, "/home/malachig/perl/ensembl_37_perl_API/ensembl/modules");
  }elsif($ensembl_api_version eq "38"){
    unshift(@INC, "/home/malachig/perl/ensembl_38_perl_API/ensembl/modules");
  }elsif($ensembl_api_version eq "39"){
    unshift(@INC, "/home/malachig/perl/ensembl_39_perl_API/ensembl/modules");
  }elsif($ensembl_api_version eq "40"){
    unshift(@INC, "/home/malachig/perl/ensembl_40_perl_API/ensembl/modules");
  }elsif($ensembl_api_version eq "41"){
    unshift(@INC, "/home/malachig/perl/ensembl_41_perl_API/ensembl/modules");
  }elsif($ensembl_api_version eq "42"){
    unshift(@INC, "/home/malachig/perl/ensembl_42_perl_API/ensembl/modules");
  }elsif($ensembl_api_version eq "43"){
    unshift(@INC, "/home/malachig/perl/ensembl_43_perl_API/ensembl/modules");
  }elsif($ensembl_api_version eq "44"){
    unshift(@INC, "/home/malachig/perl/ensembl_44_perl_API/ensembl/modules");
  }elsif($ensembl_api_version eq "45"){
    unshift(@INC, "/home/malachig/perl/ensembl_45_perl_API/ensembl/modules");
  }elsif($ensembl_api_version eq "46"){
    unshift(@INC, "/home/malachig/perl/ensembl_46_perl_API/ensembl/modules");
  }elsif($ensembl_api_version eq "47"){
    unshift(@INC, "/home/malachig/perl/ensembl_47_perl_API/ensembl/modules");
  }elsif($ensembl_api_version eq "48"){
    unshift(@INC, "/home/malachig/perl/ensembl_48_perl_API/ensembl/modules");
  }elsif($ensembl_api_version eq "49"){
    unshift(@INC, "/home/malachig/perl/ensembl_49_perl_API/ensembl/modules");
  }elsif($ensembl_api_version eq "50"){
    unshift(@INC, "/home/malachig/perl/ensembl_50_perl_API/ensembl/modules");
  }elsif($ensembl_api_version eq "51"){
    unshift(@INC, "/home/malachig/perl/ensembl_51_perl_API/ensembl/modules");
  }elsif($ensembl_api_version eq "52"){
    unshift(@INC, "/home/malachig/perl/ensembl_52_perl_API/ensembl/modules");
  }elsif($ensembl_api_version eq "53"){
    unshift(@INC, "/home/malachig/perl/ensembl_53_perl_API/ensembl/modules");
  }elsif($ensembl_api_version eq "54"){
    unshift(@INC, "/home/malachig/perl/ensembl_54_perl_API/ensembl/modules");
  }elsif($ensembl_api_version eq "55"){
    unshift(@INC, "/home/malachig/perl/ensembl_55_perl_API/ensembl/modules");
  }else{
    print RED, "\nEnsEMBL API version: $ensembl_api_version is not defined, modify script before proceeding\n\n", RESET;
    exit();
  }
}else{
  print RED, "\nEnsEMBL API version format: $ensembl_api_version not understood!\n\n", RESET;
  exit();
}
use lib "/home/malachig/perl/bioperl-1.4";    #Bioperl
#*********************************************************************************************************
require Bio::EnsEMBL::DBSQL::DBAdaptor; #Used for local connections
require Bio::EnsEMBL::Registry;  #Use for remote connections over the web


#Note: If the user specifies an API version of 33 or earlier a 'legacy' connection will be required
if ($connect_type =~ /^local$|^remote$/i){
  if ($ensembl_api_version <= 33){
    print RED, "\nEnsEMBL API version earlier than v34 require a legacy connection! - See examples\n\n", RESET;
    exit();
  }
}

#1.) Establish connections to source EnsEMBL database - either locally, or remotely over the web
my $ensembl_api;

if ($connect_type =~ /^local$/i){
  #A.) Using a local ensembl version

  $ensembl_api = 'Bio::EnsEMBL::Registry';
  $ensembl_api ->load_registry_from_db(-host=>$ensembl_server, -user=>$ensembl_user, -pass=>$ensembl_password);

}elsif ($connect_type =~ /^remote$/i){
  #B.) Connecting over the web
  $ensembl_api = 'Bio::EnsEMBL::Registry';
  $ensembl_api ->load_registry_from_db(-host=>'ensembldb.ensembl.org', -user=>'anonymous');

}elsif ($connect_type =~ /^legacy$/i){
  #C.) Legacy connection - Old connection style could be used for local OR web connections - replaced by 'registry method'
  if ($ensembl_user eq "anonymous"){
    $ensembl_password = '';
  }
  $ensembl_api = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => $ensembl_server,
  						    -user => $ensembl_user,
  						    -dbname => $ensembl_database,
						    -pass => $ensembl_password);
}


#2.) Get the regions of interest from the input file

#Check the input file
unless (-e $segment_file){
  print RED, "\nCould not find input file: $segment_file\n\n", RESET;
  exit();
}
my $result = `cat -A $segment_file`;
if ($result =~ /\^M\$/g){
  print RED, "\n\nThis file appears to be from windows... you should run dos2unix on it before proceeding!\n\n", RESET;
  exit();
}

my %segs;
open (SEGS, "$segment_file") || die "\nCould not open segment file: $segment_file\n\n";
my $seg_count = 0;
while (<SEGS>){
  $seg_count++;
  chomp($_);
  if ($_ =~ /\d+/){
    my @line = split ("\t", $_);
    $segs{$seg_count}{chr} = $line[0];
    $segs{$seg_count}{start} = $line[1];
    $segs{$seg_count}{end} = $line[2];
  }
}
close (SEGS);

my %data_lines;

#3.) For each region, get the sequence and write to a fasta file
open (FASTA, ">$fasta_file") || die "\nCould not open output fasta file: $fasta_file\n\n";
my $count = 0;
foreach my $seg_number (sort {$a <=> $b} keys %segs){

  my $chr = $segs{$seg_number}{chr};
  my $start = $segs{$seg_number}{start};;
  my $end = $segs{$seg_number}{end};;

  print BLUE, "\nGetting sequence and genes for region: chr$chr:$start-$end\n", RESET;

  #A Slice object represents a continuous region of a genome
  my $slice_adaptor;
  if ($connect_type =~ /^local$|^remote$/i){
    $slice_adaptor = $ensembl_api->get_adaptor($species, 'Core', 'Slice');

    #Check if a slice adaptor was actually found
    unless ($slice_adaptor){
      print RED, "\nCould not get slice for this species - check EnsEMBL version, species, available databases on server, etc.\n\n", RESET;
      exit();
    }

  }elsif ($connect_type =~ /^legacy$/i){
    $slice_adaptor = $ensembl_api->get_SliceAdaptor();
  }

  #Retrieve a slice with respect to a gene, with a specified flanking sequence on either side
  my $slice = $slice_adaptor->fetch_by_region('chromosome', $chr, $start, $end);

  unless(defined($slice)){
    print RED, "\nCould not retrieve a sequence slice for chr$chr:$start-$end - Are you sure that is a valid chromosome/coordinate??\n\n", RESET;
    exit();
  }

  my $sequence = $slice->seq();
  print FASTA ">$chr:$start-$end\n$sequence\n";

  my @genes = @{$slice->get_all_Genes()};

  foreach my $gene (@genes){
    $count++;
    my %gene_info = %{&feature_details('-feature_obj'=>$gene)};
    my $stable_id = $gene->stable_id();

    my $test_name = $gene->external_name();
    unless ($gene->is_known()){
      $test_name = "Unknown";
    }

    my $gene_description = $gene->description();
    unless ($gene_description){
      $gene_description = "NA";
    }

    #Get chromosome coordinates for this gene
    my $new_gene = $gene->transform('chromosome');

    my $chr_test = $new_gene->slice->seq_region_name();
    my $chr_start_test = $new_gene->start();
    my $chr_end_test = $new_gene->end();
    my $strand_test = $new_gene->strand();
  
    #Coordinates relative to input chromosome segments
    #print YELLOW, "\n\tGene: $stable_id ($test_name)\t$gene_info{$stable_id}{chromosome} ($gene_info{$stable_id}{strand}):$gene_info{$stable_id}{start} - $gene_info{$stable_id}{end})", RESET;
    #$data_lines{$count}{line} = "$stable_id\t$test_name\t$gene_info{$stable_id}{chromosome}\t$gene_info{$stable_id}{strand}\t$gene_info{$stable_id}{start}\t$gene_info{$stable_id}{end}\t$gene_description\n";
    
    #Coordinates relative to whole chromosome
    print YELLOW, "\n\tGene: $stable_id ($test_name)\t$chr_test ($strand_test):$chr_start_test - $chr_end_test)", RESET;
    $data_lines{$count}{line} = "$stable_id\t$test_name\t$chr_test\t$strand_test\t$chr_start_test\t$chr_end_test\t$gene_description\n";

  }
  print "\n";
}
close (FASTA);

open (GENE, ">$gene_file") || die "\nCould not open output gene file: $gene_file\n\n";
foreach my $record (sort {$a <=> $b} keys %data_lines){
  print GENE "$data_lines{$record}{line}";
}

close (GENE);


exit();


########################################################################
#For any basic feature, get basic details and return as hash keyed on  #
#the stable id                                                         #
########################################################################
sub feature_details {
  my %args = @_;
  my $f = $args{'-feature_obj'};

  my %feature;

  my $stable_id = $f->stable_id();

  $feature{$stable_id}{chromosome} = $f->slice->seq_region_name();
  $feature{$stable_id}{start} = $f->start();
  $feature{$stable_id}{end} = $f->end();
  $feature{$stable_id}{strand} = $f->strand();

  return (\%feature);
}
