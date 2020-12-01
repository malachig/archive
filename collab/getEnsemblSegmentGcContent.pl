#!/usr/bin/perl -w
#Written by Malachi Griffith
#The purpose of this script is to get data for a list of chromosome segments specified by coordinates

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
my $connect_type = '';
my $species = '';
my $segment_file = '';
my $output_file = '';

GetOptions ('ensembl_api_version=s'=>\$ensembl_api_version,
	    'ensembl_database=s'=>\$ensembl_database, 'ensembl_server=s'=>\$ensembl_server,
	    'ensembl_user=s'=>\$ensembl_user, 'ensembl_password=s'=>\$ensembl_password,
	    'connect_type=s'=>\$connect_type, 'species=s'=>\$species,
	    'segment_file=s'=>\$segment_file, 'output_file=s'=>\$output_file);

#Provide instruction to the user
print BLUE, "\n\nNOTE: Before using this script, make sure the correct API version is hard coded!!\n\n", RESET;
print BLUE, "\n\nBasic required parameters:", RESET;
print GREEN, "\n\tSpecify the correct EnsEMBL API version using: --ensembl_api_version (41, 42, etc.)", RESET;
print GREEN, "\n\tSpecify the desired connection type using: --connect_type (see below for examples)", RESET;
print GREEN, "\n\tSpecify the input file containing chromosome regions of interest", RESET;
print GREEN, "\n\t\tThis file should be of the form: chr1:1000-2000 (each line has: chromosome:start_position-end_position)", RESET;

print BLUE, "\n\nOption 1: Local database connection:", RESET;
print GREEN, "\n\tTo use a local EnsEMBL database specify the following:", RESET;
print GREEN, "\n\tUse: --connect_type=local", RESET;
print GREEN, "\n\tThe Ensembl Server, User and Password for access using: --ensembl_server, --ensembl_user and --ensembl_password", RESET;
print GREEN, "\n\tThe species using: --species (e.g. --species=Human or --species='Homo sapiens')", RESET;
print GREEN, "\n\t\tThe 'Core' database for this species for the EnsEMBL version you specified above will be used automatically", RESET;
print GREEN, "\n\nLocal example: getEnsemblSegmentGcContent.pl  --ensembl_api_version=35  --connect_type=local  --ensembl_server=ensembl01.bcgsc.ca  --ensembl_user=ensembl  --ensembl_password=ensembl  --species=Human  --segment_file=ExampleSegments_hg17.txt  --output_file=GC_content.txt\n", RESET;

print BLUE, "\n\nOption 2: Remote database connection over the web:", RESET;
print GREEN, "\n\tTo connect to EnsEMBL over the web specify the following:", RESET;
print GREEN, "\n\tUse: --connect_type=remote", RESET;
print GREEN, "\n\tThe species using: --species (e.g. --species=Human or --species='Homo sapiens')", RESET;
print GREEN, "\n\t\tThe 'Core' database for this species for the EnsEMBL version you specified above will be used automatically", RESET;
print GREEN, "\n\nWeb example: getEnsemblSegmentGcContent.pl  --ensembl_api_version=35  --connect_type=remote  --species=Human  --segment_file=ExampleSegments_hg17.txt  --output_file=GC_content.txt\n", RESET;

print BLUE, "\n\nOption 3: Legacy connection (required if EnsEMBL API version is X or earlier):", RESET;
print GREEN, "\n\tUse: --connect_type=legacy", RESET;

print GREEN, "\n\nlocal example: getEnsemblSegmentGcContent.pl  --ensembl_api_version=31  --connect_type=legacy  --ensembl_database=homo_sapiens_core_31_35d  --ensembl_server=ensembl01.bcgsc.ca  --ensembl_user=ensembl  --ensembl_password=ensembl  --segment_file=ExampleSegments_hg17.txt  --output_file=GC_content.txt", RESET;

print GREEN, "\n\nWeb example: getEnsemblSegmentGcContent.pl  --ensembl_api_version=31  --connect_type=legacy  --ensembl_database=homo_sapiens_core_31_35d  --ensembl_server=ensembldb.ensembl.org  --ensembl_user=anonymous  --ensembl_password=anything  --segment_file=ExampleSegments_hg17.txt  --output_file=GC_content.txt\n\n", RESET;


unless (($ensembl_api_version =~ /^\d+/) && $segment_file && $connect_type && $output_file){
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
my %segs;
open (SEGS, "$segment_file") || die "\nCould not open segment file: $segment_file\n\n";
my $seg_count = 0;
while (<SEGS>){
  $seg_count++;
  chomp($_);

  #e.g. chr9:118104113-118179376
  if ($_ =~ /^chr(\w+)\:(\d+)\-(\d+)/){

    $segs{$seg_count}{chr} = $1;
    $segs{$seg_count}{start} = $2;
    $segs{$seg_count}{end} = $3;

  }else{
    print RED, "\nLine format not understood:  $_\n\n", RESET;
    exit();
  }
}
close (SEGS);

print BLUE, "\n\nGetting sequence and GC content for regions in input list\n";

#3.) For each region, get the sequence calculate the GC content and add this column to the input file
open (OUT, ">$output_file") || die "\nCould not open output file: $output_file\n\n";
foreach my $seg_number (sort {$a <=> $b} keys %segs){

  my $chr = $segs{$seg_number}{chr};
  my $start = $segs{$seg_number}{start};;
  my $end = $segs{$seg_number}{end};;

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

  my $sequence = $slice->seq();

  # Finds G's & C's and replaces them with themselves ($1) at each letter (gi)
  #my $gc_count = $sequence =~ s/([GC])/$1/gi;

  my $uc_sequence = uc($sequence);
  my $G_count = ($uc_sequence =~ tr/G/G/);
  my $C_count = ($uc_sequence =~ tr/C/C/);
  my $gc_count = $G_count + $C_count;

  # rounds GC content calculation to single decimal place
  my $gc_content = sprintf "%.1f", (100*$gc_count / (length $sequence));

  print OUT "chr$chr:$start-$end\t$gc_content\n";
  print YELLOW, "$seg_number\tchr$chr:$start-$end\t$gc_content\n", RESET;

}

close (OUT);

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
