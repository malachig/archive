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
my $species = '';

GetOptions ('ensembl_api_version=s'=>\$ensembl_api_version,
	    'ensembl_database=s'=>\$ensembl_database, 'ensembl_server=s'=>\$ensembl_server,
	    'ensembl_user=s'=>\$ensembl_user, 'ensembl_password=s'=>\$ensembl_password,
	    'species=s'=>\$species);

#Provide instruction to the user
print BLUE, "\n\nNOTE: Before using this script, make sure the correct API version is hard coded!!\n\n", RESET;
print GREEN, "\n\nExample: ensembl_API_template.pl  --ensembl_api_version=54  --connect_type=local  --ensembl_server=ensembl01.bcgsc.ca  --ensembl_user=ensembl  --ensembl_password=ensembl  --species=Human\n", RESET;

unless (($ensembl_api_version =~ /^\d+/) && $ensembl_database && $ensembl_server && $ensembl_user && $ensembl_password && $species){
  print RED, "\nBasic option(s) missing or incorrect format\n\n", RESET;
  exit();
}


#**********************************************************************************************************
#IMPORTANT NOTE: You must have the correct Ensembl API installed locally AND bioperl 1.2 or greater!!
if ($ensembl_api_version =~ /^\d+/){
  if($ensembl_api_version eq "54"){
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


#1.) Establish connections to source EnsEMBL database - either locally, or remotely over the web
my $ensembl_api;

#A.) Using a local ensembl version
$ensembl_api = 'Bio::EnsEMBL::Registry';
$ensembl_api ->load_registry_from_db(-host=>$ensembl_server, -user=>$ensembl_user, -pass=>$ensembl_password);


exit();

