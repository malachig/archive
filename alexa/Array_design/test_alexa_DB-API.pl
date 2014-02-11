#!/usr/bin/perl -w
#Written by Malachi Griffith
#The purpose of this script is to illustrate simple code blocks for testing various methods in the ALEXA database API

use strict;
use Data::Dumper;
use Term::ANSIColor qw(:constants);
use Getopt::Long;
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
my $alexa_database = '';
my $alexa_server = '';
my $alexa_user = '';
my $alexa_password = '';

GetOptions ('alexa_database=s'=>\$alexa_database, 'alexa_server=s'=>\$alexa_server, 'alexa_user=s'=>\$alexa_user, 'alexa_password=s'=>\$alexa_password);

#Provide instruction to the user
print GREEN, "\n\nNOTE: Before using this script, uncomment the desired test or create a new one\n\n", RESET;
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the Database and Server to query using: --alexa_database and --alexa_server", RESET;
print GREEN, "\n\tSpecify the User and Password for access using: --alexa_user and --alexa_password\n", RESET;
print GREEN, "\n\nExample: test_alexa_DB-API.pl  --alexa_database=ALEXA_hs_41_36c  --alexa_server=server_name  --alexa_user=user  --alexa_password=pwd\n\n", RESET;

unless ($alexa_database && $alexa_server && $alexa_user && $alexa_password){
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  exit();
}

my $alexa_dbh = &connectDB('-database'=>$alexa_database, '-server'=>$alexa_server, '-user'=>$alexa_user, '-password'=>$alexa_password);

#1.) Test getAllGenes()
#print Dumper $alexa_dbh;
#my @gene_ids = @{getAllGenes ('-dbh'=>$alexa_dbh, '-gene_type'=>'Non-pseudo')};
#print Dumper @gene_ids;

#2.) Test convertGenomicCoordinates()
#+ve strand
#my $start = 99645975;
#my $end = 99661022;
#my $gene_id = 2;

#-ve strand
#my $start = 99689940;
#my $end = 99697939;
#my $gene_id = 1;

#my %coords = %{&convertGenomicCoordinates ('-dbh'=>$alexa_dbh, '-gene_id'=>$gene_id, '-start_pos'=>$start, '-end_pos'=>$end)};
#print Dumper %coords;

#3.) Test getExons();
#my @gene_ids = ('1','2','3');

#print Dumper @gene_ids;

#my %exons = %{getExons ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids, '-sequence'=>"no")};
#my %exons = %{getExons ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids, '-sequence'=>"yes")};

#print Dumper %exons;


#4.) Test getTranscripts()

#my @gene_ids = ('1','2','3');
#my $gene_transcripts_ref = &getTranscripts ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids, '-sequence'=>"no");

#my $gene_records = keys %{$gene_transcripts_ref};
#print "\nFound records for $gene_records genes\n\n";

#print Dumper %{$gene_transcripts_ref};


#5.) Test getIntronContent()

#my @gene_ids = ('1','2','3');
#my $genes_ref = &getIntronContent ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids);

#print Dumper $genes_ref;


#6.) Test getGeneTerms
#my $type = "EntrezGene";
#my @gene_ids = ('1','2','3');
#my $gene_ids_ref = &getGeneTerms ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids, '-id_type'=>$type);
#print Dumper $gene_ids_ref;


#7.) Test getGeneProbes()
#my @gene_ids = ('1','5','8');

#unfiltered
#my $gene_probes_ref = &getGeneProbes ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids, '-filtered'=>"no", '-silent'=>"yes");
#print Dumper $gene_probes_ref;
#print Dumper $gene_probes_ref->{total_probe_count};

#filtered
#my $gene_probes_ref = &getGeneProbes ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids, '-filtered'=>"yes", '-filter_set'=>1, '-silent'=>"yes");
#print Dumper $gene_probes_ref;
#print Dumper $gene_probes_ref->{total_probe_count};


#8.) Test getGeneIds()
#my @ensembl_g_ids = qw (ENSG00000000000 ENSG00000001167 ENSG00000001626 ENSG00000001629 ENSG00000001631);
#my @ensembl_g_ids = qw (ENSG00000000000);
#my @ensembl_g_ids = qw ();
#my $gene_ids_ref = &getGeneIds ('-dbh'=>$alexa_dbh, '-ensembl_g_ids'=>\@ensembl_g_ids);

#my @alexa_ids = qw (1 2 3);
#my @alexa_ids = qw (1);
#my @alexa_ids = qw (0);
#my $gene_ids_ref = &getGeneIds ('-dbh'=>$alexa_dbh, '-alexa_ids'=>\@alexa_ids);

#print Dumper $gene_ids_ref;


#Dump database name
my $test = $alexa_dbh->{Name};
print "\nDEBUG\n";
print Dumper $test;

#Close database connection
$alexa_dbh->disconnect();

exit();
