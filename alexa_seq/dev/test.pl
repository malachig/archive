#!/usr/bin/perl -w
#Written by Malachi Griffith

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;

use lib '/home/malachig/svn/solexa_analysis';
use utilities::utility qw(:all);
use utilities::ALEXA_DB qw(:all);

my $database = '';
my $server = '';
my $user = '';
my $password = '';
my $chr_filter = '';

GetOptions ('database=s'=>\$database,'server=s'=>\$server, 'user=s'=>\$user, 'password=s'=>\$password, 'chr_filter=s'=>\$chr_filter);

print GREEN, "\n\nExample: test.pl  --database=ALEXA_hs_53_36o  --server=jango.bcgsc.ca  --user=viewer  --password=viewer  --chr_filter='3:16:121020102-126260001'\n\n", RESET;

unless ($database && $server && $user && $password && $chr_filter){
  print RED, "\nRequired input parameter(s) missing!\n\n", RESET;
  exit();
}

my $region_number;
my $start_filter;
my $end_filter;
if ($chr_filter =~ /(.*)\:(\d+)\:(\d+)\-(\d+)/){
  $chr_filter = $1;
  $region_number = $2;
  $start_filter = $3;
  $end_filter = $4;
  unless ($end_filter > $start_filter){
    print RED, "\nStart of range must be smaller than end ($chr_filter)\n\n", RESET;
    exit();
  }
}else{
  print RED, "\nFormat of chr_filter not understood: $chr_filter (should be of the form:  Y:1:1-9999001)\n\n", RESET;
  exit();
}

my $alexa_dbh = &connectDB('-database'=>$database, '-server'=>$server, '-user'=>$user, '-password'=>$password);
my @gene_ids = @{&getAllGenes ('-dbh'=>$alexa_dbh, '-gene_type'=>'All', '-chromosome'=>$chr_filter, '-range'=>"$start_filter-$end_filter")};
my $protein_bases_ref = &getProteinBases ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids);
$alexa_dbh->disconnect();

#Test on a single set of coordinates

my $chromosome = 3;
my $chr_start = 125931903;
my $chr_end = 125946730;
my $protein_bases_count = &getProteinBaseCount ('-protein_bases_ref'=>$protein_bases_ref, '-chromosome'=>$chromosome, '-chr_start'=>$chr_start, '-chr_end'=>$chr_end);

print Dumper $protein_bases_count;


exit();




