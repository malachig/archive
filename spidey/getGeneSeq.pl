#!/usr/bin/perl

use strict;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use Data::Dumper;

use lib '/home/malachig/svn/alexa_seq/';
use utilities::ALEXA_DB qw(:all);
use utilities::utility qw(:all);

my $target_ensg_id = '';
my $alexa_db = '';
my $working_dir = '';

GetOptions( 'target_ensg_id=s' => \$target_ensg_id,
            'alexa_db=s' => \$alexa_db,
            'working_dir=s' => \$working_dir);


unless ($target_ensg_id && $alexa_db && $working_dir)
{
  print MAGENTA, "Please specify the gene ID, using the following format:\n";
  print "--target_ensg_id=ENSG00000090061 --alexa_db=ALEXA_hs_53_36o --working_dir=/home/bcho/Temp/HS04391/CCNK/\n";
  print "getGeneSeq.pl terminated\n\n", RESET;
  exit();
}


#Get an EnsEMBL gene id from the user
#my $target_ensg_id = "ENSG00000114491"; #i.e. 'UMPS'
#my $target_ensg_id = "ENSG00000090061"; #i.e. CCNK

#Push it into an array
my @ensg_id_list;
push(@ensg_id_list, $target_ensg_id);

#Create a connection to the ALEXA database
my $dbh = &connectDB('-database'=>$alexa_db, '-server'=>'jango.bcgsc.ca', '-user'=>'viewer', '-password'=>'viewer');

#Get a gene info object for the specified gene ID:
#First get the 'ALEXA Gene Id' that corresponds to the ensembl gene id
my %gene_ids = %{&getGeneIds ('-dbh'=>$dbh, '-ensembl_g_ids'=>\@ensg_id_list)};
my $gene_id = $gene_ids{$target_ensg_id}{alexa_gene_id};

#Now get the gene info object
my %gene_info = %{&getGeneInfo ('-dbh'=>$dbh, '-gene_id'=>$gene_id, '-sequence'=>"yes")};

#Close the DB connection
$dbh->disconnect();

my $gene_seq = $gene_info{$gene_id}{sequence};
my $gene_name = $gene_info{$gene_id}{gene_name};
my $chr_start = $gene_info{$gene_id}{chr_start};
my $chr_end = $gene_info{$gene_id}{chr_end};
my $chr_strand = $gene_info{$gene_id}{chr_strand};
my $chr = $gene_info{$gene_id}{chromosome};

my $gene_seq_file = "$working_dir"."$gene_name"."_seq.fa";
open (GENE_SEQ, ">$gene_seq_file");
print GENE_SEQ "\>$gene_name\n";
print GENE_SEQ "$gene_seq";
close (GENE_SEQ);

#print the other gene information to standard out
#accept these as standard in from the identifyExpressedJunctions.pl script
my $gene_info_file = "$working_dir"."gene_info.txt";
open (GENE_INFO, ">$gene_info_file");
print GENE_INFO "$gene_name\n";
print GENE_INFO "$chr_start\n";
print GENE_INFO "$chr_end\n";
print GENE_INFO "$chr_strand\n";
print GENE_INFO "$chr\n";
close (GENE_INFO);


exit();



