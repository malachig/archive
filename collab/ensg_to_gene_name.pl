#!/usr/bin/env perl

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;

use lib '/home/malachig/svn/solexa_analysis';
use utilities::utility qw(:all);
use utilities::ALEXA_DB qw(:all);

my $server = "jango.bcgsc.ca";
my $database = "ALEXA_hs_51_36m";
#my $database = "ALEXA_hs_49_36k";
my $user = "viewer";
my $pass = "viewer";
my $name_type = '';

GetOptions ('name_type=s'=>\$name_type);

unless ($name_type =~ /hugo|entrez/i){
  print RED, "\n\nSpecify the type of names you would like using: --name_type (hugo or entrez)\n\ne.g. cat genelist.txt | ensg_to_gene_name.pl --name_type=hugo", RESET;
  exit();
}

#Open database connection
my $alexa_dbh = &connectDB('-database'=>$database, '-server'=>$server, '-user'=>$user, '-password'=>$pass);

#First get EntrezGene names from an ALEXA datbase using the getGeneTerms
my $type = "EntrezGene";
my $gene_storable = "$database"."_AllGenes_GeneInfo_NoSeq.storable";
my $gene_ids_ref = &getAllGenes ('-dbh'=>$alexa_dbh, '-gene_type'=>'All');
my $gene_terms_ref = &getGeneTerms ('-dbh'=>$alexa_dbh, '-gene_ids'=>$gene_ids_ref, '-id_type'=>$type, '-silent'=>"yes");
my $genes_ref = &getGeneInfo ('-dbh'=>$alexa_dbh, '-gene_ids'=>$gene_ids_ref, '-sequence'=>"no", '-storable'=>$gene_storable, '-silent'=>"yes");

my %ensg_to_alexa;
foreach my $gene_id (keys %{$genes_ref}){
  my $ensg = $genes_ref->{$gene_id}->{ensembl_g_id};
  $ensg_to_alexa{$ensg}{gene_id} = $gene_id;
}

#Close database connection
$alexa_dbh->disconnect();

#Get Entrez-to-HGNC mapping from Entrez
#ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz
my $map_file = "/home/malachig/svn/collab/Homo_sapiens.gene_info";
my %entrez_to_hgnc;

open(MAP, "$map_file") || die "\nCould not open Entrez info file: $map_file\n\n";
while(<MAP>){
  chomp($_);
  my @line = split("\t", $_);
  my $entrez_id = $line[2];
  my $hgnc_id = $line[10];
  if ($hgnc_id eq "-"){
    next();
  }
  if ($entrez_to_hgnc{$entrez_id}{hgnc}){
    push(@{$entrez_to_hgnc{$entrez_id}{hgnc}}, $hgnc_id);
  }else{
    my @tmp;
    push (@tmp, $hgnc_id);
    $entrez_to_hgnc{$entrez_id}{hgnc} = \@tmp;
  }
}
close(MAP);

while(<STDIN>){
  chomp ($_);
  my $ensg;
  if ($_ =~ /(ENSG\S+)/){
    $ensg = $1;
  }elsif($_ =~ /(ENSG\d+)/){
    $ensg = $1;
  }else{
    print "$_\n";
    next();
  }
  my $gene_id = $ensg_to_alexa{$ensg}{gene_id};
  my @entrez_ids = keys %{$gene_terms_ref->{$gene_id}->{external_ids}};

  my %hgnc_ids;
  foreach my $entrez_id (@entrez_ids){
    if ($entrez_to_hgnc{$entrez_id}{hgnc}){
      my @mapped_ids = @{$entrez_to_hgnc{$entrez_id}{hgnc}};
      foreach my $hgnc_id (@mapped_ids){
        $hgnc_ids{$hgnc_id}{tmp} = '';
      }
    }
  }

  my $sep = $";
  my @hgnc_ids = keys %hgnc_ids;
  if ($name_type =~ /hugo/i){
    $" = ","; print "$_\t@hgnc_ids\n"; $" = $sep;
  }elsif ($name_type =~ /entrez/i){
    $" = ","; print "$_\t@entrez_ids\n"; $" = $sep;
  }else{
    print RED, "\nName type not understood\n\n", RESET;
    exit();
  }
}

exit();


