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
#my $database = "ALEXA_hs_51_36m";
#my $database = "ALEXA_hs_49_36k";
my $database = "ALEXA_hs_53_36o";
my $user = "viewer";
my $pass = "viewer";

#e.g. cat /projects/malachig/solexa/figures_and_stats/DE_SI_intersect/5FU/all_genes.txt  | ensg_to_entrez.pl

#Open database connection
my $alexa_dbh = &connectDB('-database'=>$database, '-server'=>$server, '-user'=>$user, '-password'=>$pass);


#Define type of IDs you want to convert to...
my $type = "EntrezGene";
#my $type = "HGNC";

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

while(<STDIN>){
#    print $_;
    chomp ($_);
    my $ensg;
    if ($_ =~ /(ENSG\S+)/){
      $ensg = $1;
    }elsif($_ =~ /(ENSG\d+)/){
      $ensg = $1;
    }else{
      print "$_\n", RESET;
      next();
    }
    my $gene_id = $ensg_to_alexa{$ensg}{gene_id};
    my @converted_ids = keys %{$gene_terms_ref->{$gene_id}->{external_ids}};

    foreach my $id (@converted_ids){
      print "$id\n";
    }
    #print "$_\t@converted_ids\n";

}

exit();


