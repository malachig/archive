#!/usr/bin/perl -w
#Written by Malachi Griffith
#The purpose of this script is to get gene descriptions for all Genes in EnsEMBL

use DBI;
use strict;
use Data::Dumper;
use Getopt::Long;
use Benchmark;
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
my $ensembl_api_version = ''; #Version of EnsEMBL to use
my $species = '';
my $ensembl_database = '';
my $ensembl_server = '';
my $ensembl_user = '';
my $ensembl_password = '';
my $outfile = '';

GetOptions ('ensembl_api_version=s'=>\$ensembl_api_version, 'species=s'=>\$species,
	    #Source database access
	    'ensembl_database=s'=>\$ensembl_database, 'ensembl_server=s'=>\$ensembl_server, 'ensembl_user=s'=>\$ensembl_user, 'ensembl_password=s'=>\$ensembl_password,
	    'outfile=s'=>\$outfile);

#Provide instruction to the user
print GREEN, "\n\nNOTE: Before using this script, make sure the correct API version is hard coded!!\n\n", RESET;
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the correct EnsEMBL API version using: --ensembl_api_version (41, 42, etc.)", RESET;
print GREEN, "\n\tSpecify the SOURCE Ensembl Database, Server, User and Password using: --ensembl_database  --ensembl_server  --ensembl_user  and  --ensembl_password", RESET;
print GREEN, "\n\tThe species using: --species (e.g. --species=Human or --species='Homo sapiens')", RESET;
print GREEN, "\n\t\tMake sure the species you supply matches the EnsEMBL database you supply!!", RESET;
print GREEN, "\n\tSpecify the name of an output file using: --outfile", RESET;
print GREEN, "\n\nExample: getEnsemblGeneDescriptions.pl  --ensembl_api_version=48  --species=Human  --ensembl_database=homo_sapiens_core_48_36j  --ensembl_server=ensembl01.bcgsc.ca  --ensembl_user=ensembl  --ensembl_password=ensembl  --outfile=EnsEMBL_48_GeneDescriptions.txt\n\n", RESET;

unless ($ensembl_api_version && $species && $ensembl_database && $ensembl_server && $ensembl_user && $ensembl_password && $outfile){
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  exit();
}

#**********************************************************************************************************
#IMPORTANT NOTE: You must have the correct Ensembl API installed locally AND bioperl 1.2 or greater!!
#Both the EnsEMBL core API as well as Compara are required
#Refer to the ALEXA manual for additional details on how to install these
#Then update the following paths:
if ($ensembl_api_version <= 33){
  print RED, "\nEnsEMBL API version earlier than v34 do not work with the connection types used in this script!\n\n", RESET;
  exit();
}
if ($ensembl_api_version =~ /^\d+/){
  if ($ensembl_api_version eq "34"){
    unshift(@INC, "/home/malachig/perl/ensembl_34_perl_API/ensembl/modules");
    unshift(@INC, "/home/malachig/perl/ensembl_34_perl_API/ensembl-variation/modules");
  }elsif ($ensembl_api_version eq "35"){
    unshift(@INC, "/home/malachig/perl/ensembl_35_perl_API/ensembl/modules");
    unshift(@INC, "/home/malachig/perl/ensembl_35_perl_API/ensembl-variation/modules");
  }elsif ($ensembl_api_version eq "36"){
    unshift(@INC, "/home/malachig/perl/ensembl_36_perl_API/ensembl/modules");
    unshift(@INC, "/home/malachig/perl/ensembl_36_perl_API/ensembl-variation/modules");
  }elsif ($ensembl_api_version eq "37"){
    unshift(@INC, "/home/malachig/perl/ensembl_37_perl_API/ensembl/modules");
    unshift(@INC, "/home/malachig/perl/ensembl_37_perl_API/ensembl-variation/modules");
  }elsif ($ensembl_api_version eq "38"){
    unshift(@INC, "/home/malachig/perl/ensembl_38_perl_API/ensembl/modules");
    unshift(@INC, "/home/malachig/perl/ensembl_38_perl_API/ensembl-variation/modules");
  }elsif ($ensembl_api_version eq "39"){
    unshift(@INC, "/home/malachig/perl/ensembl_39_perl_API/ensembl/modules");
    unshift(@INC, "/home/malachig/perl/ensembl_39_perl_API/ensembl-variation/modules");
  }elsif ($ensembl_api_version eq "40"){
    unshift(@INC, "/home/malachig/perl/ensembl_40_perl_API/ensembl/modules");
    unshift(@INC, "/home/malachig/perl/ensembl_40_perl_API/ensembl-variation/modules");
  }elsif ($ensembl_api_version eq "41"){
    unshift(@INC, "/home/malachig/perl/ensembl_41_perl_API/ensembl/modules");
    unshift(@INC, "/home/malachig/perl/ensembl_41_perl_API/ensembl-variation/modules");
  }elsif($ensembl_api_version eq "42"){
    unshift(@INC, "/home/malachig/perl/ensembl_42_perl_API/ensembl/modules");
    unshift(@INC, "/home/malachig/perl/ensembl_42_perl_API/ensembl-variation/modules");
  }elsif($ensembl_api_version eq "43"){
    unshift(@INC, "/home/malachig/perl/ensembl_43_perl_API/ensembl/modules");
    unshift(@INC, "/home/malachig/perl/ensembl_43_perl_API/ensembl-variation/modules");
  }elsif($ensembl_api_version eq "44"){
    unshift(@INC, "/home/malachig/perl/ensembl_44_perl_API/ensembl/modules");
    unshift(@INC, "/home/malachig/perl/ensembl_44_perl_API/ensembl-variation/modules");
  }elsif($ensembl_api_version eq "45"){
    unshift(@INC, "/home/malachig/perl/ensembl_45_perl_API/ensembl/modules");
    unshift(@INC, "/home/malachig/perl/ensembl_45_perl_API/ensembl-variation/modules");
  }elsif($ensembl_api_version eq "46"){
    unshift(@INC, "/home/malachig/perl/ensembl_46_perl_API/ensembl/modules");
    unshift(@INC, "/home/malachig/perl/ensembl_46_perl_API/ensembl-variation/modules");
  }elsif($ensembl_api_version eq "47"){
    unshift(@INC, "/home/malachig/perl/ensembl_47_perl_API/ensembl/modules");
    unshift(@INC, "/home/malachig/perl/ensembl_47_perl_API/ensembl-variation/modules");
  }elsif($ensembl_api_version eq "48"){
    unshift(@INC, "/home/malachig/perl/ensembl_48_perl_API/ensembl/modules");
    unshift(@INC, "/home/malachig/perl/ensembl_48_perl_API/ensembl-variation/modules");
  }elsif($ensembl_api_version eq "49"){
    unshift(@INC, "/home/malachig/perl/ensembl_49_perl_API/ensembl/modules");
    unshift(@INC, "/home/malachig/perl/ensembl_49_perl_API/ensembl-variation/modules");
  }elsif($ensembl_api_version eq "50"){
    unshift(@INC, "/home/malachig/perl/ensembl_50_perl_API/ensembl/modules");
    unshift(@INC, "/home/malachig/perl/ensembl_50_perl_API/ensembl-variation/modules");
  }elsif($ensembl_api_version eq "51"){
    unshift(@INC, "/home/malachig/perl/ensembl_51_perl_API/ensembl/modules");
    unshift(@INC, "/home/malachig/perl/ensembl_51_perl_API/ensembl-variation/modules");
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
require Bio::EnsEMBL::DBSQL::DBAdaptor; #Used for local connections to EnsEMBL core databases
require Bio::EnsEMBL::Variation::DBSQL::DBAdaptor; #Used for local connections to EnsEMBL variation databases


#1.) Establish connections to source EnsEMBL database and target ALEXA databases

#A.) Using a local ensembl version

##CONNECT TO ENSEMBL SERVER
my $registry = "Bio::EnsEMBL::Registry";
$registry->load_registry_from_db(-host =>$ensembl_server, -user =>$ensembl_user, -pass=>$ensembl_password, -db_version =>$ensembl_api_version);

#Get a connection to the local Ensembl CORE database
my $ensembl_core_api = Bio::EnsEMBL::Registry->get_DBAdaptor($species, "core");

#Get a connection to the local Ensembl VARIATION database
my $ensembl_variation_api = Bio::EnsEMBL::Registry->get_DBAdaptor($species, "variation");
my $vf_adaptor = $ensembl_variation_api->get_VariationFeatureAdaptor();

my $ensembl_ids_ref = &getAllEnsemblIds();
my $source = $ensembl_database;

my $gene_count = 0;
my $known_gene_count = 0;
my $unknown_gene_count = 0;
my $gene_update_count = 0;
my $chr_strand;
my $total_bases = 0;
my $masked_bases = 0;

open (OUT, ">$outfile") || die "\nCould not open output file: $outfile\n\n";
print OUT "EnsEMBL_ID\tGene_Description\n";

foreach my $ensembl_id (sort keys %{$ensembl_ids_ref}){
  my $t0 = new Benchmark;
  my %gene_info;
  my $gene_info_ref = \%gene_info;
  $gene_count++;

  #A Slice object represents a continuous region of a genome
  my $slice_adaptor = $registry->get_adaptor($species, 'Core', 'Slice');
  #my $slice_adaptor = $ensembl_core_api->get_SliceAdaptor(); #Alternate way of getting a 'slice adaptor'

  #Retrieve a slice with respect to a gene, with a specified flanking sequence on either side
  my $slice = $slice_adaptor->fetch_by_gene_stable_id($ensembl_id, 0);

  #BEFORE DOING ANYTHING: Check the strand of the target gene within this strand.  If it is on the reverse strand - invert the slice and continue
  #To simplify matters we always want the slice to be oriented so that the first exon (in terms of transcription, coding) at the beginning of the slice
  $chr_strand = 1; #Initialize to 1, switch to -1 if gene is on reverse strand.

  #Determine the coordinate system used here - should be chromosome for most genomes
  #These coordinates can be used to map the location of genes, transcripts, and exons to their chromosomal positions
  my $coord_sys = $slice->coord_system()->name();
  my $chr = $slice->seq_region_name();

  my $chr_start = $slice->start();
  my $chr_end = $slice->end();

  #Remember this slice may contain more than one gene if there are overlapping genes
  #The boundaries of the slice will only correspond exactly to the gene ID we actually used to get the slice in the first place
  # - for this gene, the slice will be from the first base of the first exon to the last base of the last exon.
  foreach my $gene (@{$slice->get_all_Genes()}){

    #Make sure the right gene is being considered and not an overlapping gene (which will be processed individually)
    my $gene_id = $gene->stable_id();
    unless ($gene_id eq $ensembl_id){
      next();
    }

    my $test_name = $gene->external_name();
    if ($test_name){
      $known_gene_count++;
    }
    unless ($test_name){
      $test_name = "Unknown";
      $unknown_gene_count++;
    }
    my $gene_description = $gene->description();
    unless ($gene_description){
      $gene_description = "NA";
    }

    print BLUE, "\n$gene_count -> Gene NAME: $test_name ($ensembl_id)", RESET;
    print OUT "$ensembl_id\t$gene_description\n";

  }

}

print BLUE, "\n\nProcessed a total of $gene_count genes: $known_gene_count are KNOWN, $unknown_gene_count are UNKNOWN (predicted only)\n", RESET;

print "\n\n";

close(OUT);

exit();


##############################################################################################
#Get all ensembl Ids from the local ensembl database that are appropriate                    #
##############################################################################################
sub getAllEnsemblIds{
  my %ids;

  #Get a normal DBI connection to the ensembl database
  my $ensembl_dbh = &connectDB('-database'=>$ensembl_database, '-server'=>$ensembl_server, '-user'=>$ensembl_user, '-password'=>$ensembl_password);

  #Get the gene sequence for this gene
  my $sql = "SELECT gene.gene_id,stable_id,biotype FROM gene,gene_stable_id WHERE gene_stable_id.gene_id = gene.gene_id;";
  my $sth = $ensembl_dbh->prepare("$sql");
  $sth->execute();

  while (my ($gene_id,$ensembl_id,$gene_type) = $sth->fetchrow_array()){
    $ids{$ensembl_id}{ensembl_internal_id} = $gene_id;
    $ids{$ensembl_id}{gene_type} = $gene_type;
  }

  $sth->finish();

  $ensembl_dbh->disconnect();

  return(\%ids);
}



