#!/usr/bin/perl -w

#Written by Malachi Griffith
#This script takes as input a tab delimeted file with EnsEMBL gene IDs in the first column and finds orthologs
#The user must also provide the species to get the orthologs from (species must be indicated by an NCBI taxon ID)
#Orthologs will be extracted from the EnsEMBL compara database using the EnsEMBL api

#This script relies on working installations of bioperl as well as the ensembl and ensembl-compara APIs
#Once you have installed those, update the code below to account for new versions of the API
#Instructions for installing EnsEMBL APIs can be found here:
#   http://www.ensembl.org/info/software/core/index.html
#   http://www.ensembl.org/info/software/compara/index.html

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);


#Initialize command line options
my $ensembl_api_version = ''; #Version of EnsEMBL to use    
my $gene_file = '';           #Tab delimited file containing ensembl gene IDs in the first column (header line expected)
my $target_taxon_id = '';     #Valid ncbi taxon id for species to get orthologs for (see http://www.ncbi.nlm.nih.gov/Taxonomy/)
my $min_percent_id = '';      #Minimum sequence identity between orthologs
my $outfile = '';             #Output file, gene_ids for orthologs will be appended as the FIRST column of this file
my $logfile = '';

GetOptions ('ensembl_api_version=s'=>\$ensembl_api_version, 'gene_file=s'=>\$gene_file, 
	    'target_taxon_id=i'=>\$target_taxon_id, 'min_percent_id=f'=>\$min_percent_id, 
            'outfile=s'=>\$outfile, 'logfile=s'=>\$logfile);


#Provide instruction to the user
print GREEN, "\n\nNOTE: You may need to adjust this script to include the desired API version\n\n", RESET;
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the correct EnsEMBL API version using: --ensembl_api_version (41, 42, etc.)", RESET;
print GREEN, "\n\tSpecify a input gene file containing EnsEMBL ids in the first column using: --gene_file", RESET;
print GREEN, "\n\tSpecify the species you wish to get orthologs from using (by taxon id): --target_taxon_id", RESET;
print GREEN, "\n\t\tMust be a valid NCBI taxon ID:", RESET;
print GREEN, "\n\t\t\tHuman=9606, Xenopus=8364, Rat=10116, Yeast=4932, Mouse=10090\n\t\t\tChimp=9598, Chicken=9031, Fly=7227, Worm=6239, Cow=9913\n\t\t\tDog=9615, Zebrafish=7955, Frog=8364\n", RESET;
print GREEN, "\n\tSpecify the minimum percent identity between orthologous sequences retrieved using: --min_percent_id", RESET;
print GREEN, "\n\tSpecify the output file using: --outfile", RESET;
print GREEN, "\n\tSpecify a logfile for output using: --logfile", RESET;
print GREEN, "\n\nExample: getEnsemblOrthologs.pl  --ensembl_api_version=41  --gene_file=/home/user/alexa/human_housekeeping_genes.txt  --target_taxon_id=10090  --min_percent_id=90.0  --outfile=/home/user/alexa/ALEXA_version/mouse_housekeeping_genes.txt  --logfile=/home/user/alexa/ALEXA_version/logs/getEnsemblOrthologs_LOG.txt\n\n", RESET;

unless ($ensembl_api_version && $gene_file && $target_taxon_id && $min_percent_id && $outfile && $logfile){
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  exit();
}

#**********************************************************************************************************
#IMPORTANT NOTE: You must have the correct Ensembl API installed locally AND bioperl 1.2 or greater!!
#Both the EnsEMBL core API as well as Compara are required
#Refer to the ALEXA manual for additional details on how to install these
#Then update the following paths:
if ($ensembl_api_version =~ /^\d+/){
  if ($ensembl_api_version eq "27"){
    unshift(@INC, "/home/alexa/perl/ensembl_27_perl_API/ensembl/modules");
    unshift(@INC, "/home/alexa/perl/ensembl_27_perl_API/ensembl-compara/modules");
  }elsif($ensembl_api_version eq "28"){
    unshift(@INC, "/home/alexa/perl/ensembl_28_perl_API/ensembl/modules");
    unshift(@INC, "/home/alexa/perl/ensembl_28_perl_API/ensembl-compara/modules");
  }elsif($ensembl_api_version eq "29"){
    unshift(@INC, "/home/alexa/perl/ensembl_29_perl_API/ensembl/modules");
    unshift(@INC, "/home/alexa/perl/ensembl_29_perl_API/ensembl-compara/modules");
  }elsif($ensembl_api_version eq "30"){
    unshift(@INC, "/home/alexa/perl/ensembl_30_perl_API/ensembl/modules");
    unshift(@INC, "/home/alexa/perl/ensembl_30_perl_API/ensembl-compara/modules");
  }elsif($ensembl_api_version eq "31"){
    unshift(@INC, "/home/alexa/perl/ensembl_31_perl_API/ensembl/modules");
    unshift(@INC, "/home/alexa/perl/ensembl_31_perl_API/ensembl-compara/modules");
  }elsif($ensembl_api_version eq "32"){
    unshift(@INC, "/home/alexa/perl/ensembl_32_perl_API/ensembl/modules");
    unshift(@INC, "/home/alexa/perl/ensembl_32_perl_API/ensembl-compara/modules");
  }elsif($ensembl_api_version eq "33"){
    unshift(@INC, "/home/alexa/perl/ensembl_33_perl_API/ensembl/modules");
    unshift(@INC, "/home/alexa/perl/ensembl_33_perl_API/ensembl-compara/modules");
  }elsif($ensembl_api_version eq "34"){
    unshift(@INC, "/home/alexa/perl/ensembl_34_perl_API/ensembl/modules");
    unshift(@INC, "/home/alexa/perl/ensembl_34_perl_API/ensembl-compara/modules");
  }elsif($ensembl_api_version eq "35"){
    unshift(@INC, "/home/alexa/perl/ensembl_35_perl_API/ensembl/modules");
    unshift(@INC, "/home/alexa/perl/ensembl_35_perl_API/ensembl-compara/modules");
  }elsif($ensembl_api_version eq "36"){
    unshift(@INC, "/home/alexa/perl/ensembl_36_perl_API/ensembl/modules");
    unshift(@INC, "/home/alexa/perl/ensembl_36_perl_API/ensembl-compara/modules");
  }elsif($ensembl_api_version eq "37"){
    unshift(@INC, "/home/alexa/perl/ensembl_37_perl_API/ensembl/modules");
    unshift(@INC, "/home/alexa/perl/ensembl_37_perl_API/ensembl-compara/modules");
  }elsif($ensembl_api_version eq "38"){
    unshift(@INC, "/home/alexa/perl/ensembl_38_perl_API/ensembl/modules");
    unshift(@INC, "/home/alexa/perl/ensembl_38_perl_API/ensembl-compara/modules");
  }elsif($ensembl_api_version eq "39"){
    unshift(@INC, "/home/alexa/perl/ensembl_39_perl_API/ensembl/modules");
    unshift(@INC, "/home/alexa/perl/ensembl_39_perl_API/ensembl-compara/modules");
  }elsif($ensembl_api_version eq "40"){
    unshift(@INC, "/home/alexa/perl/ensembl_40_perl_API/ensembl/modules");
    unshift(@INC, "/home/alexa/perl/ensembl_40_perl_API/ensembl-compara/modules");
  }elsif($ensembl_api_version eq "41"){
    unshift(@INC, "/home/alexa/perl/ensembl_41_perl_API/ensembl/modules");
    unshift(@INC, "/home/alexa/perl/ensembl_41_perl_API/ensembl-compara/modules");
  }elsif($ensembl_api_version eq "42"){
    unshift(@INC, "/home/alexa/perl/ensembl_42_perl_API/ensembl/modules");
    unshift(@INC, "/home/alexa/perl/ensembl_42_perl_API/ensembl-compara/modules");
  }elsif($ensembl_api_version eq "43"){
    unshift(@INC, "/home/alexa/perl/ensembl_43_perl_API/ensembl/modules");
    unshift(@INC, "/home/alexa/perl/ensembl_43_perl_API/ensembl-compara/modules");
  }elsif($ensembl_api_version eq "44"){
    unshift(@INC, "/home/alexa/perl/ensembl_44_perl_API/ensembl/modules");
    unshift(@INC, "/home/alexa/perl/ensembl_44_perl_API/ensembl-compara/modules");
  }elsif($ensembl_api_version eq "45"){
    unshift(@INC, "/home/alexa/perl/ensembl_45_perl_API/ensembl/modules");
    unshift(@INC, "/home/alexa/perl/ensembl_45_perl_API/ensembl-compara/modules");
  }elsif($ensembl_api_version eq "46"){
    unshift(@INC, "/home/alexa/perl/ensembl_46_perl_API/ensembl/modules");
    unshift(@INC, "/home/alexa/perl/ensembl_46_perl_API/ensembl-compara/modules");
  }else{
    print RED, "\nEnsEMBL API version: $ensembl_api_version is not defined, modify script before proceeding\n\n", RESET;
    exit();
  }
}else{
  print RED, "\nEnsEMBL API version format: $ensembl_api_version not understood!\n\n", RESET;
  exit();
}
use lib "/home/alexa/perl/bioperl-1.4";    #Bioperl
#*********************************************************************************************************
#print Dumper @INC;
require Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
require Bio::EnsEMBL::Registry;

open (LOG, ">$logfile") || die "\nCould not open logfile: $logfile\n\n";

print LOG "\nUser Specified the following options:\ngene_file = $gene_file\ntarget_taxon_id = $target_taxon_id\nmin_percent_id = $min_percent_id\noutfile = $outfile\nlogfile = $logfile\n\n";

#Parse input file
my $header;
my %genes = &parseGene_file('-infile'=>$gene_file);
my $new_species;

#Create a connection to the ensembl database online - Note that the following queries may be quite slow ...
#By default, only databases corresponding to the API version specified above are loaded.
Bio::EnsEMBL::Registry->load_registry_from_db(-host => 'ensembldb.ensembl.org', -user => 'anonymous', -port => 3306);

my $total_homologs_found = 0;
my $genes_with_homologs = 0;
my $genes_without_homologs = 0;
my $gene_count = 0;

foreach my $gene_id (keys %genes){
  my $homologs_found = 0;
  $gene_count++;

  #First get a Member object. In the case of homology this is a gene, in the case of family it can be a gene or a protein
  my $member_adaptor = Bio::EnsEMBL::Registry->get_adaptor('Multi', 'compara', 'Member');
  my $member = $member_adaptor->fetch_by_source_stable_id('ENSEMBLGENE', $gene_id);

  unless($member){
    print YELLOW, "\nCould not create Member object for EnsEMBL ID: $gene_id!\n\n", RESET;
    print LOG "\nCould not create Member object for EnsEMBL ID: $gene_id!\n\n";
    next();
  }

  my $taxon = $member->taxon;
  my $common_name = $taxon->common_name;
  my $binomial_species = $taxon->binomial;

  print BLUE, "\n\n$gene_count\tProcessing: $gene_id.  This gene comes from species: $common_name ($binomial_species)", RESET;
  print LOG "\n\n$gene_count\tProcessing: $gene_id.  This gene comes from species: $common_name ($binomial_species)";

  #Now get the homologies where the member is involved
  my $homology_adaptor = Bio::EnsEMBL::Registry->get_adaptor('Multi', 'compara', 'Homology');
  my $homologies = $homology_adaptor->fetch_by_Member($member);

  foreach my $homology (@{$homologies}){

    my $description = $homology->description;
    my $subtype = $homology->subtype;


    #Each homology relation has 2 or more members
    foreach my $member_attribute (@{$homology->get_all_Member_Attribute}){

      my ($member, $attribute) = @{$member_attribute};

      my $stable_id = $member->stable_id;
      my $taxon_id = $member->taxon_id;
      my $percent_identity = $attribute->perc_id;

      #my $source_name = $member->source_name;

      #Only consider those homologies involving the target species specified by the user - Also limit to greater than min_percent_identity
      if ($taxon_id == $target_taxon_id && $percent_identity >= $min_percent_id){
	$homologs_found++;
	$total_homologs_found++;

	#Make note of the common name associated with this homologous sequence
	my $taxon = $member->taxon;
	my $common_name = $taxon->common_name;
        my $binomial_species = $taxon->binomial;

	unless ($common_name){
	  $common_name = "na";
	}

	$new_species = $common_name;

        my $new_gene_name = '';

        my $slice_adaptor = Bio::EnsEMBL::Registry->get_adaptor($binomial_species, 'core', 'Slice');

        my $slice = $slice_adaptor->fetch_by_gene_stable_id($stable_id, 0);

        foreach my $gene (@{$slice->get_all_Genes()}){
          my $gene_id = $gene->stable_id();
          unless ($gene_id eq $stable_id){
            next();
          }

          my $test_name = $gene->external_name();
          if ($test_name){
            $new_gene_name = $test_name;
          }
          unless ($test_name){
            $new_gene_name = "Unknown";
          }
        }

	print BLUE, "\n\t\tHomology description: $description\tSubtype: $subtype", RESET;
	print BLUE, "\n\t\t\tID: $stable_id\tName: $new_gene_name\tSpecies: $common_name\tPercent_identity: $percent_identity", RESET;
	print LOG "\n\t\tHomology description: $description\tSubtype: $subtype";
	print LOG "\n\t\t\tID: $stable_id\tName: $new_gene_name\tSpecies: $common_name\tTaxon: $taxon_id\tPercent_identity: $percent_identity";

	my $homologs_ref = $genes{$gene_id}{homologs};
	$homologs_ref->{$stable_id}->{percent_id} = $percent_identity;
	$homologs_ref->{$stable_id}->{homolog_type} = $description;
        $homologs_ref->{$stable_id}->{gene_name} = $new_gene_name;
      }
    }
  }
  if ($homologs_found > 0){
    $genes_with_homologs++;
    print BLUE, "\n\t\tFound $homologs_found homologs for this gene", RESET;
    print LOG "\n\t\tFound $homologs_found homologs for this gene";
  }else{
    $genes_without_homologs++;
    print YELLOW, "\n\t\tCould not find any homologs for this gene", RESET;
    print LOG "\n\t\tCould not find any homologs for this gene";
  }
}

print BLUE, "\nFound a total of $total_homologs_found homologs for the $gene_count genes supplied\n\n", RESET;
print LOG "\nFound a total of $total_homologs_found homologs for the $gene_count genes supplied\n\n";

print BLUE, "\nFound $genes_with_homologs genes with at least one homolog and $genes_without_homologs genes without any homologs\n\n", RESET;
print LOG "\nFound $genes_with_homologs genes with at least one homolog and $genes_without_homologs genes without any homologs\n\n";

#Finally print out an appended file with the orthologs found
print BLUE, "\nPrinting appended output file to: $outfile\n\n", RESET;
print LOG "\nPrinting appended output file to: $outfile\n\n";

open (OUT, ">$outfile") || die "\nCould not open output file: $outfile\n\n";

$new_species =~ tr/ /\_/;

my $new_column = "ensembl_g_id"."_"."$new_species";
print OUT "$new_column\tpercent_identity\thomolog_type\t$header\n";

foreach my $gene_id (sort {$a cmp $b} keys %genes){

  my $homologs_ref = $genes{$gene_id}{homologs};

  foreach my $homolog (sort keys %{$homologs_ref}){
    print OUT "$homolog\t$homologs_ref->{$homolog}->{gene_name}\t$homologs_ref->{$homolog}->{percent_id}\t$homologs_ref->{$homolog}->{homolog_type}\t$genes{$gene_id}{line_record}\n";
  }
}

close (OUT);
close (LOG);

exit();


###################################################################################################
#Parse input file                                                                                 #
###################################################################################################
sub parseGene_file{
  my %args = @_;
  my $infile = $args{'-infile'};

  open (GENES, "$infile") || die "\nCould not open input file: $infile\n\n";

  print BLUE, "\nParsing the first column of $infile for EnsEMBL gene Ids\n\n", RESET;
  print LOG "\nParsing the first column of $infile for EnsEMBL gene Ids\n\n";

  my $first_line = 1;
  my %genes;

  while (<GENES>){
    chomp($_);

    if ($first_line == 1){
      $header = $_;
      $first_line = 0;
      next();
    }

    my @line = split ("\t", $_);

    $genes{$line[0]}{line_record} = $_;
    my %homologs;
    $genes{$line[0]}{homologs} = \%homologs;
  }

  close (GENES);

  my $gene_count = keys %genes;
  print BLUE, "\nFound $gene_count EnsEMBL gene Ids\n\n", RESET;
  print LOG "\nFound $gene_count EnsEMBL gene Ids\n\n";

  return(%genes);
}
