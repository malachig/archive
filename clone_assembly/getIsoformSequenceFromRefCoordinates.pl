#!/usr/bin/perl -w

#Written by Malachi Griffith
#The purpose of this script is to build mRNA isoform sequences from the reference human genome using coordinates of individual exons
#These coordinates are taken from a BLAT psl file from clones being mapped to the genome
#The exon boundaries found in this alignment will be used to build a hypothetical isoform mRNA sequence from the reference human genome

use strict;
use Data::Dumper;
use Term::ANSIColor qw(:constants);

unshift(@INC, "/home/malachig/perl/ensembl_49_perl_API/ensembl/modules");
unshift(@INC, "/home/malachig/perl/ensembl_49_perl_API/ensembl-variation/modules");

use lib "/home/malachig/perl/bioperl-1.4"; 
require Bio::EnsEMBL::DBSQL::DBAdaptor; #Used for local connections to EnsEMBL core databases
require Bio::EnsEMBL::Variation::DBSQL::DBAdaptor; #Used for local connections to EnsEMBL variation databases

my $infile = "/home/malachig/AlternativeSplicing/UMPS_Clone_Assemblies/Isoform_Reference_Coordinates.txt";

my $clones_ref = &parseExonCoordFile('-infile'=>$infile);
#print Dumper %clones;

##CONNECT TO ENSEMBL SERVER
my $registry = "Bio::EnsEMBL::Registry";
$registry->load_registry_from_db(-host =>'ensembl01.bcgsc.ca', -user =>'ensembl', -pass=>'ensembl', -db_version =>'49');
my $species = "Homo sapiens";

#Get a connection to the local Ensembl CORE database
my $ensembl_core_api = Bio::EnsEMBL::Registry->get_DBAdaptor($species, "core");
my $slice_adaptor = $registry->get_adaptor($species, 'Core', 'Slice');

foreach my $clone_count (sort {$a <=> $b} keys %{$clones_ref}){

  my @starts = @{$clones_ref->{$clone_count}->{chr_starts_a}};
  my @sizes = @{$clones_ref->{$clone_count}->{block_sizes_a}};

  my $chr = $clones_ref->{$clone_count}->{chr};
  my $chromosome;
  if ($chr =~ /^chr(.*)/){
    $chromosome = $1;
  }else{
    print RED, "\nChrmosome format: $chr not understood!\n\n", RESET;
    exit();
  }

  #Build an mRNA sequence for the clone from the reference using the chromosome coordinates of each exon boundary
  my $mrna_seq = '';
  foreach my $start (@starts){
    my $size = shift (@sizes);
    my $end = $start+$size;

    #Get a slice for this exon
    my $slice = $slice_adaptor->fetch_by_region('chromosome', $chromosome, $start+1, $end);
    my $exon_seq = $slice->seq();

    $mrna_seq = "$mrna_seq"."$exon_seq";
  }

  print ">$clones_ref->{$clone_count}->{name}\n$mrna_seq\n";

}





exit();


##############################################################################################################
#Parse input file                                                                                            #
##############################################################################################################
sub parseExonCoordFile{
  my %args = @_;
  my $infile = $args{'-infile'};

  open (IN, "$infile") || die "\nCould not open input file: $infile\n\n";

  my %clones;

  my $clone_count = 0;
  while(<IN>){
    #Skip comment lines
    if ($_ =~ /^\#/){
      next();
    }
    chomp($_);
    my @line = split ("\t", $_);

    #Skip empty lines
    unless ($line[0]){
      next();
    }
    $clone_count++;

    my $block_sizes = $line[19];
    my $chr_starts = $line[21];

    my @block_sizes = split(",", $block_sizes);
    my @chr_starts = split(",", $chr_starts);


    $clones{$clone_count}{name} = $line[0];
    $clones{$clone_count}{chr} = $line[14];
    $clones{$clone_count}{block_sizes} = $line[19];
    $clones{$clone_count}{block_sizes_a} = \@block_sizes;
    $clones{$clone_count}{chr_starts} = $line[21];
    $clones{$clone_count}{chr_starts_a} = \@chr_starts;

  }

  close(IN);

  return(\%clones);
}





