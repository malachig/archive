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
my $position_file = '';
my $connect_type = '';
my $species = '';
my $out_file = '';
my $log_file = '';

GetOptions ('ensembl_api_version=s'=>\$ensembl_api_version,
	    'ensembl_database=s'=>\$ensembl_database, 'ensembl_server=s'=>\$ensembl_server,
	    'ensembl_user=s'=>\$ensembl_user, 'ensembl_password=s'=>\$ensembl_password,
	    'connect_type=s'=>\$connect_type, 'species=s'=>\$species,
	    'position_file=s'=>\$position_file, 'out_file=s'=>\$out_file, 'log_file=s'=>\$log_file);

#Provide instruction to the user
print BLUE, "\n\nNOTE: Before using this script, make sure the correct API version is hard coded!!\n\n", RESET;
print GREEN, "\n\tThis script takes a list of coordinates as input, gets the chromosomal sequences and writes them to a fasta file", RESET;
print GREEN, "\n\tA list of all genes overlaping these coordinates is also written", RESET;
print BLUE, "\n\nBasic required parameters:", RESET;
print GREEN, "\n\tSpecify the correct EnsEMBL API version using: --ensembl_api_version (41, 42, etc.)", RESET;
print GREEN, "\n\tSpecify the desired connection type using: --connect_type (see below for examples)", RESET;
print GREEN, "\n\tSpecify the input file containing mutation positions of interest", RESET;
print GREEN, "\n\t\tThis file contain mutation positions of the form: EnsEMBL_Gene_Name\tOther_Gene_Name\tchromosome:position\tMutation", RESET;
print GREEN, "\n\t\tEach line should look something like this: ABCB11\tABCB11\t2:169497197\tR>H", RESET;
print GREEN, "\n\t\tOnly the position and EnsEMBL gene name are really critical here, the other two columns are just descriptive", RESET;
print GREEN, "\n\tSpecify output file using: --out_file", RESET;
print GREEN, "\n\tSpecify a log file using: --log_file", RESET;
print GREEN, "\n\nExample: chrCoordsToOrfCoords.pl  --ensembl_api_version=54  --connect_type=local  --ensembl_server=ensembl01.bcgsc.ca  --ensembl_user=ensembl  --ensembl_password=ensembl  --species=Human  --position_file=Nature_paper_mutations.txt  --out_file=Nature_paper_mutations_ORF.txt  --log_file=Nature_paper_mutations_LOG.txt\n", RESET;

unless (($ensembl_api_version =~ /^\d+/) && $position_file && $connect_type && $out_file && $log_file){
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
  }elsif($ensembl_api_version eq "49"){
    unshift(@INC, "/home/malachig/perl/ensembl_49_perl_API/ensembl/modules");
  }elsif($ensembl_api_version eq "50"){
    unshift(@INC, "/home/malachig/perl/ensembl_50_perl_API/ensembl/modules");
  }elsif($ensembl_api_version eq "51"){
    unshift(@INC, "/home/malachig/perl/ensembl_51_perl_API/ensembl/modules");
  }elsif($ensembl_api_version eq "52"){
    unshift(@INC, "/home/malachig/perl/ensembl_52_perl_API/ensembl/modules");
  }elsif($ensembl_api_version eq "53"){
    unshift(@INC, "/home/malachig/perl/ensembl_53_perl_API/ensembl/modules");
  }elsif($ensembl_api_version eq "54"){
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


#*********************************************************************************************************
#2.) Get the regions of interest from the input file

#Check the input file
unless (-e $position_file){
  print RED, "\nCould not find input file: $position_file\n\n", RESET;
  exit();
}
my $result = `cat -A $position_file`;
if ($result =~ /\^M\$/g){
  print RED, "\n\nThis file appears to be from windows... you should run dos2unix on it before proceeding!\n\n", RESET;
  exit();
}

my %positions;
open (POS, "$position_file") || die "\nCould not open segment file: $position_file\n\n";
my $count = 0;
while (<POS>){
  $count++;
  chomp($_);
  if ($_ =~ /\d+/){
    my @line = split ("\t", $_);
    $positions{$count}{gene_name} = $line[0];
    $positions{$count}{original_name} = $line[1];
    if ($line[2] =~ /(\w+)\:(\d+)/){
      $positions{$count}{chr} = $1;
      $positions{$count}{start} = $2;
    } 
    $positions{$count}{mutation} = $line[3];
  }
}
close (POS);

open (OUT, ">$out_file") || die "\nCould not open output file: $out_file\n\n";
open (LOG, ">$log_file") || die "\nCould not open log file: $log_file\n\n";
print BLUE "\n\nPrinting output to screen and to log file: $log_file\n\n", RESET;

print OUT "gene_name\tgenomic_position\tmutation\ttrans_id\tprotein_length\tamino_acid_position\tpeptide_1\tpeptide_3\tpeptide_11\n";

foreach my $pos (sort {$a <=> $b} keys %positions){
  my $chr = $positions{$pos}{chr};
  my $position_chr_start = $positions{$pos}{start};;
  my $source_name = $positions{$pos}{gene_name};
  my $original_name = $positions{$pos}{original_name};

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
  my $slice = $slice_adaptor->fetch_by_region('chromosome', $chr, $position_chr_start, $position_chr_start);
  my @genes = @{$slice->get_all_Genes()};

  print LOG "\nTarget Genomic Coord: $chr:$position_chr_start (Fixed Name: $source_name  Original Name: $original_name).  Mutation: $positions{$pos}{mutation}";
  print BLUE, "\nTarget Genomic Coord: $chr:$position_chr_start (Fixed Name: $source_name  Original Name: $original_name).  Mutation: $positions{$pos}{mutation}", RESET;

  foreach my $gene (@genes){

    $count++;
    my %gene_info = %{&feature_details('-feature_obj'=>$gene)};
    my $stable_id = $gene->stable_id();

    #Get a new slice using the Gene ID
    $slice = $slice_adaptor->fetch_by_gene_stable_id($stable_id, 0);
    my @genes = @{$slice->get_all_Genes()};

    foreach my $gene (@genes){
      my $test_name = $gene->external_name();

      unless ($gene->is_known()){
        $test_name = "Unknown";
      }

      #Only continue if the found gene name matches the target gene name
      unless ($test_name eq $source_name){
        next();
      }

      my $gene_description = $gene->description();
      unless ($gene_description){
        $gene_description = "NA";
      }

      #Get chromosome coordinates for this gene
      my $new_gene = $gene->transform('chromosome');

      my $chr_test = $new_gene->slice->seq_region_name();
      my $chr_start_test = $new_gene->start();
      my $chr_end_test = $new_gene->end();
      my $strand_test = $new_gene->strand();
 
      #Get the coordinate of the mutation within the gene
      my $genic_start = ($position_chr_start - $chr_start_test)+1;
      print LOG "\n\tGene-level coordinate for target genomic coordinate: $genic_start";
      print YELLOW, "\n\tGene-level coordinate for target genomic coordinate: $genic_start", RESET;
     
      #Coordinates relative to whole chromosome
      print LOG "\n\tGene found: $stable_id ($test_name)\t$chr_test ($strand_test):$chr_start_test - $chr_end_test)";
      print LOG "\n\tGene Coords:", $gene->start(), " - ", $gene->end();
      print YELLOW, "\n\tGene found: $stable_id ($test_name)\t$chr_test ($strand_test):$chr_start_test - $chr_end_test)", RESET;
      print YELLOW, "\n\tGene Coords:", $gene->start(), " - ", $gene->end(), RESET;

      #Get transcript 
      my %transcripts;
      foreach my $trans (@{$gene->get_all_Transcripts()}){
        my $trans_id = $trans->stable_id();
        my $trans_start = $trans->start();
        my $trans_end = $trans->end();
        my $coding_region_start = $trans->coding_region_start();
        my $coding_region_end = $trans->coding_region_end();

        #Make sure the target coordinate is within this transcript!
        unless ($genic_start >= $coding_region_start && $genic_start <= $coding_region_end){
          next();
        }
        #It must also be within an exon of this transcript!
        my $exon_overlap = 0;
        foreach my $exon (@{$trans->get_all_Exons()}){
          my $exon_id = $exon->stable_id();
          my $exon_details_ref = &feature_details('-feature_obj'=>$exon);
          my $exon_start = $exon_details_ref->{$exon_id}->{start};
          my $exon_end = $exon_details_ref->{$exon_id}->{end};
          if ($genic_start >= $exon_start && $genic_start <= $exon_end){
            $exon_overlap = 1;
          }
        }
        unless ($exon_overlap == 1){
          next();
        }

        print LOG "\n\t\tTranscript: $trans_id\tCodingRegion: $coding_region_start - $coding_region_end";
        print YELLOW, "\n\t\tTranscript: $trans_id\tCodingRegion: $coding_region_start - $coding_region_end", RESET;
        my $trmapper = Bio::EnsEMBL::TranscriptMapper->new($trans);
        my ($pep_coords) = $trans->genomic2pep($genic_start,$genic_start, $strand_test);
        my $pep_coord = $pep_coords->start;
        my $peptide_seq = $trans->translate()->seq();
        my $peptide_1 = substr($peptide_seq, ($pep_coord-1), 1);
        my $peptide_3 = substr($peptide_seq, ($pep_coord-2), 3);
        my $peptide_11 = substr($peptide_seq, ($pep_coord-6), 11);
        my $protein_length = length($peptide_seq);
        print LOG "\n\t\t\tProtein length = $protein_length";
        print LOG "\n\t\t\tAmino acid position: $pep_coord ($peptide_1) ($peptide_3) ($peptide_11)";
        print YELLOW, "\n\t\t\tProtein length = $protein_length", RESET;
        print YELLOW, "\n\t\t\tAmino acid position: $pep_coord ($peptide_1) ($peptide_3) ($peptide_11)", RESET;
 
        print OUT "$original_name\tchr$chr:$position_chr_start\t$positions{$pos}{mutation}\t$trans_id\t$protein_length\t$pep_coord\t$peptide_1\t$peptide_3\t$peptide_11\n";

      }
    }
  }
  print LOG "\n";
  print "\n";

}
close(OUT);
close(LOG);

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
