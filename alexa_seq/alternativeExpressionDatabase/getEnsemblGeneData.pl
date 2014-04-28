#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#The purpose of this script is to get gene/transcript/exon data from ensembl so that it can be populated to an ALEXA database for probe design
#The Ensembl perl API is used for all database interactions
#Ensembl genes will be retrieved as slices with respect to the genes themselves (and have coordinates will be relative to this genomic slice containing the gene)
#Each ensembl gene may have multiple transcripts.  These transcripts will consist of combinations of exons (some which overlap multiple transcripts and others
#that are unique to a single transcript).  Some exons may overlap exons from another transcript for the same gene.
#The full gene sequence will start with the first exon (from any transcript) and go until the end of the last exon (from any transcript)
#All genes from Ensembl will be populated to the database.  Further analysis will have to consider whether to include pseudogenes, mtRNA, miRNA, etc.

#NOTE: Before running this script you must have access to a local copy of an EnsEMBL data you which to use for probe design
#      - You must also have access to an empty ALEXA database which has been created.
#      - A new ALEXA database should be created for each species or EnsEMBL version to be used for probe design.

use DBI;
use strict;
use Data::Dumper;
use Getopt::Long;
use Benchmark;
use Term::ANSIColor qw(:constants);
use File::Basename;

#Load the ALEXA modules
BEGIN {
  use Cwd 'abs_path';
  if (abs_path($0) =~ /(.*)\/.*\/.*\.pl/){
    push (@INC, $1);
  }
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
my $alexa_database = '';
my $alexa_server = '';
my $alexa_user = '';
my $alexa_password = '';
my $single_id = '';
my $all_ids = '';
my $populate_database = '';
my $logfile = '';
my $verbose = '';

GetOptions ('ensembl_api_version=s'=>\$ensembl_api_version, 'species=s'=>\$species,
	    #Source database access
	    'ensembl_database=s'=>\$ensembl_database, 'ensembl_server=s'=>\$ensembl_server, 'ensembl_user=s'=>\$ensembl_user, 'ensembl_password=s'=>\$ensembl_password,
	    #Target database access
	    'alexa_database=s'=>\$alexa_database, 'alexa_server=s'=>\$alexa_server, 'alexa_user=s'=>\$alexa_user, 'alexa_password=s'=>\$alexa_password,
	    'single_id=s'=>\$single_id, 'all_ids=i'=>\$all_ids, 'populate_database=s'=>\$populate_database, 'logfile=s'=>\$logfile,
	    'verbose=s'=>\$verbose);

#Provide instruction to the user
print GREEN, "\n\nNOTE: Before using this script, make sure the correct API version is hard coded!!\n\n", RESET;
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the correct EnsEMBL API version using: --ensembl_api_version (41, 42, etc.)", RESET;
print GREEN, "\n\tSpecify the SOURCE Ensembl Database, Server, User and Password using: --ensembl_database  --ensembl_server  --ensembl_user  and  --ensembl_password", RESET;
print GREEN, "\n\tThe species using: --species (e.g. --species=Human or --species='Homo sapiens')", RESET;
print GREEN, "\n\t\tMake sure the species you supply matches the EnsEMBL database you supply!!", RESET;
print GREEN, "\n\tSpecify the TARGET Database and Server to query using: --alexa_database and --alexa_server", RESET;
print GREEN, "\n\tSpecify the User and Password for access using: --alexa_user and --alexa_password\n", RESET;
print GREEN, "\n\tTo test, specify a single id using: --single_id (e.g. --single_id=ENSG00000171163  OR  --single_id=\"ENSG00000225151,ENSG00000240886\")", RESET;
print GREEN, "\n\tTo process all genes use: --all_ids=1", RESET;
print GREEN, "\n\tAfter testing, import data to specified ALEXA database using: --populate_database=yes", RESET;
print GREEN, "\n\tSpecify a logfile for output using: --logfile", RESET;
print GREEN, "\n\nExample: getEnsemblGeneData.pl  --ensembl_api_version=53  --species=Human  --ensembl_database=homo_sapiens_core_53_36o  --ensembl_server=source_server  --ensembl_user=ensembl_user  --ensembl_password=pwd  --alexa_database=ALEXA_hs_53_36o  --alexa_server=target_server  --alexa_user=user  --alexa_password=pwd  --all_ids=1  --populate_database=no  --logfile=/home/user/alexa/ALEXA_version/logs/database_population/getEnsemblGeneData_LOG.txt\n\n", RESET;

unless ($ensembl_api_version && $species && $ensembl_database && $ensembl_server && $ensembl_user && defined($ensembl_password) && $alexa_database && $alexa_server && $alexa_user && $alexa_password && $populate_database && $logfile){
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  exit();
}

#Load the specified ensembl API
&loadEnsemblApi('-api'=>$ensembl_api_version);

unless ($single_id || $all_ids){
  print RED, "\n\tSpecify --single_id or --all_ids\n\n", RESET;
  print RED, "\tTo test on a single ensembl ID: use --single_id=valid_ensembl_id\n", RESET;
  print RED, "\te.g. --single_id=ENSG00000171163\n", RESET;
  print RED, "\tUse --all_ids=1 to process all ensembl genes\n\n", RESET;
  exit();
}
unless ($populate_database eq "yes"){
  print YELLOW, "Once you have completed testing: Use -all_ids=1 and --populate_database=yes to process all genes and insert to database\n", RESET;
}

#Open logfile for output
open (LOG, ">$logfile") || die "\nCould not open logfile: $logfile\n\n";

print LOG "\nUser Specified the following options:\nspecies = $species\nalexa_database = $alexa_database\nsingle_id = $single_id\nall_ids = $all_ids\npopulate_database = $populate_database\nlogfile = $logfile\n\n";

#1.) Establish connections to source EnsEMBL database and target ALEXA databases

#A.) Using a local ensembl version

##CONNECT TO ENSEMBL SERVER
my $registry = "Bio::EnsEMBL::Registry";
$registry->load_registry_from_db(-host =>$ensembl_server, -user =>$ensembl_user, -pass=>$ensembl_password, -db_version =>$ensembl_api_version);

#Get a connection to the local Ensembl CORE database
my $ensembl_core_api = Bio::EnsEMBL::Registry->get_DBAdaptor($species, "core");

my $legacy = 0;
unless($ensembl_core_api){
  print YELLOW, "\nLoad registry failed.  Try legacy DBAdaptor\n", RESET;
  $ensembl_core_api = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => $ensembl_server,
  						    -user => $ensembl_user,
  						    -dbname => $ensembl_database,
						    -pass => $ensembl_password);
  $legacy=1;
}



#Get a connection to the local Ensembl VARIATION database
#my $ensembl_variation_api = Bio::EnsEMBL::Registry->get_DBAdaptor($species, "variation");
#my $vf_adaptor = $ensembl_variation_api->get_VariationFeatureAdaptor();

#B.) Alternatively you might try getting a connection to the Ensembl server itself over the internet
#    - This will be much slower but you can try something like this:
#  require Bio::EnsEMBL::Registry;  #Use for remote connections over the web
#  $ensembl_core_api = 'Bio::EnsEMBL::Registry';
#  $ensembl_core_api ->load_registry_from_db(-host=>'ensembldb.ensembl.org', -user=>'anonymous');

#2.) Establish connection with the Alternative Splicing Expression database using details provided by user
my $alexa_dbh = &connectDB('-database'=>$alexa_database, '-server'=>$alexa_server, '-user'=>$alexa_user, '-password'=>$alexa_password);

#Make sure the max allowed packet size is large enough to allow large genes
&checkMaxPacketSize('-dbh'=>$alexa_dbh);

my $ensembl_ids_ref;
if ($single_id){
  my @ids = split(",", $single_id);
  foreach my $id (@ids){
    $ensembl_ids_ref->{$id}->{gene_type} = "TESTING";
  }
}elsif($all_ids eq "1"){
  #Get all ensembl IDs for processing
  $ensembl_ids_ref = &getAllEnsemblIds();
}

my $source = $ensembl_database;

my $gene_count = 0;
my $haplotype_gene_count = 0;
my $known_gene_count = 0;
my $unknown_gene_count = 0;
my $gene_update_count = 0;
my $chr_strand;
my $total_bases = 0;
my $masked_bases = 0;
my $multi_gene_exons = 0;

foreach my $ensembl_id (sort keys %{$ensembl_ids_ref}){
  my $t0 = new Benchmark;
  my %gene_info;
  my $gene_info_ref = \%gene_info;

  #A Slice object represents a continuous region of a genome
  my $slice_adaptor;
  if ($legacy){
    $slice_adaptor = $ensembl_core_api->get_SliceAdaptor(); #Alternate way of getting a 'slice adaptor'
  }else{
    $slice_adaptor = $registry->get_adaptor($species, 'Core', 'Slice');
  }

  #Retrieve a slice with respect to a gene, with a specified flanking sequence on either side
  my $slice = $slice_adaptor->fetch_by_gene_stable_id($ensembl_id, 0);

  #BEFORE DOING ANYTHING: Check the strand of the target gene within this strand.  If it is on the reverse strand - invert the slice and continue
  #To simplify matters we always want the slice to be oriented so that the first exon (in terms of transcription, coding) at the beginning of the slice
  $chr_strand = 1; #Initialize to 1, switch to -1 if gene is on reverse strand.
  $slice = check_strand('-slice_obj'=>$slice, '-ensembl_id'=>$ensembl_id);

  #Determine the coordinate system used here - should be chromosome for most genomes
  #These coordinates can be used to map the location of genes, transcripts, and exons to their chromosomal positions
  my $coord_sys = $slice->coord_system()->name();
  my $chr = $slice->seq_region_name();

  #NOTE:  exclude genes from the 'haplotype' or 'alternate' chromosomes.  These are redundant annotations of those on the reference chromosomes and complicate downstream analyses!!
  if ($chr eq "c6_COX" || $chr eq "c6_QBL" || $chr eq "c5_H2"){
    #NCBI Build 36
    print YELLOW, "\n\tFeature is defined on a haplotype/alternate chromosome ($chr) - skipping", RESET;
    print LOG "\n\tFeature is defined on a haplotype/alternate chromosome ($chr) - skipping";
    $haplotype_gene_count++;
    next();
  }elsif($chr eq "HSCHR6_MHC_APD" || $chr eq "HSCHR6_MHC_COX" || $chr eq "HSCHR6_MHC_DBB" || $chr eq "HSCHR6_MHC_MANN" || $chr eq "HSCHR6_MHC_MCF" || $chr eq "HSCHR6_MHC_QBL" || $chr eq "HSCHR6_MHC_SSTO" || $chr eq "HSCHR4_1" || $chr eq "HSCHR17_1"){

    #Genome Reference Consortium - GRCh37
    print YELLOW, "\n\tFeature is defined on a haplotype/alternate chromosome ($chr) - skipping", RESET;
    print LOG "\n\tFeature is defined on a haplotype/alternate chromosome ($chr) - skipping";
    $haplotype_gene_count++;
    next();
  }


  #Get chromosome coordinates - if the gene is on a 'supercontig' (e.g. NT_XXXXXX) then these coordinates will be relative to that
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
    $gene_count++;

    my $test_name = $gene->external_name();
    if ($gene->is_known()){
      $known_gene_count++;
    }else{
      $test_name = "Unknown";
      $unknown_gene_count++;
    }
    my $gene_description = $gene->description();
    unless ($gene_description){
      $gene_description = "NA";
    }
    #Clean up description string.  Some heretic decided it was a good idea to put things like ' in these...
    $gene_description =~ s/'/\\'/g;
  
    print BLUE, "\n$gene_count -> Gene NAME: $test_name ($ensembl_id)", RESET;
    print LOG "\n$gene_count -> Gene NAME: $test_name ($ensembl_id)";


    #Transform the coordinates to chromosome system
    if (my $new_gene = $gene->transform('chromosome')){

      my $chr_test = $new_gene->slice->seq_region_name();
      my $chr_start_test = $new_gene->start();
      my $chr_end_test = $new_gene->end();
      my $strand_test = $new_gene->strand();

      print BLUE, "\n\tChr Coordinates: $strand_test $chr_test ($chr_start_test - $chr_end_test)", RESET;
      print LOG "\n\tChr Coordinates: $strand_test $chr_test ($chr_start_test - $chr_end_test)";

    }else{
      print YELLOW, "\n\tFeature is not defined in chromosome coordinate system", RESET;
      print LOG "\n\tFeature is not defined in chromosome coordinate system";
    }

    #Get all useful gene info for this ensembl id (HUGO, refseq, GO, OMIM, etc.)
    &getBasicGeneInfo('-ensembl_gene_obj'=>$gene, '-ensembl_id'=>$ensembl_id, '-gene_object'=>$gene_info_ref);

    #Get the entire sequence for this gene (from start of 1st exon to end of last exon, including introns)
    my $full_gene_sequence = $slice->seq();
    $gene_info_ref->{$ensembl_id}->{sequence} = $full_gene_sequence;

    #Get the entire sequence in hard-masked form (N's representing all bases identified by RepeatMasker, Tandem Repeat Finder, and Dust)
    my $hard_masked_seq_object = $slice->get_repeatmasked_seq();
    my $hard_masked_seq = $hard_masked_seq_object->seq();
    my $upper_masked_seq = uc($hard_masked_seq);
    $gene_info_ref->{$ensembl_id}->{masked_sequence} = $upper_masked_seq;
    $total_bases += length($upper_masked_seq);
    $masked_bases += ($upper_masked_seq =~ tr/N/N/);

    #Sanity check of unmasked vs. masked sequences
    unless (length($full_gene_sequence) == length($upper_masked_seq)){
      print RED, "\nLength of masked and unmasked sequences are not equal!!\n\n", RESET;
      exit();
    }

    #Get feature details for this gene
    my $gene_details_ref = &feature_details('-feature_obj'=>$gene);
    my $gene_start = $gene_details_ref->{$gene_id}->{start};
    my $gene_end = $gene_details_ref->{$gene_id}->{end};
    $gene_info_ref->{$ensembl_id}->{chr_strand} = $chr_strand;
    $gene_info_ref->{$ensembl_id}->{corrected_strand} = $gene_details_ref->{$gene_id}->{strand};
    $gene_info_ref->{$ensembl_id}->{chromosome} = $gene_details_ref->{$gene_id}->{chromosome};
    $gene_info_ref->{$ensembl_id}->{chr_start} = $chr_start;
    $gene_info_ref->{$ensembl_id}->{chr_end} = $chr_end;
    $gene_info_ref->{$ensembl_id}->{gene_start} = $gene_start;
    $gene_info_ref->{$ensembl_id}->{gene_end} = $gene_end;
    $gene_info_ref->{$ensembl_id}->{description} = $gene_description;

    #Get all transcripts for this gene
    my %transcripts;
    my @trans_list = @{$gene->get_all_Transcripts()};
    my $trans_count = scalar(@trans_list);
    print BLUE, "\n\tFound $trans_count transcripts for this gene", RESET;

    #foreach my $trans (@{$gene->get_all_Transcripts()}){
    foreach my $trans (@trans_list){

      my $trans_id = $trans->stable_id();

      my $trans_details_ref = &feature_details('-feature_obj'=>$trans);
      $transcripts{$trans_id}{start} = $trans_details_ref->{$trans_id}->{start};
      $transcripts{$trans_id}{end} = $trans_details_ref->{$trans_id}->{end};

      my %exons;
      foreach my $exon (@{$trans->get_all_Exons()}){
	my $exon_id = $exon->stable_id();
        my $exon_details_ref = &feature_details('-feature_obj'=>$exon);
        my $exon_start = $exon_details_ref->{$exon_id}->{start};
        my $exon_end = $exon_details_ref->{$exon_id}->{end};
        &checkCoords('-g_start'=>$gene_start ,'-g_end'=>$gene_end, '-f_start'=>$exon_start, '-f_end'=>$exon_end);
        $exons{$exon_id}{start} = $exon_start;
	$exons{$exon_id}{end} = $exon_end;
      }
      $transcripts{$trans_id}{exons} = \%exons;

      #Get all useful IDs for every transcript of every gene
      my %external_ids;
      foreach my $link (@{$trans->get_all_DBLinks()}) {
	my $link_id = $link->display_id;
	my $link_database = $link->database;

	#Eliminate redundancies if they exist (shouldnt be redundant at transcript level)
	$external_ids{$link_id}{type} = $link_database;
      }

      $transcripts{$trans_id}{external_ids} = \%external_ids;

      #Get CDS start and end info for each transcript (get genomic coordinates, not those relative to the cDNA)
      my $coding_region_start = $trans->coding_region_start();
      my $coding_region_end = $trans->coding_region_end();

      #Pseudogenes and RNA genes will not have a cDNA start and end.  Set to 'na'
      my @cds_starts;
      my @cds_ends;
      if ($coding_region_start && $coding_region_end){
	
	my $cdna_start = $trans->cdna_coding_start();
	my $cdna_end = $trans->cdna_coding_end();

	my $trmapper = Bio::EnsEMBL::TranscriptMapper->new($trans);
	my @genomic_coords = $trmapper->cdna2genomic($cdna_start,$cdna_end);

	foreach my $gc (@genomic_coords){
	  my $gc_start = $gc->start();
	  my $gc_end = $gc->end();
	  push(@cds_starts,$gc_start);
	  push(@cds_ends, $gc_end);
	}

      }else{
	$coding_region_start = 'na';
	$coding_region_end = 'na';
	push(@cds_starts, 'na');
	push(@cds_ends, 'na');
      }

      $transcripts{$trans_id}{coding_region_start} = $coding_region_start;
      $transcripts{$trans_id}{coding_region_end} = $coding_region_end;
      $transcripts{$trans_id}{cds_starts} = \@cds_starts;
      $transcripts{$trans_id}{cds_ends} = \@cds_ends;

      #Get protein feature info for each transcript (things like protein domains, transmembrane domains, etc.)
      my %protein_features;
      my $pf_count = 0;

      #Note that not all transcripts will have a translation (pseudogenes, rna genes, etc.).
      my $translation = $trans->translation();

      #If a translation is not defined then do not search for protein features
      if($translation){
	my $protein_features = $translation->get_all_ProteinFeatures();

	foreach my $pf (@$protein_features){
	  $pf_count++;
	  my $logic_name = $pf->analysis()->logic_name();
	  my $pf_start = $pf->start();
	  my $pf_end = $pf->end();
	  my $interpro_ac = $pf->interpro_ac();     #not always defined
	  my $idesc = $pf->idesc();                 #not always defined
	  my $program = $pf->analysis()->program(); #not always defined

	  unless ($interpro_ac){$interpro_ac = 'na'};
	  unless ($idesc){$idesc = 'na'};
	  unless ($program){$program = 'na'};

	  #Convert protein coordinates to gene coordinates for storage
	  #This will return coordinate for each exon involved in this protein domain.
	  #The start and end will be stored, but all coordinates will also be stored for convenience
	  my $trmapper = Bio::EnsEMBL::TranscriptMapper->new($trans);
	  my @genomic_coords = $trmapper->pep2genomic($pf_start,$pf_end);

	  my $pf_gene_start = $genomic_coords[0]->start();
	  my $pf_gene_end = $genomic_coords[0]->end();
	  my @pf_starts;
	  my @pf_ends;
	  foreach my $gc (@genomic_coords){
	    my $gc_start = $gc->start();
	    my $gc_end = $gc->end();

	    if ($gc_start < $pf_gene_start){
	      $pf_gene_start = $gc_start;
	    }
	    if ($gc_end > $pf_gene_end){
	      $pf_gene_end = $gc_end;
	    }
	    push(@pf_starts,$gc_start);
	    push(@pf_ends, $gc_end);
	  }

	  $protein_features{$pf_count}{logic_name} = $logic_name;
	  $protein_features{$pf_count}{interpro_ac} = $interpro_ac;
	  $protein_features{$pf_count}{idesc} = $idesc;
	  $protein_features{$pf_count}{program} = $program;
	  $protein_features{$pf_count}{start} = $pf_gene_start;
	  $protein_features{$pf_count}{end} = $pf_gene_end;
	  $protein_features{$pf_count}{start_coords} = \@pf_starts;
	  $protein_features{$pf_count}{end_coords} = \@pf_ends;
	}
      }else{
	print YELLOW, "\n\t\tNOTE: Trans ($trans_id) of Gene ($ensembl_id) does not have a translation", RESET;
	print LOG "\n\t\tNOTE: Gene $ensembl_id does not have a translation";
      }
      $transcripts{$trans_id}{protein_features} = \%protein_features;
    }
    $gene_info_ref->{$ensembl_id}->{trans} = \%transcripts;

  }

  #Populate the database with this gene record if the user specified the -d flag
  if ($populate_database eq "yes"){
    $gene_update_count++;
    &populateDatabase('-gene_object'=>$gene_info_ref);
    my $t1 = new Benchmark;
    my $td1 = timediff($t1, $t0);
    print BLUE, "\n\n$gene_update_count Genes Inserted to Database\tElapsed Time:",timestr($td1),"\n", RESET;
    print LOG "\n\n$gene_update_count Genes Inserted to Database\tElapsed Time:",timestr($td1),"\n";
  }
}

#Close database connection
$alexa_dbh->disconnect();

my $percent_masked = sprintf("%.3f",(($masked_bases/$total_bases)*100));
print BLUE, "\n\nEnsEMBL reports the following as being identified by RepeatMasker, TRF or Dust:", RESET;
print BLUE, "\nFound a total of $masked_bases masked bases of $total_bases total bases ($percent_masked%)", RESET;
print LOG "\n\nEnsEMBL reports the following as being identified by RepeatMasker, TRF or Dust:";
print LOG "\nFound a total of $masked_bases masked bases of $total_bases total bases ($percent_masked%)";

print BLUE, "\n\nProcessed a total of $gene_count genes: $known_gene_count are KNOWN, $unknown_gene_count are UNKNOWN (predicted only)\n", RESET;
print LOG "\n\nProcessed a total of $gene_count genes: $known_gene_count are KNOWN, $unknown_gene_count are UNKNOWN (predicted only)\n";

my $total_gene_count = $gene_count + $haplotype_gene_count;
print BLUE, "\n\nDatabase should contain $gene_count genes.  $haplotype_gene_count genes from haplotype chromosome were skipped.  Total EnsEMBL genes considered = $total_gene_count (check against online database)\n", RESET;
print LOG "\n\nDatabase should contain $gene_count genes.  $haplotype_gene_count genes from haplotype chromosome were skipped.  Total EnsEMBL genes considered = $total_gene_count (check against online database)\n";

print BLUE, "\n\nFound $multi_gene_exons exons that appear to be defined for multiple genes but have the same EnsEMBL exon ID\n", RESET;
print LOG "\n\nFound $multi_gene_exons exons that appear to be defined for multiple genes but have the same EnsEMBL exon ID\n";

print "\n\n";
print LOG "\n\n";
close (LOG);

exit();


##############################################################################################
#Get all ensembl Ids from the local ensembl database that are appropriate                    #
##############################################################################################
sub getAllEnsemblIds{
  my %ids;

  #Get a normal DBI connection to the ensembl database
  my $ensembl_dbh = &connectDB('-database'=>$ensembl_database, '-server'=>$ensembl_server, '-user'=>$ensembl_user, '-password'=>$ensembl_password);

  #Get the gene sequence for this gene
  my $sql = "SELECT gene.gene_id,gene_stable_id.stable_id,gene.biotype FROM gene,gene_stable_id WHERE gene_stable_id.gene_id = gene.gene_id;";
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


###################################################################################
#For a given target gene, determine the strand of this gene in its genomic slice  #
#If it is on the reverse strand, invert the entire slice and return it            #
###################################################################################
sub check_strand{
  my %args = @_;
  my $slice = $args{'-slice_obj'};
  my $ensembl_id = $args{'-ensembl_id'};

  foreach my $gene (@{$slice->get_all_Genes()}){

    #Make sure the right gene is being considered and not an overlapping gene (which will be processed individually)
    my $gene_id = $gene->stable_id();

    unless ($gene_id eq $ensembl_id){
      next();
    }

    #If the gene is on the reverse strand with respect to the slice, invert the slice
    my $strand_test = $gene->strand();

    if ($strand_test eq "-1"){
      #print "\nFeature is on reverse strand - inverting\n\n";
      $slice = $slice->invert();
      $chr_strand = -1;
    }
  }
  return ($slice);
}

##################################################################################
#Get external descriptive information if a gene is known                         #
##################################################################################
sub getBasicGeneInfo{
  my %args = @_;

  my $gene = $args{'-ensembl_gene_obj'};
  my $ensembl_id = $args{'-ensembl_id'};
  my $gene_object_ref = $args{'-gene_object'};

  my $gene_name = $gene->external_name();

  unless ($gene_name){
    $gene_name = "Unknown";
  }

  #Get the gene type
  my $gene_type = $gene->biotype();

  #Get the gene version
  my $ensembl_version = $gene->version();

  $gene_object_ref->{$ensembl_id}->{source} = $source;
  $gene_object_ref->{$ensembl_id}->{gene_type} = $gene_type;
  $gene_object_ref->{$ensembl_id}->{gene_name} = $gene_name;
  $gene_object_ref->{$ensembl_id}->{ensembl_version} = $ensembl_version;

  if ($gene->is_known()) {
    $gene_object_ref->{$ensembl_id}->{evidence} = "Known Gene";
  }else{
    $gene_object_ref->{$ensembl_id}->{evidence} = "Unknown Gene";
    print YELLOW, "\n\tNOTE: Gene " . $gene->stable_id() . " is not a known gene", RESET;
    print LOG "\n\tNOTE: Gene " . $gene->stable_id() . " is not a known gene";
  }

  return();
}

######################################################################
#Populate the database with Gene, Transcript, Exon info              #
######################################################################
sub populateDatabase{
  my %args = @_;
  my $gene_object = $args{'-gene_object'};

  foreach my $ensg_id (sort keys %{$gene_object}){

    #Clean up gene_name string.  Some heretic decided it was a good idea to put things like ' in these...
    my $gene_name = $gene_object->{$ensg_id}->{gene_name};
    $gene_name =~ s/'/\\'/g;

    #1a.)Insert gene into Gene table
    my $sql_gene = "INSERT INTO Gene (ensembl_g_id,ensembl_version,source,gene_type,gene_name,description,evidence,sequence,chr_strand,corrected_strand,chromosome,chr_start,chr_end,gene_start,gene_end) VALUES (\'$ensg_id\',\'$gene_object->{$ensg_id}->{ensembl_version}\',\'$gene_object->{$ensg_id}->{source}\',\'$gene_object->{$ensg_id}->{gene_type}\',\'$gene_name\',\'$gene_object->{$ensg_id}->{description}\',\'$gene_object->{$ensg_id}->{evidence}\',\'$gene_object->{$ensg_id}->{sequence}\',\'$gene_object->{$ensg_id}->{chr_strand}\',\'$gene_object->{$ensg_id}->{corrected_strand}\',\'$gene_object->{$ensg_id}->{chromosome}\',\'$gene_object->{$ensg_id}->{chr_start}\',\'$gene_object->{$ensg_id}->{chr_end}\',\'$gene_object->{$ensg_id}->{gene_start}\',\'$gene_object->{$ensg_id}->{gene_end}\');";

    my $sth_gene = $alexa_dbh->prepare("$sql_gene");
    $sth_gene->execute() || die "Couldn't execute sql_gene statement:\n\n" . $sth_gene->errstr;
    $sth_gene->finish();

    my $gene_id = $alexa_dbh->{'mysql_insertid'};

    #Sanity check, the returned gene ID (from auto increment) should match the current valid gene count
    unless ($gene_id == $gene_count){
      print RED, "\nGene ID ($gene_id) reported by mysql_insertid does not match that expected according to current count ($gene_count)\n\n", RESET;
      exit();
    }

    print BLUE, "\n->New Gene Record: $gene_id", RESET;
    print LOG "\n->New Gene Record: $gene_id";

    #1b.)Insert masked gene sequence into MaskedGene table
    my $masked_seq = $gene_object->{$ensg_id}->{masked_sequence};
    my $sql_masked = "INSERT INTO MaskedGene (fk_Gene__id, sequence) VALUES ($gene_id,\'$masked_seq\');";

    #Actually insert the record
    my $sth_masked = $alexa_dbh->prepare("$sql_masked");
    $sth_masked->execute() || die "Couldn't execute sql_masked statement:\n\n" . $sth_masked->errstr;
    $sth_masked->finish();
    my $masked_gene_id = $alexa_dbh->{'mysql_insertid'};

    print BLUE, "\n\t->New Masked Gene Record: $masked_gene_id", RESET;
    print LOG "\n\t->New Masked Gene Record: $masked_gene_id";

    #2.) Insert each transcript for this gene and relate it to the parent gene
    my $transcripts_ref = $gene_object->{$ensg_id}->{trans};
    foreach my $enst_id (sort keys %{$transcripts_ref}){

      #Make sure this transcript has not already been entered (look for one transcript linked to more than one ensembl gene)
      my $sql_test1 = "SELECT id FROM Transcript WHERE ensembl_t_id = '$enst_id'";
      my $sth_test1 = $alexa_dbh->prepare("$sql_test1");
      $sth_test1->execute() || die "Couldn't execute sql_test1 statement:\n\n" . $sth_test1->errstr;
      my $t_id = $sth_test1->fetchrow_array();
      $sth_test1->finish();

      if ($t_id){
	print RED, "\nWARNING: Transcript: $t_id was already in the database!\n\n", RESET;
	print LOG "\nWARNING: Transcript: $t_id was already in the database!\n\n";
      }

      my $sql_trans = "INSERT INTO Transcript (fk_Gene__id,ensembl_t_id,start,end,coding_region_start,coding_region_end,cds_start_coords,cds_end_coords) VALUES (\'$gene_id\',\'$enst_id\',\'$transcripts_ref->{$enst_id}->{start}\',\'$transcripts_ref->{$enst_id}->{end}\',\'$transcripts_ref->{$enst_id}->{coding_region_start}\',\'$transcripts_ref->{$enst_id}->{coding_region_end}\',\'@{$transcripts_ref->{$enst_id}->{cds_starts}}\',\'@{$transcripts_ref->{$enst_id}->{cds_ends}}\');";

      my $sth_trans = $alexa_dbh->prepare("$sql_trans");
      $sth_trans->execute() || die "Couldn't execute sql_trans statement:\n\n" . $sth_trans->errstr;
      $sth_trans->finish();

      my $trans_id = $alexa_dbh->{'mysql_insertid'};

      print BLUE, "\n\t->New Transcript Record: $trans_id", RESET;
      print LOG "\n\t->New Transcript Record: $trans_id";

      #3) Insert all the external IDs for this transcript (HUGO, Refseq, etc. etc.)
      my $external_id_count = 0;

      my $external_ids_ref = $transcripts_ref->{$enst_id}->{external_ids};

      #Get the number of external ID entries needed for this transcript
      my $num_ids = keys %{$external_ids_ref};

      my $sql1 = "INSERT INTO External_id (fk_Transcript__id,type,external_id) VALUES ";

      foreach my $external_id (sort keys %{$external_ids_ref}){
	$external_id_count++;
	my $external_id_type = $external_ids_ref->{$external_id}->{type};
        $external_id =~ s/'/\\'/g;

	#If this is the last external ID for this transcript, do the multiple insert
	if ($external_id_count == $num_ids){
	  #Add insert record onto end, but omit the comma since this is the last one of this block
	  $sql1 = "$sql1"."($trans_id,\'$external_id_type\',\'$external_id\');";

	  #Actually do the insert!
	  my $sth1 = $alexa_dbh->prepare("$sql1");
	  $sth1->execute() || die "Couldn't execute sql1 statement:\n\n" . $sth1->errstr;
	  $sth1->finish();

	  print BLUE, "\n\t\t->Entered $external_id_count link Ids", RESET;
	  print LOG "\n\t\t->Entered $external_id_count link Ids";

	}else{
	  $sql1 = "$sql1"."($trans_id,\'$external_id_type\',\'$external_id\'),";
	}
      }

      #4.) Insert all of the protein features for this transcript
      my $protein_feature_count = 0;

      my $protein_features_ref = $transcripts_ref->{$enst_id}->{protein_features};

      #Get the number of protein feature entries needed for this transcript
      my $num_features = keys %{$protein_features_ref};

      my $sql_pf = "INSERT INTO Protein_feature (fk_Transcript__id,type,interpro_ac,name,program,start,end,start_coords,end_coords) VALUES ";

      foreach my $pf_count (sort {$a <=> $b} keys %{$protein_features_ref}){
	$protein_feature_count++;
        my $logic_name = $protein_features_ref->{$pf_count}->{logic_name};
        my $idesc = $protein_features_ref->{$pf_count}->{idesc};
        $logic_name =~ s/'/\\'/g;
        $idesc =~ s/'/\\'/g;

	#If this is the last feature for this transcript, do the multiple insert
	if ($protein_feature_count == $num_features){
	  #Add insert record onto end, but omit the comma since this is the last one of this block

          #print Dumper $sql_pf;

	  $sql_pf = "$sql_pf"."($trans_id,\'$logic_name\',\'$protein_features_ref->{$pf_count}->{interpro_ac}\',\'$idesc\',\'$protein_features_ref->{$pf_count}->{program}\',\'$protein_features_ref->{$pf_count}->{start}\',\'$protein_features_ref->{$pf_count}->{end}\',\'@{$protein_features_ref->{$pf_count}->{start_coords}}\',\'@{$protein_features_ref->{$pf_count}->{end_coords}}\');";

	  #Actually do the insert!
	  my $sth_pf = $alexa_dbh->prepare("$sql_pf");
	  $sth_pf->execute() || die "Couldn't execute sql_pf statement:\n\n" . $sth_pf->errstr;
	  $sth_pf->finish();

	  print BLUE, "\n\t\t->Entered $protein_feature_count protein features", RESET;
	  print LOG "\n\t\t->Entered $protein_feature_count protein features";

	}else{
	  $sql_pf = "$sql_pf"."($trans_id,\'$logic_name\',\'$protein_features_ref->{$pf_count}->{interpro_ac}\',\'$idesc\',\'$protein_features_ref->{$pf_count}->{program}\',\'$protein_features_ref->{$pf_count}->{start}\',\'$protein_features_ref->{$pf_count}->{end}\',\'@{$protein_features_ref->{$pf_count}->{start_coords}}\',\'@{$protein_features_ref->{$pf_count}->{end_coords}}\'),";
	}
      }

      #5.) Insert each exon for this transcript and relate it to the parent transcript
      #Make sure it hasn't been entered already for a different transcript of this gene - this is expected
      my $exons_ref = $transcripts_ref->{$enst_id}->{exons};
      my $num_exons = keys %{$exons_ref};

      my $trans_exon_count = 0;

      my $sql2 = "INSERT INTO TranscriptExon (fk_Transcript__id,fk_Exon__id) VALUES ";

      foreach my $ense_id (sort keys %{$exons_ref}){
	$trans_exon_count++;
        my $exon_start = $exons_ref->{$ense_id}->{start};
        my $exon_end = $exons_ref->{$ense_id}->{end};

	#Check if this exon has already been entered for another transcript
        my $sql_test2 = "SELECT id,start,end FROM Exon WHERE ensembl_e_id = '$ense_id'";
	my $sth_test2 = $alexa_dbh->prepare("$sql_test2");
	$sth_test2->execute() || die "Couldn't execute sql_test2 statement:\n\n" . $sth_test2->errstr;
	my ($e_id, $stored_start, $stored_end) = $sth_test2->fetchrow_array();
	$sth_test2->finish();

	#If the exon is already there, simply create a new relationship to it's parent transcript
        #Make sure it is not the same ensembl_e_id being defined for a different gene!  This will cause coordinate problems
        #If this occurs, you need to make a new exon record with the same ensembl_e_id
        #As long as the start and end coordinates are the same for the stored exon and the current one - it should be okay
	if ($e_id){

          if ($exon_start == $stored_start && $exon_end == $stored_end){
	    #If this is the last exon for this transcript, finish up the sql multiple query (otherwise just add it on):
	    if($trans_exon_count == $num_exons){
	      $sql2 = "$sql2"."($trans_id,$e_id);";
	    }else{
	      $sql2 = "$sql2"."($trans_id,$e_id),";
	    }
          }else{
            print RED, "\n\t\t-> Found an EnsEMBL exon ID defined for multiple genes", RESET;
            $multi_gene_exons++;

            #Create a new exon record for this duplicate exon
	    my $sql_exon = "INSERT INTO Exon (ensembl_e_id,start,end) VALUES (\'$ense_id\',\'$exon_start\',\'$exon_end\');";
	    my $sth_exon = $alexa_dbh->prepare("$sql_exon");
	    $sth_exon->execute() || die "Couldn't execute sql_exon statement:\n\n" . $sth_exon->errstr;
	    $sth_exon->finish();
	    my $new_exon_id = $alexa_dbh->{'mysql_insertid'};
	    print BLUE, "\n\t\t->New Exon Record: $new_exon_id", RESET;
	    print LOG "\n\t\t->New Exon Record: $new_exon_id";
	    #Again, if this is the last exon for this transcript, finish up the sql multiple query (otherwise just add it on):
	    if($trans_exon_count == $num_exons){
	      $sql2 = "$sql2"."($trans_id,$new_exon_id);";
	    }else{
	      $sql2 = "$sql2"."($trans_id,$new_exon_id),";
	    }

          }

	}else{
	  #If the exon is not already entered, enter it now.
	  my $sql_exon = "INSERT INTO Exon (ensembl_e_id,start,end) VALUES (\'$ense_id\',\'$exon_start\',\'$exon_end\');";
	  my $sth_exon = $alexa_dbh->prepare("$sql_exon");
	  $sth_exon->execute() || die "Couldn't execute sql_exon statement:\n\n" . $sth_exon->errstr;
	  $sth_exon->finish();
	  my $new_exon_id = $alexa_dbh->{'mysql_insertid'};
	  print BLUE, "\n\t\t->New Exon Record: $new_exon_id", RESET;
	  print LOG "\n\t\t->New Exon Record: $new_exon_id";
	  #Again, if this is the last exon for this transcript, finish up the sql multiple query (otherwise just add it on):
	  if($trans_exon_count == $num_exons){
	    $sql2 = "$sql2"."($trans_id,$new_exon_id);";
	  }else{
	    $sql2 = "$sql2"."($trans_id,$new_exon_id),";
	  }

	}
      }
      #Actually insert all the TranscriptExon relationships for this transcript before going on to the next one.
      my $sth2 = $alexa_dbh->prepare("$sql2");
      $sth2->execute() || die "Couldn't execute sql2 statement:\n\n" . $sth2->errstr;
      $sth2->finish();

      print BLUE, "\n\t\t->Entered $trans_exon_count TranscriptExon Relationships", RESET;
      print LOG "\n\t\t->Entered $trans_exon_count TranscriptExon Relationships";
    }
  }
  return();
}


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

############################################################################
#checkCoords - make sure a featue is defined within the bounds of the gene #
############################################################################
sub checkCoords{
  my %args = @_;
  my $g_start = $args{'-g_start'};
  my $g_end = $args{'-g_end'};
  my $f_start = $args{'-f_start'};
  my $f_end = $args{'-f_end'};

  unless (($f_start >= $g_start) && ($f_start <= $g_end) && ($f_end >= $g_start) && ($f_end <= $g_end)){
    print RED, "\n\nFound a feature with coordinates that are outside the region of the gene!", RESET;
    exit();
  }
  
  return();
}


