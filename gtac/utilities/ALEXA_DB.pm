=head1 NAME

ALEXA_DB.pm - Basic methods for accessing the Alternative Splicing Expression Analysis Database (ALEXA)

=head1 SYNOPSIS

use ALEXA_DB qw(:all);

=head2 NOTE

Currently located in '~/Array_design/utilities'

=head2 RECENT CHANGES

Various modifications.  Last modified 21 December 2006

=head1 DESCRIPTION

Generic utility for accessing data in the Alternative Splicing Expression Database Analysis.

=head1 EXAMPLES

use lib './';

use utilities::ALEXA_DB qw(:all);

=head1 SEE ALSO

None

=head1 BUGS

Contact author via email

=head1 AUTHOR

Written by Malachi Griffith (malachig@bcgsc.ca)

=head1 ACKNOWLEDGEMENTS

University of British Columbia Graduate Studies

Michael Smith Foundation for Health Research

Natural Sciences and Engineering Research Council

Genome British Columbia

=head1 AFFLIATIONS

Malachi Griffith is supervised by Marco A. Marra

Genome Sciences Centre, BC Cancer Research Centre, BC Cancer Agency, UBC Faculty of Medicine - Medical Genetics

=head1 SUBROUTINES

=cut

package utilities::ALEXA_DB;
require Exporter;

@ISA = qw(Exporter);
@EXPORT = qw();

@EXPORT_OK = qw(&checkMaxPacketSize &getAllGenes &getGeneIds &getGeneInfo &getTranscripts &getProbeInfo &getGeneProbes &getNcProbes &getExons
		&getMaskedExons &getExonContent &getIntronContent &getMaskedGene &getGeneTerms &getProteinFeatures &junctionProbeCombinations &displayGeneStats
		&convertGeneCoordinates &convertGenomicCoordinates &relativeGenePosition &getMicroarrayProbes &selectFilteredSet);

%EXPORT_TAGS = (all => [qw(&checkMaxPacketSize &getAllGenes &getGeneIds &getGeneInfo &getTranscripts &getProbeInfo &getGeneProbes &getNcProbes &getExons
			   &getMaskedExons &getExonContent &getIntronContent &getMaskedGene &getGeneTerms &getProteinFeatures &junctionProbeCombinations &displayGeneStats 
			   &convertGeneCoordinates &convertGenomicCoordinates &relativeGenePosition &getMicroarrayProbes &selectFilteredSet)]);

use strict;
use Data::Dumper;
use DBI;
use Term::ANSIColor qw(:constants);


=head2 checkMaxPacketSize()

=over 3

=item Function:

Check the database max_allowed_packet variable to ensure that large gene entries will work

=item Return:

N/A - Will exit if packet size is not sufficient

=item Args:

'-dbh' => database handle

=item Example(s):

&checkMaxPacketSize('-dbh'=>$alexa_dbh);

=back

=cut

#############################################################################################
#checkMaxPacketSize                                                                         #
#############################################################################################
sub checkMaxPacketSize{
  my %args = @_;
  my $dbh = $args{'-dbh'};

  my $sql = "SHOW GLOBAL variables LIKE 'max_allowed_packet'";
  my $sth = $dbh->prepare("$sql");
  $sth->execute();
  my $packet_size = $sth->fetchrow_array();
  unless ($packet_size > 5000000){
    print RED, "\nmax_allowed_packet = $packet_size bytes", RESET;
    print RED, "\nWARNING: max_allowed_packet variable is too small to allow insertion of large genes", RESET;
    print RED, "\nGet a DB Admin to execute the following: SET GLOBAL max_allowed_packet=10000000\n\n", RESET;
    exit();
  }

  print BLUE, "\nFound a suitable packet size ($packet_size) for the target ALEXA database\n\n", RESET;
  $sth->finish();
  return();
}


=head2 getAllGenes()

=over 3

=item Function:

Get a complete list of gene IDs from the database (with options)

=item Return:

Array of ALEXA gene ids

=item Args:

'-dbh' => database handle

'-gene_type' => 'All' or 'Non-pseudo' or 'protein_coding'

'-evidence' => 'Known Gene'  Use this to limit your pool of genes to only known genes

Other valid gene types are: miRNA, miRNA_pseudogene, misc_RNA, misc_RNA_pseudogene, Mt_rRNA, Mt_tRNA, Mt_tRNA_pseudogene, pseudogene, rRNA, rRNA_pseudogene, scRNA, scRNA_pseudogene, snoRNA, snoRNA_pseudogene, snRNA, snRNA_pseudogene, tRNA_pseudogene

=item Example(s):

my @gene_ids = @{&getAllGenes ('-dbh'=>$dbh, '-gene_type'=>'protein_coding', '-evidence'=>"Known Gene")};

=back

=cut

#######################
#getAllGenes          #
#######################
sub getAllGenes{
  my %args = @_;

  my $dbh = $args{'-dbh'};
  my $gene_type = $args{'-gene_type'};
  my $evidence = $args{'-evidence'};

  my @gene_ids;

  #Make sure the user provided a $dbh
  unless ($dbh){
    print "\nCannot get genes without a valid -dbh\n\n";
    exit();
  }
  if ($evidence){
    unless (($evidence eq "Known Gene") || ($evidence eq "Unknown Gene")){
      print "\nEvidence type: $evidence not understood by getAllGenes()\n\n";
      exit();
    }
  }

  if ($gene_type eq "All"){

    #Get all gene ids in the Gene table
    my $sql;
    if ($evidence){
      $sql = "SELECT id FROM Gene WHERE evidence = \"$evidence\"";
    }else{
      $sql = "SELECT id FROM Gene";
    }
    my $sth = $dbh->prepare("$sql");
    $sth->execute();

    while (my ($gene_id) = $sth->fetchrow_array()){
      push (@gene_ids, $gene_id);
    }
    $sth->finish();

    return(\@gene_ids);

  }elsif ($gene_type eq "Non-pseudo"){

    #Get all gene ids except those corresponding to pseudogenes in the Gene table (still includes miRNA, mtRNA, etc.)
    my $sql;
    if ($evidence){
      $sql = "SELECT id,gene_type FROM Gene WHERE evidence = \"$evidence\"";
    }else{
      $sql = "SELECT id,gene_type FROM Gene";
    }
    my $sth = $dbh->prepare("$sql");
    $sth->execute();

    while (my ($gene_id,$gt) = $sth->fetchrow_array()){

      unless ($gt =~ /pseudo/){
	push (@gene_ids, $gene_id);
      }
    }
    $sth->finish();

    return(\@gene_ids);

  }else{

    #Get only the gene ids corresponding to standard ensembl coding genes in the Gene table
    my $sql;
    if ($evidence){
      $sql = "SELECT id FROM Gene WHERE gene_type = \"$gene_type\" AND evidence = \"$evidence\"";
    }else{
      $sql = "SELECT id FROM Gene WHERE gene_type = \"$gene_type\"";
    }
    my $sth = $dbh->prepare("$sql");
    $sth->execute();

    while (my ($gene_id) = $sth->fetchrow_array()){
      push (@gene_ids, $gene_id);
    }
    $sth->finish();
  }
  return(\@gene_ids);
}

=head2 getGeneIds()

=over 3

=item Function:

Get the ALEXA gene id for a particular ensembl gene ID

=item Return:

Hash of gene ids (keyed on input gene ID!) - for those actually found in the database

=item Args:

'-dbh' => database handle

'-ensembl_g_id' => '$ensembl_g_id'

'-alexa_id' => '$alexa_id'

=item Example(s):

my %gene_ids = %{&getGeneIds ('-dbh'=>$dbh, '-ensembl_g_ids'=>\@ensembl_g_ids)};

my %gene_ids = %{&getGeneIds ('-dbh'=>$dbh, '-alexa_ids'=>\@alexa_ids)};

=back

=cut


#######################
#getGeneIds           #
#######################
sub getGeneIds{

  my %args = @_;
  my $dbh = $args{'-dbh'};

  my @ensembl_g_ids;
  my @alexa_ids;
  if ($args{'-ensembl_g_ids'}){
    @ensembl_g_ids = @{$args{'-ensembl_g_ids'}};
  }
  if ($args{'-alexa_ids'}){
    @alexa_ids = @{$args{'-alexa_ids'}};
  }

  #Make sure the user provided a $dbh
  unless ($dbh){
    print "\nCannot get gene ID without a valid -dbh\n\n";
    exit();
  }
  unless ($args{'-ensembl_g_ids'} || $args{'-alexa_ids'}){
    print "\n\nMust provide a valid list of ensembl or alexa gene ids for getGeneID()\n";
    exit();
  }

  my $input_gene_count;
  my %gene_ids;

  #Build the query depending on which IDs were provided
  if ($args{'-ensembl_g_ids'}){
    #List of EnsEMBL IDs - Get ALEXA IDs

    $input_gene_count = @ensembl_g_ids;
    unless ($input_gene_count > 0){
      print "\nList of EnsEMBL gene Ids does not seem valid! - getGeneIds()\n\n";
      exit();
    }
    my $id_string = join("','", @ensembl_g_ids);
    $id_string =~ s/^\,//;
    $id_string = "'"."$id_string"."'";

    #Get the ALEXA gene id corresponding to this ensembl gene id 
    my $sql = "SELECT id,ensembl_g_id FROM Gene WHERE ensembl_g_id IN ($id_string)";

    my $sth = $dbh->prepare("$sql");
    $sth->execute();

    while (my ($alexa_id, $ensembl_g_id) = $sth->fetchrow_array()){
      $gene_ids{$ensembl_g_id}{alexa_gene_id} = $alexa_id;
    }
    $sth->finish();

  }else{
    #List of ALEXA Ids - get EnsEMBL IDs

    $input_gene_count = @alexa_ids;
    unless ($input_gene_count > 0){
      print "\nList of ALEXA gene Ids does not seem valid! - getGeneIds()\n\n";
      exit();
    }
    my $id_string = join(",", @alexa_ids);
    $id_string =~ s/^\,//;

    #Get the ALEXA gene id corresponding to this ensembl gene id 
    my $sql = "SELECT id,ensembl_g_id FROM Gene WHERE id IN ($id_string)";
    my $sth = $dbh->prepare("$sql");
    $sth->execute();

    while (my ($alexa_id, $ensembl_g_id) = $sth->fetchrow_array()){
      $gene_ids{$alexa_id}{ensembl_g_id} = $ensembl_g_id;
    }
    $sth->finish();
  }

  my $genes_found = keys %gene_ids;
  unless ($genes_found > 0){
    print "\ngetGeneIds() Could not find a valid gene ID matching the input list of $input_gene_count genes\n\n";
    exit();
  }

  return(\%gene_ids)
}


=head2 getGeneInfo()

=over 3

=item Function:

Get a complete gene object for a particular gene ID

=item Return:

Gene object as a hash

=item Args:

'-dbh' => database handle

'-gene_id' => '$gene_id'

'-sequence' => 1 to include sequence

'-silent' => yes (prevent output to standard out)


=item Example(s):

my %gene = %{&getGeneInfo ('-dbh'=>$dbh, '-gene_id'=>$gene_id, '-sequence'=>"no")};

my %genes = %{&getGeneInfo ('-dbh'=>$dbh, '-gene_ids'=>\@gene_ids,'-sequence'=>"no")};

=back

=cut


#######################
#getGeneInfo          #
#######################
sub getGeneInfo{

  my %args = @_;
  my $dbh = $args{'-dbh'};
  my $gene_id = $args{'-gene_id'};
  my $sequence_flag = $args{'-sequence'};
  my $species = $args{'-species'};
  my $silent = $args{'-silent'};

  unless($silent){
    $silent = "no";
  }

  my @gene_ids;
  if ($args{'-gene_ids'}){
    unless ($silent eq "yes"){
      print "\nProcessing multiple Gene IDs supplied as a list\n";
    }
    @gene_ids = @{$args{'-gene_ids'}};
  }

  my $gene_count = @gene_ids;

  #Make sure the user provided a $dbh
  unless ($dbh){
    print "\nCannot get genes without a valid -dbh\n\n";
    exit();
  }
  unless ($gene_id || $gene_count > 0){
    print "\n\nMust provide a valid ALEXA gene id or array of gene IDs for getGeneInfo()\n";
    exit();
  }
  unless ($sequence_flag){
    print "\n\nMust specify whether you want sequence data or not for getGeneInfo()\n";
    exit();
  }
  unless ($sequence_flag eq "yes" || $sequence_flag eq "no"){
    print "\n\nSequence flag not understood by getGeneInfo()\n";
    exit();
  }

  #1.) For multiple Gene IDs provided as a list
  my $genes_processed = 0;
  my @sub_gene_list;
  if ($gene_count > 0){
    my %genes;

    #Process 100 genes at a time
    foreach my $gene (@gene_ids){
      $genes_processed++;
      push (@sub_gene_list, $gene);

      #If the sub gene list contains 100 genes or the last gene has been reached submit the query and reset variables
      my $current_gene_count = @sub_gene_list;

      if ($current_gene_count == 100 || $genes_processed == $gene_count){
	#Disable print buffer to get these print statements to display immediately instead of on return()
	unless ($silent eq "yes"){
	  $| = 1;
	  print ".";
	  $| = 0;
	}

	my $id_string = join(",", @sub_gene_list);
	$id_string =~ s/^\,//;

	#Submit query.
	my $sql_gene = "SELECT id,ensembl_g_id,ensembl_version,source,gene_type,gene_name,evidence,sequence,chr_strand,corrected_strand,chromosome,chr_start,chr_end,gene_start,gene_end FROM Gene WHERE id IN ($id_string)";

	my $sth_gene = $dbh->prepare("$sql_gene");
	$sth_gene->execute();

	while (my ($g_id,$ensembl_g_id,$ensembl_version,$source,$gene_type,$gene_name,$evidence,$sequence,$chr_strand,$corrected_strand,$chromosome,$chr_start,$chr_end,$gene_start,$gene_end) = $sth_gene->fetchrow_array()){

	  #Populate hash
	  $genes{$g_id}{ensembl_g_id} = $ensembl_g_id;
	  $genes{$g_id}{ensembl_version} = $ensembl_version;
	  $genes{$g_id}{source} = $source;
	  $genes{$g_id}{gene_type} = $gene_type;
	  $genes{$g_id}{gene_name} = $gene_name;
	  $genes{$g_id}{evidence} = $evidence;
	  $genes{$g_id}{chr_strand} = $chr_strand;
	  $genes{$g_id}{corrected_strand} = $corrected_strand;
	  $genes{$g_id}{chr_start} = $chr_start;
	  $genes{$g_id}{chr_end} = $chr_end;
	  $genes{$g_id}{gene_start} = $gene_start;
	  $genes{$g_id}{gene_end} = $gene_end;

	  if ($sequence_flag eq "yes"){
	    $genes{$g_id}{sequence} = $sequence;
	  }
	  #Determine sequence length and store for convenience
	  $genes{$g_id}{seq_length} = length($sequence);

	  ###############################################################################################
	  #CHROMOSOME NAME DB SPECIFIC CODE
	  #Provide a formatted chromosome name in case roman numbers, etc. were used
	  #Attempt to maximize compatibility between UCSC and EnsEMBL.
	  #UCSC seems to generally use numerals or single letters ...
	  #This will be hardcoded for particular databases
	  my $test_name = $dbh->{Name};
	  #database=ALEXA_sc_41_1d;host=....
	  my $database = '';
	  if ($test_name =~ /database\=(.*)\;/){
	    $database = $1;
	  }
	  if ($database eq "ALEXA_cf_42_2"){
	    if ($chromosome eq "MT"){$chromosome = "M";}
	  #}elsif($database eq "ALEXA_dr_42_6c"){
	  #  if ($chromosome =~ /^Zv6_NA\d+/){$chromosome = "NA_random";}
	  #Doesn't work.  Coordinates seem to be handles differently for EnsEMBL vs. UCSC for these
	  }elsif($database eq "ALEXA_gg_42_2"){
	    if ($chromosome eq "MT"){$chromosome = "M";}
	  }elsif($database eq "ALEXA_hs_41_36c"){
	    if ($chromosome eq "MT"){$chromosome = "M";}
	    if ($chromosome eq "c5_H2"){$chromosome = "5_h2_hap1";}
	    if ($chromosome eq "c6_COX"){$chromosome = "6_cox_hap1";}
	    if ($chromosome eq "c6_QBL"){$chromosome = "6_qbl_hap2";}
	  }elsif($database eq "ALEXA_mm_41_36b"){
	    if ($chromosome eq "MT"){$chromosome = "M";}
	  }elsif($database eq "ALEXA_pt_41_21"){
	    if ($chromosome eq "MT"){$chromosome = "M";}
	  }elsif($database eq "ALEXA_rn_43_34m"){
	    if ($chromosome eq "MT"){$chromosome = "M";}
	  }elsif ($database eq "ALEXA_sc_41_1d"){
	    if ($chromosome eq "Mito"){$chromosome = "M";}
	    if ($chromosome eq "I"){$chromosome = "1";}
	    if ($chromosome eq "II"){$chromosome = "2";}
	    if ($chromosome eq "III"){$chromosome = "3";}
	    if ($chromosome eq "IV"){$chromosome = "4";}
	    if ($chromosome eq "V"){$chromosome = "5";}
	    if ($chromosome eq "VI"){$chromosome = "6";}
	    if ($chromosome eq "VII"){$chromosome = "7";}
	    if ($chromosome eq "VIII"){$chromosome = "8";}
	    if ($chromosome eq "IX"){$chromosome = "9";}
	    if ($chromosome eq "X"){$chromosome = "10";}
	    if ($chromosome eq "XI"){$chromosome = "11";}
	    if ($chromosome eq "XII"){$chromosome = "12";}
	    if ($chromosome eq "XIII"){$chromosome = "13";}
	    if ($chromosome eq "XIV"){$chromosome = "14";}
	    if ($chromosome eq "XV"){$chromosome = "15";}
	    if ($chromosome eq "XVI"){$chromosome = "16";}
	  }

	  $genes{$g_id}{chromosome} = $chromosome;
	  ###################################################################################################

	}
	$sth_gene->finish();

	#Reset variables
	@sub_gene_list = ();
      }
    }
    return(\%genes);

  }else{
    #2.) For a single Gene ID
    my %gene;

    #Use the gene id to get the related gene info
    my $sql_gene = "SELECT ensembl_g_id,ensembl_version,source,gene_type,gene_name,evidence,sequence,chr_strand,corrected_strand,chromosome,chr_start,chr_end,gene_start,gene_end FROM Gene WHERE id = '$gene_id'";

    my $sth_gene = $dbh->prepare("$sql_gene");
    $sth_gene->execute();

    my ($ensembl_g_id,$ensembl_version,$source,$gene_type,$gene_name,$evidence,$sequence,$chr_strand,$corrected_strand,$chromosome,$chr_start,$chr_end,$gene_start,$gene_end) = $sth_gene->fetchrow_array();
    $sth_gene->finish();

    $gene{$gene_id}{ensembl_g_id} = $ensembl_g_id;
    $gene{$gene_id}{ensembl_version} = $ensembl_version;
    $gene{$gene_id}{source} = $source;
    $gene{$gene_id}{gene_type} = $gene_type;
    $gene{$gene_id}{gene_name} = $gene_name;
    $gene{$gene_id}{evidence} = $evidence;
    $gene{$gene_id}{sequence} = $sequence;
    $gene{$gene_id}{chr_strand} = $chr_strand;
    $gene{$gene_id}{corrected_strand} = $corrected_strand;
    $gene{$gene_id}{chr_start} = $chr_start;
    $gene{$gene_id}{chr_end} = $chr_end;
    $gene{$gene_id}{gene_start} = $gene_start;
    $gene{$gene_id}{gene_end} = $gene_end;

    #Provide a formatted chromosome name in case roman numbers, etc. were used
    #Attempt to maximize compatibility between UCSC and EnsEMBL.
    #UCSC seems to generally use numerals or single letters ...
    #This will be hardcoded for particular databases
    my $test_name = $dbh->{Name};
    #database=ALEXA_sc_41_1d;host=....
    my $database = '';
    if ($test_name =~ /database\=(.*)\;/){
      $database = $1;
    }
    if ($database eq "ALEXA_sc_41_1d"){
      if ($chromosome eq "Mito"){$chromosome = "M";}
      if ($chromosome eq "I"){$chromosome = "1";}
      if ($chromosome eq "II"){$chromosome = "2";}
      if ($chromosome eq "III"){$chromosome = "3";}
      if ($chromosome eq "IV"){$chromosome = "4";}
      if ($chromosome eq "V"){$chromosome = "5";}
      if ($chromosome eq "VI"){$chromosome = "6";}
      if ($chromosome eq "VII"){$chromosome = "7";}
      if ($chromosome eq "VIII"){$chromosome = "8";}
      if ($chromosome eq "IX"){$chromosome = "9";}
      if ($chromosome eq "X"){$chromosome = "10";}
      if ($chromosome eq "XI"){$chromosome = "11";}
      if ($chromosome eq "XII"){$chromosome = "12";}
      if ($chromosome eq "XIII"){$chromosome = "13";}
      if ($chromosome eq "XIV"){$chromosome = "14";}
      if ($chromosome eq "XV"){$chromosome = "15";}
      if ($chromosome eq "XVI"){$chromosome = "16";}
    }
    $gene{$gene_id}{chromosome} = $chromosome;

    return(\%gene);
  }
}


=head2 getTranscripts()

=over 3

=item Function:

Get the transcripts and corresponding exons for a particular gene ID

=item Return:

Genes object (keyed on genes ids) as a hash reference

Each gene has a reference to a transcript object which in turn contains a reference to an exon object

=item Args:

'-dbh' => database handle

'-gene_ids' => \@gene_ids

'-sequence' => "yes" or "no"

=item Example(s):

my $gene_transcripts_ref = &getTranscripts('-dbh'=>$dbh, '-gene_ids'=>\@gene_ids, '-sequence'=>"no");

=back

=cut


#######################
#getTranscripts       #
#######################
sub getTranscripts{

  my %args = @_;
  my $dbh = $args{'-dbh'};
  my @gene_ids = @{$args{'-gene_ids'}};
  my $sequence_flag = $args{'-sequence'};
  my $silent = $args{'-silent'};

  unless($silent){
    $silent = "no";
  }

  unless ($silent eq "yes"){
    print "\nProcessing multiple Gene IDs supplied as a list\n";
  }

  my $gene_count = @gene_ids;

  #Make sure the user provided a $dbh
  unless ($dbh){
    print "\nCannot get transcripts without a valid -dbh\n\n";
    exit();
  }
  unless ($gene_count > 0){
    print "\n\nMust provide at least one valid ALEXA gene id as an array for getTranscripts()\n";
    exit();
  }
  unless ($sequence_flag){
    print "\n\nMust specify whether you want sequence data or not for getTranscripts()\n";
    exit();
  }
  unless ($sequence_flag eq "yes" || $sequence_flag eq "no"){
    print "\n\nSequence flag not understood by getTranscripts()\n";
    exit();
  }


  #Multiple Gene IDs provided as a list
  my $genes_processed = 0;
  my @sub_gene_list;

  my %gene_transcripts;


  #Process 100 genes at a time
  foreach my $gene (@gene_ids){
    $genes_processed++;
    push (@sub_gene_list, $gene);

    #If the sub gene list contains 100 genes or the last gene has been reached submit the query and reset variables
    my $current_gene_count = @sub_gene_list;

    if ($current_gene_count == 100 || $genes_processed == $gene_count){

      #Disable print buffer to get these print statements to display immediately instead of on return()
      unless ($silent eq "yes"){
	$| = 1;
	print ".";
	$| = 0;
      }
      my $id_string = join(",", @sub_gene_list);
      $id_string =~ s/^\,//;

      #1.) Now use the gene ids to get all the transcripts associated with each gene
      my $sql_transcript;
      my $test = $dbh->{Name};

      if ($test =~ /ALEXA_cf_42_2|ALEXA_dm_42_43|ALEXA_dr_42_6c|ALEXA_gg_42_2|ALEXA_hs_41_36c|ALEXA_mm_41_36b|ALEXA_pt_41_21|ALEXA_rn_43_34m|ALEXA_sc_41_1d/){
	$sql_transcript = "SELECT fk_Gene__id, id, ensembl_t_id, start, end, coding_region_start, coding_region_end FROM Transcript WHERE fk_Gene__id IN ($id_string)";

	my $sth_transcript = $dbh->prepare("$sql_transcript");
	$sth_transcript->execute();

	while (my ($fk_gene_id,$transcript_id,$ensembl_t_id,$transcript_start,$transcript_end,$cds_start,$cds_end) = $sth_transcript->fetchrow_array()){

	  #If this gene has already been observed, simply add the new transcript record
	  if ($gene_transcripts{$fk_gene_id}{transcripts}){
	    my $transcripts_ref = $gene_transcripts{$fk_gene_id}{transcripts};
	    $transcripts_ref->{$transcript_id}->{ensembl_t_id} = $ensembl_t_id;
	    $transcripts_ref->{$transcript_id}->{transcript_start} = $transcript_start;
	    $transcripts_ref->{$transcript_id}->{transcript_end} = $transcript_end;
	    $transcripts_ref->{$transcript_id}->{cds_start} = $cds_start;
	    $transcripts_ref->{$transcript_id}->{cds_end} = $cds_end;

	  }else{
	    #If this is the first transcript of this gene, create a new transcript object
	    my %transcripts;
	    $transcripts{$transcript_id}{ensembl_t_id} = $ensembl_t_id;
	    $transcripts{$transcript_id}{transcript_start} = $transcript_start;
	    $transcripts{$transcript_id}{transcript_end} = $transcript_end;
	    $transcripts{$transcript_id}{cds_start} = $cds_start;
	    $transcripts{$transcript_id}{cds_end} = $cds_end;
	    $gene_transcripts{$fk_gene_id}{transcripts} = \%transcripts;
	  }
	}
	$sth_transcript->finish();

      }else{
	$sql_transcript = "SELECT fk_Gene__id, id, ensembl_t_id, start, end, coding_region_start, coding_region_end, cds_start_coords, cds_end_coords FROM Transcript WHERE fk_Gene__id IN ($id_string)";

	my $sth_transcript = $dbh->prepare("$sql_transcript");
	$sth_transcript->execute();

	while (my ($fk_gene_id,$transcript_id,$ensembl_t_id,$transcript_start,$transcript_end,$cds_start,$cds_end,$cds_start_coords,$cds_end_coords) = $sth_transcript->fetchrow_array()){

	  my @cds_start_coords = split(" ", $cds_start_coords);
	  my @cds_end_coords = split(" ", $cds_end_coords);

	  #If this gene has already been observed, simply add the new transcript record
	  if ($gene_transcripts{$fk_gene_id}{transcripts}){
	    my $transcripts_ref = $gene_transcripts{$fk_gene_id}{transcripts};
	    $transcripts_ref->{$transcript_id}->{ensembl_t_id} = $ensembl_t_id;
	    $transcripts_ref->{$transcript_id}->{transcript_start} = $transcript_start;
	    $transcripts_ref->{$transcript_id}->{transcript_end} = $transcript_end;
	    $transcripts_ref->{$transcript_id}->{cds_start} = $cds_start;
	    $transcripts_ref->{$transcript_id}->{cds_end} = $cds_end;
	    $transcripts_ref->{$transcript_id}->{cds_start_coords} = \@cds_start_coords;
	    $transcripts_ref->{$transcript_id}->{cds_end_coords} = \@cds_end_coords;

	    if ($cds_start =~ /\d+/ && $cds_end =~ /\d+/){
	      if ($cds_start < $gene_transcripts{$fk_gene_id}{grand_cds_start}){
		$gene_transcripts{$fk_gene_id}{grand_cds_start} = $cds_start;
	      }
	      if ($cds_end > $gene_transcripts{$fk_gene_id}{grand_cds_end}){
		$gene_transcripts{$fk_gene_id}{grand_cds_end} = $cds_end;
	      }
	    }

	  }else{
	    #If this is the first transcript of this gene, create a new transcript object
	    my %transcripts;
	    $transcripts{$transcript_id}{ensembl_t_id} = $ensembl_t_id;
	    $transcripts{$transcript_id}{transcript_start} = $transcript_start;
	    $transcripts{$transcript_id}{transcript_end} = $transcript_end;
	    $transcripts{$transcript_id}{cds_start} = $cds_start;
	    $transcripts{$transcript_id}{cds_end} = $cds_end;
	    $transcripts{$transcript_id}{cds_start_coords} = \@cds_start_coords;
	    $transcripts{$transcript_id}{cds_end_coords} = \@cds_end_coords;

	    #Keep track of the 5' most CDS start coord and 3' most CDS end coord (in the case of multiple transcripts)
	    $gene_transcripts{$fk_gene_id}{grand_cds_start} = $cds_start;
	    $gene_transcripts{$fk_gene_id}{grand_cds_end} = $cds_end;

	    $gene_transcripts{$fk_gene_id}{transcripts} = \%transcripts;
	  }
	}
	$sth_transcript->finish();

      }


      #2.) Next get all the exons associated with each transcript
      my $sql_exon = "SELECT Transcript.fk_Gene__id,Transcript.id,Exon.id,ensembl_e_id,Exon.start,Exon.end FROM Exon,TranscriptExon,Transcript WHERE Transcript.id=TranscriptExon.fk_Transcript__id AND TranscriptExon.fk_Exon__id=Exon.id AND Transcript.fk_Gene__id IN ($id_string)";
      #print "\n\nDEBUG $sql_exon\n\n";

      my $sth_exon = $dbh->prepare("$sql_exon");

      $sth_exon->execute();

      while (my ($fk_gene_id,$transcript_id,$exon_id, $ensembl_e_id, $exon_start, $exon_end) = $sth_exon->fetchrow_array()){

	#The gene and transcript objects should have already been created above, simply need to create the exon object for each transcript
	my $transcripts_ref = $gene_transcripts{$fk_gene_id}{transcripts};

	#If this transcript has already been seen before, simply add the exon record
	if ($transcripts_ref->{$transcript_id}->{exons}){
	  my $exons_ref = $transcripts_ref->{$transcript_id}->{exons};
	  $exons_ref->{$exon_id}->{ensembl_e_id} = $ensembl_e_id;
	  $exons_ref->{$exon_id}->{exon_start} = $exon_start;
	  $exons_ref->{$exon_id}->{exon_end} = $exon_end;
	}else{
	  #If this is the first exon of this transcript, create a new exon record
	  my %exons;
	  $exons{$exon_id}{ensembl_e_id} = $ensembl_e_id;
	  $exons{$exon_id}{exon_start} = $exon_start;
	  $exons{$exon_id}{exon_end} = $exon_end;
	  $transcripts_ref->{$transcript_id}->{exons} = \%exons;
	}

      }
      $sth_exon->finish();

      #3.) Get the sequence for each exon if the user specified to do so.
      if ($sequence_flag eq "yes"){
	
	#Get the complete gene sequence for this gene.
	my $sql_seq = "SELECT id,sequence from Gene WHERE id IN ($id_string)";

	my $sth_seq = $dbh->prepare("$sql_seq");
	$sth_seq->execute();

	while (my ($gene_id,$sequence) = $sth_seq->fetchrow_array()){

	  #Get the exon sequence for each exon, using the start and stop positions
	  my $transcripts_ref = $gene_transcripts{$gene_id}{transcripts};

	  foreach my $trans_id (keys %{$transcripts_ref}){

	    my $exons_ref = $transcripts_ref->{$trans_id}->{exons};

	    foreach my $exon_id (keys %{$exons_ref}){
	      my $start = $exons_ref->{$exon_id}->{exon_start};
	      my $end = $exons_ref->{$exon_id}->{exon_end};

	      my $length = $end-$start;

	      my $exon_seq = substr($sequence,$start-1,$length+1);

	      $exons_ref->{$exon_id}->{sequence} = $exon_seq;
	    }
	  }
	}
	$sth_seq->finish();
      }
      #Reset variables
      @sub_gene_list = ();
    }
  }
  return(\%gene_transcripts);
}


=head2 getProbeInfo()

=over 3

=item Function:

Get information for a specific ALEXA probe ID

Note that each probe has two IDs:

One which is assured to match that from generated probe files

The second is an auto-increment ID, which will only match if the probe ID in probe files if ALL probes were imported in order

If in doubt, use the '-probe_count_id' option

=item Return:

Probe object (keyed on various properties) as a hash

=item Args:

'-dbh' => database handle

'-probe_id' => '$probe_id'

'-probe_count_id' => '$probe_id'

=item Example(s):

my %probe = %{&getProbeInfo ('-dbh'=>$dbh, '-probe_id'=>$probe_id)};

my %probe = %{&getProbeInfo ('-dbh'=>$dbh, '-probe_count_id'=>$probe_id)};

$probe_pair_id = $probe{probe_pair_id};

Other possible values: sequence,type,offset,tm,exons_skipped,gene_id,unit1_start,unit1_end,unit2_start,unit2_end

Note that for negative control probes, all values related to the gene target will be 'na'

=back

=cut


#######################
#getProbeInfo         #
#######################
sub getProbeInfo{

  my %args = @_;
  my $dbh = $args{'-dbh'};
  my $probe_id = $args{'-probe_id'};
  my $probe_count_id = $args{'-probe_count_id'};

  #Make sure the user provided a $dbh
  unless ($dbh){
    print "\nCannot get probe info without a valid -dbh\n\n";
    exit();
  }
  unless ($probe_id || $probe_count_id){
    print "\n\nMust provide a valid probe_id or probe_count_id for getProbeInfo()\n";
    exit();
  }

  unless($probe_id){
    $probe_id = 'na'
  }
  unless($probe_count_id){
    $probe_count_id = 'na';
  }

  my $probe_table = "Probe";
  my $geneprobe_table = "GeneProbe";

  my %probe;

  #1.) First get info from the 'Probe' table
  my $sql1;

  #Build query for each case depending on the type of probe ID provided
  if ($probe_count_id =~ /^\d+$/){
    $sql1 = "SELECT id,probe_set_id,sequence,probe_length,type,tm,exons_skipped FROM $probe_table WHERE probe_count = '$probe_count_id'";
  }elsif($probe_id =~ /^\d+$/){
    $sql1 = "SELECT id,probe_set_id,sequence,probe_length,type,tm,exons_skipped FROM $probe_table WHERE id = '$probe_id'";
  }

  my $sth1 = $dbh->prepare("$sql1");
  $sth1->execute();

  my ($id,$probe_set_id,$sequence,$probe_length,$type,$tm,$exons_skipped) = $sth1->fetchrow_array();

  #Build probe hash object
  $probe{db_id} = $id;
  $probe{probe_set_id} = $probe_set_id;
  $probe{sequence} = $sequence;
  $probe{probe_length} = $probe_length;
  $probe{type} = $type;
  $probe{tm} = $tm;
  $probe{exons_skipped} = $exons_skipped;
  $sth1->finish();

  #2.) Get info from the 'GeneProbe' table - Note that negative control probe do not have a target ID and therefore these values are 'na'
  my $sql2;

  if ($probe_count_id =~ /^\d+$/){
    $sql2 = "SELECT $geneprobe_table.fk_Gene__id,$geneprobe_table.unit1_start,$geneprobe_table.unit1_end,$geneprobe_table.unit2_start,$geneprobe_table.
unit2_end FROM $probe_table,$geneprobe_table WHERE $probe_table.id=$geneprobe_table.fk_Probe__id AND $probe_table.probe_count = '$probe_count_id'";
  }elsif($probe_id =~ /^\d+$/){
    $sql2 = "SELECT $geneprobe_table.fk_Gene__id,$geneprobe_table.unit1_start,$geneprobe_table.unit1_end,$geneprobe_table.unit2_start,$geneprobe_table.
unit2_end FROM $probe_table,$geneprobe_table WHERE $probe_table.id=$geneprobe_table.fk_Probe__id AND $probe_table.id = '$probe_id'";
  }

  my $sth2 = $dbh->prepare("$sql2");
  $sth2->execute();

  my ($gene_id,$unit1_start,$unit1_end,$unit2_start,$unit2_end) = $sth2->fetchrow_array();

  if ($gene_id){
    $probe{gene_id} = $gene_id;
    $probe{unit1_start} = $unit1_start;
    $probe{unit1_end} = $unit1_end;
    $probe{unit2_start} = $unit2_start;
    $probe{unit2_end} = $unit2_end;
  }else{
    $probe{gene_id} = 'na';
    $probe{unit1_start} = 'na';
    $probe{unit1_end} = 'na';
    $probe{unit2_start} = 'na';
    $probe{unit2_end} = 'na';
  }

  $sth2->finish();

  return(\%probe);
}


=head2 getGeneProbes()

=over 3

=item Function:

Get information for a specific ALEXA probe ID

Note that each probe has two IDs:

One which is assured to match that from generated probe files

The second is an auto-increment ID, which will only match if the probe ID in probe files if ALL probes were imported in order

If in doubt, use the '-probe_count_id' option

=item Return:

Probe object (keyed on various properties) as a hash

=item Args:

'-dbh' => database handle

'-gene_ids' => '\@gene_ids'

'-silent' => 'yes' or 'no'

=item Example(s):

my %gene_probes = %{&getGeneProbes ('-dbh'=>$dbh, '-gene_ids'=>\@gene_ids, '-filtered'=>"no", '-silent'=>"yes")};

my %gene_probes = %{&getGeneProbes ('-dbh'=>$dbh, '-gene_ids'=>\@gene_ids, '-filtered'=>"yes", '-filter_set'=>1, '-silent'=>"yes")};

my $probes_ref = $gene_probes{$gene_id}{probes};

Possible values for a particular probe ID: type,offset,tm,exons_skipped,gene_id,unit1_start,unit1_end,unit2_start,unit2_end,mrna_thl,est_thl,$enst_thl

Note that for negative control probes, all values related to the gene target will be 'na'

=back

=cut


#######################
#getGeneProbes         #
#######################
sub getGeneProbes{

  my %args = @_;
  my $dbh = $args{'-dbh'};
  my @gene_ids = @{$args{'-gene_ids'}};
  my $filtered = $args{'-filtered'};
  my $filter_set = $args{'-filter_set'};
  my $silent = $args{'-silent'};

  unless($silent){
    $silent = "no";
  }
  unless ($silent eq "yes"){
    print "\nProcessing multiple Gene IDs supplied as a list\n";
  }

  my $gene_count = @gene_ids;

  #Make sure the user provided a $dbh
  unless ($dbh){
    print "\nCannot get gene probe info without a valid -dbh\n\n";
    exit();
  }
  unless ($gene_count > 0){
    print "\n\nMust provide at least one valid ALEXA gene id as an array for getGeneProbes()\n";
    exit();
  }
  unless ($filtered){
    print "\n\nMust specify whether you want filtered probe only or not for getGeneProbes()\n";
    exit();
  }
  unless ($filtered eq "yes" || $filtered eq "no"){
    print "\n\nFiltered flag not understood by getGeneProbes()\n";
    exit();
  }
  if ($filtered eq "yes"){
    unless ($filter_set =~ /\d+/){
      print "\n\nIf you want filtered probes only you must specify a filter set for getGeneProbes()\n";
      exit();
    }
  }

  #Process multiple Gene IDs provided as a list
  my $genes_processed = 0;
  my @sub_gene_list;

  my %gene_probes;
  $gene_probes{total_probe_count} = 0;

  #Process 100 genes at a time
  foreach my $gene_id (@gene_ids){
    $genes_processed++;
    push (@sub_gene_list, $gene_id);

    #Initialize object
    $gene_probes{$gene_id}{probe_count} = 0;


    #If the sub gene list contains 100 genes or the last gene has been reached submit the query and reset variables
    my $current_gene_count = @sub_gene_list;

    if ($current_gene_count == 100 || $genes_processed == $gene_count){

      #Disable print buffer to get these print statements to display immediately instead of on return()
      unless ($silent eq "yes"){
	$| = 1;
	print ".";
	$| = 0;
      }

      #Build a string containing the list of genes in this batch
      my $id_string = join(",", @sub_gene_list);
      $id_string =~ s/^\,//;

      #Generate the SQL query for this batch of genes
      my $sql_probes;
      if ($filtered eq "no"){
	#All probes
	$sql_probes = "SELECT fk_Gene__id,probe_set_id,probe_count,type,sequence,probe_length,tm,unit1_start,unit1_end,unit1_start_chr,unit1_end_chr,unit2_start,unit2_end,unit2_start_chr,unit2_end_chr,mrna_target_hit_length,est_target_hit_length,enst_target_hit_length,exons_skipped FROM Probe,GeneProbe WHERE GeneProbe.fk_Probe__id=Probe.id AND GeneProbe.fk_Gene__id IN ($id_string)";

      }else{
	#Filtered probes only
	$sql_probes = "SELECT fk_Gene__id,probe_set_id,probe_count,type,sequence,probe_length,tm,unit1_start,unit1_end,unit1_start_chr,unit1_end_chr,unit2_start,unit2_end,unit2_start_chr,unit2_end_chr,mrna_target_hit_length,est_target_hit_length,enst_target_hit_length,exons_skipped FROM Probe,GeneProbe,ProbeProbe_set WHERE GeneProbe.fk_Probe__id=Probe.id AND Probe.id=ProbeProbe_set.fk_Probe__id AND ProbeProbe_set.fk_Probe_set__id = $filter_set AND GeneProbe.fk_Gene__id IN ($id_string)";
      }

      #Submit the query and add the results to the gene_probes object, update the number of probes found for each gene
      my $sth_probes = $dbh->prepare("$sql_probes");
      $sth_probes->execute();

      while (my ($gene_id,$probe_set_id,$probe_count,$probe_type,$sequence,$probe_length,$tm_celsius,$unit1_start,$unit1_end,$unit1_start_chr,$unit1_end_chr,$unit2_start,$unit2_end,$unit2_start_chr,$unit2_end_chr,$mrna_thl,$est_thl,$enst_thl,$exons_skipped) = $sth_probes->fetchrow_array()){

	unless ($mrna_thl){$mrna_thl="na";}
	unless ($est_thl){$est_thl="na";}
	unless ($enst_thl){$enst_thl="na";}

	$gene_probes{total_probe_count}++;
	$gene_probes{$gene_id}{probe_count}++;

	if ($gene_probes{$gene_id}{probesets}){
	  #Gene has already been observed
	  my $probesets_ref = $gene_probes{$gene_id}{probesets};
	
	  if ($probesets_ref->{$probe_set_id}){
	    #Probeset has also already been observed
	    my $probes_ref = $probesets_ref->{$probe_set_id}->{probes};
	    $probes_ref->{$probe_count}->{probe_type} = $probe_type;
	    $probes_ref->{$probe_count}->{sequence} = $sequence;
	    $probes_ref->{$probe_count}->{probe_length} = $probe_length;
	    $probes_ref->{$probe_count}->{tm_celsius} = $tm_celsius;
	    $probes_ref->{$probe_count}->{unit1_start} = $unit1_start;
	    $probes_ref->{$probe_count}->{unit1_end} = $unit1_end;
	    $probes_ref->{$probe_count}->{unit2_start} = $unit2_start;
	    $probes_ref->{$probe_count}->{unit2_end} = $unit2_end;
	    $probes_ref->{$probe_count}->{unit1_start_chr} = $unit1_start_chr;
	    $probes_ref->{$probe_count}->{unit1_end_chr} = $unit1_end_chr;
	    $probes_ref->{$probe_count}->{unit2_start_chr} = $unit2_start_chr;
	    $probes_ref->{$probe_count}->{unit2_end_chr} = $unit2_end_chr;
	    $probes_ref->{$probe_count}->{mrna_thl} = $mrna_thl;
	    $probes_ref->{$probe_count}->{est_thl} = $est_thl;
	    $probes_ref->{$probe_count}->{enst_thl} = $enst_thl;
	    $probes_ref->{$probe_count}->{exons_skipped} = $exons_skipped;
	    $probes_ref->{$probe_count}->{status} = "Pass";
	    $probes_ref->{$probe_count}->{sequence_support} = "Unknown";

	  }else{
	    #first probe for this probeset
	    my %probes;
	    $probes{$probe_count}{probe_type} = $probe_type;
	    $probes{$probe_count}{sequence} = $sequence;
	    $probes{$probe_count}{probe_length} = $probe_length;
	    $probes{$probe_count}{tm_celsius} = $tm_celsius;
	    $probes{$probe_count}{unit1_start} = $unit1_start;
	    $probes{$probe_count}{unit1_end} = $unit1_end;
	    $probes{$probe_count}{unit2_start} = $unit2_start;
	    $probes{$probe_count}{unit2_end} = $unit2_end;
	    $probes{$probe_count}{unit1_start_chr} = $unit1_start_chr;
	    $probes{$probe_count}{unit1_end_chr} = $unit1_end_chr;
	    $probes{$probe_count}{unit2_start_chr} = $unit2_start_chr;
	    $probes{$probe_count}{unit2_end_chr} = $unit2_end_chr;
	    $probes{$probe_count}{mrna_thl} = $mrna_thl;
	    $probes{$probe_count}{est_thl} = $est_thl;
	    $probes{$probe_count}{enst_thl} = $enst_thl;
	    $probes{$probe_count}{exons_skipped} = $exons_skipped;
	    $probes{$probe_count}{status} = "Pass";
	    $probes{$probe_count}{sequence_support} = "Unknown";
	    $probesets_ref->{$probe_set_id}->{probe_type} = $probe_type;
	    $probesets_ref->{$probe_set_id}->{probes} = \%probes;
	  }
	}else{
	  #First probe for this gene
	  my %probesets;
	  my %probes;
	  $probes{$probe_count}{probe_type} = $probe_type;
	  $probes{$probe_count}{sequence} = $sequence;
	  $probes{$probe_count}{probe_length} = $probe_length;
	  $probes{$probe_count}{tm_celsius} = $tm_celsius;
	  $probes{$probe_count}{unit1_start} = $unit1_start;
	  $probes{$probe_count}{unit1_end} = $unit1_end;
	  $probes{$probe_count}{unit2_start} = $unit2_start;
	  $probes{$probe_count}{unit2_end} = $unit2_end;
	  $probes{$probe_count}{unit1_start_chr} = $unit1_start_chr;
	  $probes{$probe_count}{unit1_end_chr} = $unit1_end_chr;
	  $probes{$probe_count}{unit2_start_chr} = $unit2_start_chr;
	  $probes{$probe_count}{unit2_end_chr} = $unit2_end_chr;
	  $probes{$probe_count}{mrna_thl} = $mrna_thl;
	  $probes{$probe_count}{est_thl} = $est_thl;
	  $probes{$probe_count}{enst_thl} = $enst_thl;
	  $probes{$probe_count}{exons_skipped} = $exons_skipped;
	  $probes{$probe_count}{status} = "Pass";
	  $probes{$probe_count}{sequence_support} = "Unknown";
	  $probesets{$probe_set_id}{probe_type} = $probe_type;
	  $probesets{$probe_set_id}{probes} = \%probes;
	  $gene_probes{$gene_id}{probesets} = \%probesets;
	}
      }

      $sth_probes->finish();

      #Reset variables
      @sub_gene_list = ();
    }
  }
  return(\%gene_probes);
}


=head2 getNcProbes()

=over 3

=item Function:

Get all filtered/unfiltered Negative Control probe from an ALEXA DB

Note that each probe has two IDs:

One which is assured to match that from generated probe files

The second is an auto-increment ID, which will only match if the probe ID in probe files if ALL probes were imported in order

If in doubt, use the '-probe_count_id' option

=item Return:

Probe object (keyed on various properties) as a hash

=item Args:

'-dbh' => database handle

'-filtered' => yes/no flag

=item Example(s):

my %probes = %{&getNcProbes ('-dbh'=>$dbh, '-filtered'=>"yes", '-filtered_set'=>$filtered_set_id)};

my %probes = %{&getNcProbes ('-dbh'=>$dbh, '-filtered'=>"no")};

Note that for negative control probes, many probe info values will be 'na'

=back

=cut


#######################
#getNc Probes         #
#######################
sub getNcProbes{

  my %args = @_;
  my $dbh = $args{'-dbh'};
  my $filtered = $args{'-filtered'};
  my $filtered_set = $args{'-filtered_set'};

  #Make sure the user provided a $dbh
  unless ($dbh){
    print "\nCannot get probe info without a valid -dbh\n\n";
    exit();
  }
  unless ($filtered =~ /^yes$|^no$/i){
    print "\n\nMust specify whether to retrieve filtered or unfiltered probes for getNcProbes()\n\n";
    exit();
  }

  if ($filtered =~ /^yes$/i){
    unless ($filtered_set =~ /\d+/){
      print "\n\nIf you want filtered probe you must specify a valid filtered set ID for getNcProbes()\n\n";
      exit();
    }
  }

  my %probes;
  my $sql;
  my $type = "Control-Negative";

  if ($filtered =~ /^yes$/i){
    $sql = "SELECT probe_count,probe_set_id,sequence,probe_length,tm FROM Probe,ProbeProbe_set WHERE type='$type' AND ProbeProbe_set.fk_Probe__id=Probe.id AND ProbeProbe_set.fk_Probe_set__id='$filtered_set'";
  }else{
    $sql = "SELECT probe_count,probe_set_id,sequence,probe_length,tm FROM Probe WHERE type='$type'";
  }

  my $sth = $dbh->prepare("$sql");
  $sth->execute();

  while (my ($probe_count,$probe_set_id,$sequence,$probe_length,$tm) = $sth->fetchrow_array()){
    $probes{$probe_count}{probe_set_id} = $probe_set_id;
    $probes{$probe_count}{sequence} = $sequence;
    $probes{$probe_count}{probe_length} = $probe_length;
    $probes{$probe_count}{tm_celsius} = $tm;
    $probes{$probe_count}{status} = "Pass";
  }
  $sth->finish();

  return(\%probes);
}


=head2 getExons()

=over 3

=item Function:

Get the exons for a particular gene ID (all exons from all transcripts)

NOTE: THESE ARE NOT NECCESSARILY UNIQUE!!!  If there are multiple transcripts you will often get duplicate exons with the same coordinates!

=item Return:

Exons object (keyed on exon ids or gene_id if multiple genes IDs were supplied) as a hash

=item Args:

'-dbh' => database handle

'-gene_ids' => \@gene_ids

=item Example(s):

my $gene_exons_ref = &getExons ('-dbh'=>$dbh, '-gene_ids'=>\@gene_ids, '-sequence'=>"no");

=back

=cut


#######################
#getExons             #
#######################
sub getExons{

  my %args = @_;
  my $dbh = $args{'-dbh'};
  my @gene_ids = @{$args{'-gene_ids'}};
  my $sequence_flag = $args{'-sequence'};
  my $silent = $args{'-silent'};

  unless($silent){
    $silent = "no";
  }
  unless ($silent eq "yes"){
    print "\nProcessing multiple Gene IDs supplied as a list\n";
  }

  my $gene_count = @gene_ids;

  #Make sure the user provided a $dbh
  unless ($dbh){
    print "\nCannot get genes without a valid -dbh\n\n";
    exit();
  }
  unless ($gene_count > 0){
    print "\n\nMust provide a valid ALEXA gene id or array of gene IDs for getExons()\n";
    exit();
  }
  unless ($sequence_flag){
    print "\n\nMust specify whether you want sequence data or not for getExons()\n";
    exit();
  }
  unless ($sequence_flag eq "yes" || $sequence_flag eq "no"){
    print "\n\nSequence flag not understood by getExons()\n";
    exit();
  }

  #Process multiple Gene IDs provided as a list
  my $genes_processed = 0;
  my @sub_gene_list;

  my %gene_exons;

  #Process 100 genes at a time
  foreach my $gene (@gene_ids){
    $genes_processed++;
    push (@sub_gene_list, $gene);

    #If the sub gene list contains 100 genes or the last gene has been reached submit the query and reset variables
    my $current_gene_count = @sub_gene_list;

    if ($current_gene_count == 100 || $genes_processed == $gene_count){

      #Disable print buffer to get these print statements to display immediately instead of on return()
      unless ($silent eq "yes"){
	$| = 1;
	print ".";
	$| = 0;
      }

      my $id_string = join(",", @sub_gene_list);
      $id_string =~ s/^\,//;

      #Submit query.
      my $sql_gene = "SELECT Transcript.fk_Gene__id,Exon.id,ensembl_e_id,Exon.start,Exon.end FROM Exon,TranscriptExon,Transcript WHERE Transcript.id=TranscriptExon.fk_Transcript__id AND TranscriptExon.fk_Exon__id=Exon.id AND Transcript.fk_Gene__id IN ($id_string)";

      my $sth_gene = $dbh->prepare("$sql_gene");
      $sth_gene->execute();

      while (my ($gene_id,$exon_id,$ensembl_e_id,$exon_start,$exon_end) = $sth_gene->fetchrow_array()){

	#Populate hash
	if ($gene_exons{$gene_id}{exons}){
	  my $exons_ref = $gene_exons{$gene_id}{exons};
	  $exons_ref->{$exon_id}->{ensembl_e_id} = $ensembl_e_id;
	  $exons_ref->{$exon_id}->{exon_start} = $exon_start;
	  $exons_ref->{$exon_id}->{exon_end} = $exon_end;
	}else{
	  my %exons;
	  $exons{$exon_id}{ensembl_e_id} = $ensembl_e_id;
	  $exons{$exon_id}{exon_start} = $exon_start;
	  $exons{$exon_id}{exon_end} = $exon_end;
	  $gene_exons{$gene_id}{exons} = \%exons;
	}
      }
      $sth_gene->finish();

      if ($sequence_flag eq "yes"){
	
	#Get the complete gene sequence for this gene.
	my $sql_seq = "SELECT id,sequence from Gene WHERE id IN ($id_string)";

	my $sth_seq = $dbh->prepare("$sql_seq");
	$sth_seq->execute();

	while (my ($gene_id,$sequence) = $sth_seq->fetchrow_array()){

	  #Get the exon sequence for each exon, using the start and stop positions
	  my $exons_ref = $gene_exons{$gene_id}{exons};

	  foreach my $exon_id (keys %{$exons_ref}){
	    my $start = $exons_ref->{$exon_id}->{exon_start};
	    my $end = $exons_ref->{$exon_id}->{exon_end};

	    my $length = $end-$start;

	    my $exon_seq = substr($sequence,$start-1,$length+1);

	    $exons_ref->{$exon_id}->{sequence} = $exon_seq;
	  }
	}
	$sth_seq->finish();
      }
      #Reset variables
      @sub_gene_list = ();
    }
  }

  return(\%gene_exons);
}


=head2 getMaskedExons()

=over 3

=item Function:

Get the UNIQUE masked exons for a particular gene ID (all exons from all transcripts)

Sequence is derived from a repeat masked meta-transcript - Sequences corresponding to repetitive elements are 

replaced with Ns

=item Return:

Exons object (keyed on exon ids) as a hash

=item Args:

'-dbh' => database handle

'-gene_id' => $gene_id

=item Example(s):

my %maskedExons = %{&getMaskedExons ('-dbh'=>$dbh, '-gene_id'=>$gene_id)};

=back

=cut


#######################
#getMasked Exons      #
#######################
sub getMaskedExons{

  my %args = @_;
  my $dbh = $args{'-dbh'};
  my $gene_id = $args{'-gene_id'};

  my %maskedExons;

  #Make sure the user provided a $dbh
  unless ($dbh){
    print "\nCannot get genes without a valid -dbh\n\n";
    exit();
  }
  unless ($gene_id){
    print "\n\nMust provide a valid ALEXA gene id for getMaskedExons()\n";
    exit();
  }

  #Get only the gene ids corresponding to standard ensembl coding genes in the Gene table
  my $sql = "SELECT Exon.id, Exon.ensembl_e_id FROM Exon,TranscriptExon,Transcript WHERE Transcript.fk_Gene__id = '$gene_id' AND Transcript.id=TranscriptExon.fk_Transcript__id AND TranscriptExon.fk_Exon__id=Exon.id";
  my $sth = $dbh->prepare("$sql");
  $sth->execute();

  while (my ($exon_id, $ensembl_e_id) = $sth->fetchrow_array()){
    $maskedExons{$ensembl_e_id}{exon_id} = $exon_id;
  }
  $sth->finish();

  #Get the masked sequence for this exon from the Exon table for each exon
  foreach my $ensembl_e_id (sort keys %maskedExons){
    my $exon_id = $maskedExons{$ensembl_e_id}{exon_id};
    my $sql_seq = "SELECT masked_sequence from Exon WHERE id = '$exon_id'";

    my $sth_seq = $dbh->prepare("$sql_seq");
    $sth_seq->execute();

    my ($sequence) = $sth_seq->fetchrow_array();

    $sth_seq->finish();
    $maskedExons{$ensembl_e_id}{masked_sequence} = $sequence;
  }
  return(\%maskedExons);
}


=head2 getMaskedGene()

=over 3

=item Function:

Get the complete masked gene sequence for a particular gene ID

Complete gene sequences were masked by RepeatMasker and stored in a seperate table 'MaskedGene' for performance reasons

Bases that correspond to repetitive elements as defined by a library of repeats are replaced with Ns

=item Return:

Gene object (keyed on gene id) as a hash

=item Args:

'-dbh' => database handle

'-gene_id' => \@gene_ids

=item Example(s):

my %maskedGene = %{&getMaskedGene ('-dbh'=>$dbh, '-gene_ids'=>\@gene_ids)};

=back

=cut


#######################
#getMaskedGene        #
#######################
sub getMaskedGene{

  my %args = @_;
  my $dbh = $args{'-dbh'};
  my @gene_ids = @{$args{'-gene_ids'}};

  print "\nProcessing multiple Gene IDs supplied as a list\n";

  my $gene_count = @gene_ids;

  #Make sure the user provided a $dbh
  unless ($dbh){
    print "\nCannot get genes without a valid -dbh\n\n";
    exit();
  }
  unless ($gene_count > 0){
    print "\n\nMust provide a valid ALEXA gene id or array of gene IDs for getMaskedGene()\n";
    exit();
  }

  my $genes_processed = 0;
  my @sub_gene_list;

  my %maskedGene;

  foreach my $gene (@gene_ids){
    $genes_processed++;
    push (@sub_gene_list, $gene);

    #If the sub gene list contains 100 genes or the last gene has been reached submit the query and reset variables
    my $current_gene_count = @sub_gene_list;

    if ($current_gene_count == 100 || $genes_processed == $gene_count){

      #Disable print buffer to get these print statements to display immediately instead of on return()
      $| = 1;
      print ".";
      $| = 0;

      my $id_string = join(",", @sub_gene_list);
      $id_string =~ s/^\,//;

      #Submit query.
      my $sql_gene = "SELECT id, fk_Gene__id, sequence FROM MaskedGene WHERE fk_Gene__id IN ($id_string)";

      my $sth_gene = $dbh->prepare("$sql_gene");
      $sth_gene->execute();

      while (my ($id,$gene_id,$masked_sequence) = $sth_gene->fetchrow_array()){
	$maskedGene{$gene_id}{masked_gene_id} = $id;
	$maskedGene{$gene_id}{sequence} = $masked_sequence;
      }
      $sth_gene->finish();

      #Reset variables
      @sub_gene_list = ();
    }
  }
  return(\%maskedGene);
}


=head2 getGeneTerms()

=over 3

=item Function:

Get the associated external IDs associated with a particular gene ID (all GO terms for example)

=item Return:

Hash of ALEXA gene IDs, each associated with a hash of external IDs of the type specified

=item Note:

This will collapse terms to the gene level and remove duplicates (transcript specific info is lost)

=item Args:

'-dbh' => database handle

'-gene_ids' => \@gene_ids

'-id_type' => $type

Example 'id_types' include: CCDS, EMBL, EntrezGene, GO, HUGO, IPI, MIM, PDB, protein_id, Refseq_dna, Refseq_dna_predicted, Refseq_peptide, Refseq_peptide_predicted, UniGene, Uniprot/SPTREMBL, Uniprot/SWISSPROT

=item Example(s):

my $gene_ids_ref = &getGeneTerms ('-dbh'=>$dbh, '-gene_ids'=>\@gene_ids, '-id_type'=>$type);

=back

=cut


#######################
#getGeneTerms         #
#######################
sub getGeneTerms{

  my %args = @_;
  my $dbh = $args{'-dbh'};
  my @gene_ids = @{$args{'-gene_ids'}};
  my $type = $args{'-id_type'};


  print "\nProcessing multiple Gene IDs supplied as a list\n";

  my $gene_count = @gene_ids;

  #Make sure the user provided a $dbh
  unless ($dbh){
    print "\nCannot get genes without a valid -dbh\n\n";
    exit();
  }
  unless ($gene_count > 0){
    print "\n\nMust provide a valid ALEXA gene id or array of gene IDs for getGeneTerms()\n";
    exit();
  }
  unless ($type){
    print "\n\nMust provide a valid External ID type for getGeneTerms()\n";
    exit();
  }

  my $genes_processed = 0;
  my @sub_gene_list;

  my %gene_external_ids;


  foreach my $gene (@gene_ids){
    $genes_processed++;
    push (@sub_gene_list, $gene);

    #If the sub gene list contains 100 genes or the last gene has been reached submit the query and reset variables
    my $current_gene_count = @sub_gene_list;

    if ($current_gene_count == 100 || $genes_processed == $gene_count){

      #Disable print buffer to get these print statements to display immediately instead of on return()
      $| = 1;
      print ".";
      $| = 0;

      my $id_string = join(",", @sub_gene_list);
      $id_string =~ s/^\,//;

      #Submit query.  Since this is being done at the GENE level, this will give you duplicates for some cases where a gene has multiple transcripts
      #Each transcript shares many IDs with other transcripts, but they may also have unique ones.
      my $test = $dbh->{Name};
      my $db;
      my $sql_gene;
      if ($test =~ /ALEXA_hs_35_35h/){
	$sql_gene = "SELECT Gene.id,External_id.external_id FROM Gene,External_id WHERE External_id.fk_Gene__id=Gene.id AND Gene.id IN ($id_string) AND type = '$type';";
      }else{
	$sql_gene = "SELECT Gene.id,External_id.external_id FROM Gene,Transcript,External_id WHERE External_id.fk_Transcript__id=Transcript.id AND Transcript.fk_Gene__id=Gene.id AND Gene.id IN ($id_string) AND type = '$type'";
      }

      my $sth_gene = $dbh->prepare("$sql_gene");
      $sth_gene->execute();

      while (my ($gene_id,$external_id) = $sth_gene->fetchrow_array()){

	if ($gene_external_ids{$gene_id}){
	  my $ids_ref = $gene_external_ids{$gene_id}{external_ids};
	  $ids_ref->{$external_id}->{type} = $type;
	}else{
	  my %ids;
	  $ids{$external_id}{type} = $type;
	  $gene_external_ids{$gene_id}{external_ids} = \%ids;
	}
      }
      $sth_gene->finish();

      #Reset variables
      @sub_gene_list = ();
    }
  }

  return(\%gene_external_ids);
}


=head2 getProteinFeatures

=over 3

=item Function:

Get a list of protein features (domanin, TMDs, SPs, etc.) and their coordinates for a gene

=item Return:

Gene object as a hash

=item Args:

'-dbh' => database handle
'-gene_id' => '$gene_id'

=item Example(s):

my %gene_protein_features = %{&getProteinFeatures ('-dbh'=>$dbh, '-gene_ids'=>\@gene_ids)};

=back

=cut


#######################
#getProteinFeatures   #
#######################
sub getProteinFeatures{

  my %args = @_;
  my $dbh = $args{'-dbh'};
  my @gene_ids = @{$args{'-gene_ids'}};

  #Make sure the user provided a $dbh
  unless ($dbh){
    print "\nCannot get genes without a valid -dbh\n\n";
    exit();
  }
  my $gene_count = @gene_ids;
  unless ($gene_count > 0){
    print "\n\nMust provide a valid ALEXA gene id or array of gene IDs for getGeneInfo()\n";
    exit();
  }

  print "\nProcessing multiple Gene IDs supplied as a list\n";

  my $genes_processed = 0;
  my @sub_gene_list;

  my %gene_pfs;

  #Process 100 genes at a time
  foreach my $gene_id (@gene_ids){
    $genes_processed++;
    push (@sub_gene_list, $gene_id);

    #If the sub gene list contains 100 genes or the last gene has been reached submit the query and reset variables
    my $current_gene_count = @sub_gene_list;

    if ($current_gene_count == 100 || $genes_processed == $gene_count){
      #Disable print buffer to get these print statements to display immediately instead of on return()
      $| = 1;
      print ".";
      $| = 0;

      my $id_string = join(",", @sub_gene_list);
      $id_string =~ s/^\,//;

      #Submit query.
      my $sql_gene = "SELECT Gene.id,Protein_feature.id,Protein_feature.fk_Transcript__id,Protein_feature.type,Protein_feature.interpro_ac,Protein_feature.name,Protein_feature.program,Protein_feature.start,Protein_feature.end,Protein_feature.start_coords,Protein_feature.end_coords FROM Gene,Transcript,Protein_feature WHERE Protein_feature.fk_Transcript__id=Transcript.id AND Transcript.fk_Gene__id=Gene.id AND Gene.id IN ($id_string)";

      my $sth_gene = $dbh->prepare("$sql_gene");
      $sth_gene->execute();

      while (my ($gid,$id,$trans_id,$type,$interpro_ac,$name,$program,$start,$end,$start_coords,$end_coords) = $sth_gene->fetchrow_array()){

	my @start_coords = split(" ", $start_coords);
	my @end_coords = split(" ", $end_coords);

	if ($gene_pfs{$gid}){
	  my $trans_ref = $gene_pfs{$gid}{transcripts};

	  if ($trans_ref->{$trans_id}){
	    my $pfs_ref = $trans_ref->{$trans_id}->{protein_features};
	    $pfs_ref->{$id}->{trans_id} = $trans_id;
	    $pfs_ref->{$id}->{type} = $type;
	    $pfs_ref->{$id}->{interpro_ac} = $interpro_ac;
	    $pfs_ref->{$id}->{name} = $name;
	    $pfs_ref->{$id}->{program} = $program;
	    $pfs_ref->{$id}->{start} = $start;
	    $pfs_ref->{$id}->{end} = $end;
	    $pfs_ref->{$id}->{start_coords} = \@start_coords;
	    $pfs_ref->{$id}->{end_coords} = \@end_coords;

	  }else{
	    my %pfs;
	    $pfs{$id}{trans_id} = $trans_id;
	    $pfs{$id}{type} = $type;
	    $pfs{$id}{interpro_ac} = $interpro_ac;
	    $pfs{$id}{name} = $name;
	    $pfs{$id}{program} = $program;
	    $pfs{$id}{start} = $start;
	    $pfs{$id}{end} = $end;
	    $pfs{$id}{start_coords} = \@start_coords;
	    $pfs{$id}{end_coords} = \@end_coords;

	    $trans_ref->{$trans_id}->{protein_features} = \%pfs;
	  }
	}else{
	  my %transcripts;
	  my %pfs;
	  $pfs{$id}{trans_id} = $trans_id;
	  $pfs{$id}{type} = $type;
	  $pfs{$id}{interpro_ac} = $interpro_ac;
	  $pfs{$id}{name} = $name;
	  $pfs{$id}{program} = $program;
	  $pfs{$id}{start} = $start;
	  $pfs{$id}{end} = $end;
	  $pfs{$id}{start_coords} = \@start_coords;
	  $pfs{$id}{end_coords} = \@end_coords;

	  $transcripts{$trans_id}{protein_features} = \%pfs;
	  $gene_pfs{$gid}{transcripts} = \%transcripts;

	  #Create a list of non-redundant protein features for this gene
	  my @nr_pfs;
	  $gene_pfs{$gid}{nr_pfs} = \@nr_pfs;
	}

      }
      $sth_gene->finish();

      #Reset variables
      @sub_gene_list = '';
    }
  }

  #Before returning the object, identify a list of non-redundant PFs
  #Multiple overlaping transcripts from a single gene will have redundant motifs if their reading frame is the same
  my %pf_list;
  foreach my $gene_id (sort {$a <=> $b} keys %gene_pfs){
    my $trans_ref = $gene_pfs{$gene_id}{transcripts};
    foreach my $trans_id (sort {$a <=> $b} keys %{$trans_ref}){
      my $pfs_ref = $trans_ref->{$trans_id}->{protein_features};

      foreach my $pf_id (sort {$a <=> $b} keys %{$pfs_ref}){
	my $type = $pfs_ref->{$pf_id}->{type};
	my @starts = @{$pfs_ref->{$pf_id}->{start_coords}};
	my @ends = @{$pfs_ref->{$pf_id}->{end_coords}};

	my $starts_string = join('',@starts);
	my $ends_string = join('',@ends);

	unless ($pf_list{$type}{$starts_string}{$ends_string}){
	  $pf_list{$type}{$starts_string}{$ends_string}{tmp} = '';
	  push (@{$gene_pfs{$gene_id}{nr_pfs}}, $pf_id)
	}
      }
    }
  }

  return(\%gene_pfs);
}


=head2 getExonContent()

=over 3

=item Function:

Get the start/stop positions for the exon content of a gene (all expressed bases from any transcript of a gene)

Does not neccessarily represent any actual transcript

=item Return:

Exon Content object (keyed on expressed region count) as a hash

=item Args:

'-dbh' => database handle

'-gene_ids' => \@gene_ids

=item Example(s):

my $genes_ref = &getExonContent ('-dbh'=>$dbh, '-gene_ids'=>\@gene_ids);

my $exonContent_ref = $genes_ref->{$gene_id}->{exon_content}; #Get exon content for each gene

=back

=cut


#######################
#getExonContent       #
#######################
sub getExonContent{

  my %args = @_;
  my $dbh = $args{'-dbh'};
  my @gene_ids = @{$args{'-gene_ids'}};
  my $silent = $args{'-silent'};

  my %genes;

  my $gene_count = @gene_ids;

  #Make sure the user provided a $dbh
  unless ($dbh){
    print "\nMust provide a valid -dbh for getExonContent()\n\n";
    exit();
  }
  unless ($gene_count > 0){
    print "\n\nMust provide a valid ALEXA gene id or array of gene IDs for getExonContent()\n";
    exit();
  }
  unless($silent){
    $silent = "no";
  }

  #1.) Get the neccessary gene and exon data
  unless ($silent eq "yes"){
    print "\nProcessing multiple Gene IDs supplied as a list - Getting gene and exon coordinates\n";
  }

  my @sub_gene_list;
  my $genes_processed = 0;
  #Process 50 genes at a time
  foreach my $gene (@gene_ids){
    $genes_processed++;
    push (@sub_gene_list, $gene);

    #If the sub gene list contains 100 genes or the last gene has been reached submit the query and reset variables
    my $current_gene_count = @sub_gene_list;

    if ($current_gene_count == 100 || $genes_processed == $gene_count){

      #Disable print buffer to get these print statements to display immediately instead of on return()
      unless ($silent eq "yes"){
	$| = 1;
	print ".";
	$| = 0;
      }

      my $id_string = join(",", @sub_gene_list);
      $id_string =~ s/^\,//;

      #Submit query.

      #1-A.) First get the genomic boundaries for every gene
      my $sql_coords = "SELECT id,chromosome,chr_strand,chr_start,chr_end,gene_start,gene_end FROM Gene WHERE id IN ($id_string)";
      my $sth_coords = $dbh->prepare("$sql_coords");
      $sth_coords->execute();

      while (my ($id,$chromosome,$chr_strand,$chr_start,$chr_end,$gene_start,$gene_end) = $sth_coords->fetchrow_array()){

	$genes{$id}{chromosome} = $chromosome;
	$genes{$id}{chr_strand} = $chr_strand;
	$genes{$id}{chr_start} = $chr_start;
	$genes{$id}{chr_end} = $chr_end;
	$genes{$id}{gene_start} = $gene_start;
	$genes{$id}{gene_end} = $gene_end;

      }
      $sth_coords->finish();

      #1-B.) Next, get all exons for all transcripts for this particular gene id

      my %exons;
      my $sql_exon = "SELECT Transcript.fk_Gene__id,Exon.id,Exon.start,Exon.end FROM Exon,TranscriptExon,Transcript WHERE Transcript.id=TranscriptExon.fk_Transcript__id AND TranscriptExon.fk_Exon__id=Exon.id AND Transcript.fk_Gene__id IN ($id_string)";
      my $sth_exon = $dbh->prepare("$sql_exon");
      $sth_exon->execute();

      while (my ($id,$exon_id,$exon_start,$exon_end) = $sth_exon->fetchrow_array()){

	#Populate hash
	if ($genes{$id}{exons}){
	  my $exons_ref = $genes{$id}{exons};
	  $exons_ref->{$exon_id}->{exon_start} = $exon_start;
	  $exons_ref->{$exon_id}->{exon_end} = $exon_end;
	}else{
	  my %exons;
	  $exons{$exon_id}{exon_start} = $exon_start;
	  $exons{$exon_id}{exon_end} = $exon_end;
	  $genes{$id}{exons} = \%exons;
	}
      }
      $sth_exon->finish();

      #Reset variables
      @sub_gene_list = ();
    }
  }

  #Now go through each gene and use the data retrieved above to determine the exon content
  unless ($silent eq "yes"){
    print "\nDetermining exon content coordinates\n";
  }
  my $process_count = 0;
  foreach my $gene_id (keys %genes){
    $process_count++;

    if ($process_count == 100){
      unless ($silent eq "yes"){
	$| = 1;
	print ".";
	$| = 0;
      }
      $process_count = 0;
    }

    #Create a hash of the exon content positions for this gene
    my %exonContent;

    #2.) Build an array the same length as the number of bases in this gene (genomic sequence)
    #And initialize all positions in this array to '0'
    my $gene_start = $genes{$gene_id}{gene_start};
    my $gene_end = $genes{$gene_id}{gene_end};
    my $gene_length = ($gene_end - $gene_start)+1;
    my @positions;
    for (my $i = 1; $i <= $gene_length; $i++){
      push(@positions, 0);
    }

    #3.) Now go through each exon for this gene, and mark all the position in the gene that are expressed in any exon of any transcript
    my $exons_ref = $genes{$gene_id}{exons};
    foreach my $e_id (sort {$exons_ref->{$a}->{exon_start} <=> $exons_ref->{$b}{exon_start}} keys %{$exons_ref}){

      my $exon_length = ($exons_ref->{$e_id}->{exon_end} - $exons_ref->{$e_id}->{exon_start})+1;

      for (my $i=$exons_ref->{$e_id}->{exon_start}; $i <= $exons_ref->{$e_id}->{exon_end}; $i++){
	$positions[$i-1] = 1;
      }
    }

    #4.) Now go through the array of expressed positions and get the start/stop positions
    my $current_pos = 0;
    my $within_exon = 0;
    my $exon_end;

    my @starts;
    my @ends;

    #Get the start stop positions
    foreach my $pos (@positions){
      $current_pos++;

      if ($pos == 1 && $within_exon == 0){
	push(@starts, $current_pos);
	$within_exon = 1;

      }elsif ($pos == 1 && $within_exon == 1){
	$exon_end = $current_pos;

      }elsif ($pos == 0 && $within_exon == 1){
	$exon_end = $current_pos-1;
	push(@ends, $exon_end);
	$within_exon = 0;
      }
    }
    #If the final end position did not get added, add it now
    my $starts = @starts;
    my $ends = @ends;

    if ($starts == $ends){
      #do nothing
    }elsif($starts == ($ends + 1)){
      push (@ends, $exon_end);
    }else{
      print "\nIncorrect number of start/end coordinates!!\n\n";
      exit();
    }

    #5.) Now build an ordered hash of start/end positions
    my $exon_count = 0;
    foreach my $start (@starts){
      $exon_count++;
      my $end = shift @ends;
      $exonContent{$exon_count}{start} = $start;
      $exonContent{$exon_count}{end} = $end;
    }

    #Store the exonContent information hash for each gene
    $genes{$gene_id}{exon_content} = \%exonContent;
  }

  return(\%genes)
}


=head2 getIntronContent()

=over 3

=item Function:

Get the start/stop positions for the intron content of a gene (all NON-expressed bases from any transcript of a gene)

Each intron does not neccessarily represent one from a single actual transcript

=item Return:

Gene Intron Content object (keyed on expressed region count) as a hash

=item Args:

'-dbh' => database handle

'-gene_ids' => \@gene_ids

=item Example(s):

my $genes_ref = &getIntronContent ('-dbh'=>$dbh, '-gene_ids'=>\@gene_ids);

my $intronContent_ref = $genes_ref->{$gene_id}->{intron_content}; #Get exon content for each gene

=back

=cut


#######################
#getIntronContent       #
#######################
sub getIntronContent{

  my %args = @_;
  my $dbh = $args{'-dbh'};
  my @gene_ids = @{$args{'-gene_ids'}};

  my %genes;

  my $gene_count = @gene_ids;

  #Make sure the user provided a $dbh
  unless ($dbh){
    print "\nMust provide a valid -dbh for getIntronContent()\n\n";
    exit();
  }
  unless ($gene_count > 0){
    print "\n\nMust provide a valid ALEXA gene id or array of gene IDs for getIntronContent()\n";
    exit();
  }

  #1.) Get the neccessary gene and exon data
  print "\nProcessing multiple Gene IDs supplied as a list - Getting gene and exon coordinates\n";
  my @sub_gene_list;
  my $genes_processed = 0;
  #Process 50 genes at a time
  foreach my $gene (@gene_ids){
    $genes_processed++;
    push (@sub_gene_list, $gene);

    #If the sub gene list contains 100 genes or the last gene has been reached submit the query and reset variables
    my $current_gene_count = @sub_gene_list;

    if ($current_gene_count == 100 || $genes_processed == $gene_count){

      #Disable print buffer to get these print statements to display immediately instead of on return()
      $| = 1;
      print ".";
      $| = 0;

      my $id_string = join(",", @sub_gene_list);
      $id_string =~ s/^\,//;

      #Submit query.

      #1-A.) First get the genomic boundaries for every gene
      my $sql_coords = "SELECT id,chromosome,chr_strand,chr_start,chr_end,gene_start,gene_end FROM Gene WHERE id IN ($id_string)";
      my $sth_coords = $dbh->prepare("$sql_coords");
      $sth_coords->execute();

      while (my ($id,$chromosome,$chr_strand,$chr_start,$chr_end,$gene_start,$gene_end) = $sth_coords->fetchrow_array()){

	$genes{$id}{chromosome} = $chromosome;
	$genes{$id}{chr_strand} = $chr_strand;
	$genes{$id}{chr_start} = $chr_start;
	$genes{$id}{chr_end} = $chr_end;
	$genes{$id}{gene_start} = $gene_start;
	$genes{$id}{gene_end} = $gene_end;

      }
      $sth_coords->finish();

      #1-B.) Next, get all exons for all transcripts for this particular gene id

      my %exons;
      my $sql_exon = "SELECT Transcript.fk_Gene__id,Exon.id,Exon.start,Exon.end FROM Exon,TranscriptExon,Transcript WHERE Transcript.id=TranscriptExon.fk_Transcript__id AND TranscriptExon.fk_Exon__id=Exon.id AND Transcript.fk_Gene__id IN ($id_string)";
      my $sth_exon = $dbh->prepare("$sql_exon");
      $sth_exon->execute();

      while (my ($id,$exon_id,$exon_start,$exon_end) = $sth_exon->fetchrow_array()){

	#Populate hash
	if ($genes{$id}{exons}){
	  my $exons_ref = $genes{$id}{exons};
	  $exons_ref->{$exon_id}->{exon_start} = $exon_start;
	  $exons_ref->{$exon_id}->{exon_end} = $exon_end;
	}else{
	  my %exons;
	  $exons{$exon_id}{exon_start} = $exon_start;
	  $exons{$exon_id}{exon_end} = $exon_end;
	  $genes{$id}{exons} = \%exons;
	}
      }
      $sth_exon->finish();

      #Reset variables
      @sub_gene_list = ();
    }
  }

  #Now go through each gene and use the data retrieved above to determine the exon content
  print "\nDetermining exon content coordinates\n";
  my $process_count = 0;
  foreach my $gene_id (keys %genes){
    $process_count++;

    if ($process_count == 100){
      print ".";
      $process_count = 0;
    }

    #Create a hash of the exon content positions for this gene
    my %intronContent;

    #2.) Build an array the same length as the number of bases in this gene (genomic sequence)
    #And initialize all positions in this array to '1'
    my $gene_start = $genes{$gene_id}{gene_start};
    my $gene_end = $genes{$gene_id}{gene_end};
    my $gene_length = ($gene_end - $gene_start)+1;
    my @positions;
    for (my $i = 1; $i <= $gene_length; $i++){
      push(@positions, 1);
    }

    #3.) Now go through each exon for this gene, and mark all the position in the gene that are expressed in any exon of any transcript
    my $exons_ref = $genes{$gene_id}{exons};
    foreach my $e_id (sort {$exons_ref->{$a}->{exon_start} <=> $exons_ref->{$b}{exon_start}} keys %{$exons_ref}){

      my $exon_length = ($exons_ref->{$e_id}->{exon_end} - $exons_ref->{$e_id}->{exon_start})+1;

      for (my $i=$exons_ref->{$e_id}->{exon_start}; $i <= $exons_ref->{$e_id}->{exon_end}; $i++){
	$positions[$i-1] = 0;
      }
    }

    #4.) Now go through the array of expressed positions and get the start/stop positions
    my $current_pos = 0;
    my $within_intron = 0;
    my $intron_end;

    my @starts;
    my @ends;

    #Get the start stop positions
    foreach my $pos (@positions){
      $current_pos++;

      if ($pos == 1 && $within_intron == 0){
	push(@starts, $current_pos);
	$within_intron = 1;

      }elsif ($pos == 1 && $within_intron == 1){
	$intron_end = $current_pos;

      }elsif ($pos == 0 && $within_intron == 1){
	$intron_end = $current_pos-1;
	push(@ends, $intron_end);
	$within_intron = 0;
      }
    }
    #If the final end position did not get added, add it now
    my $starts = @starts;
    my $ends = @ends;

    if ($starts == $ends){
      #do nothing
    }elsif($starts == ($ends + 1)){
      push (@ends, $intron_end);
    }else{
      print "\nIncorrect number of start/end coordinates!!\n\n";
      exit();
    }

    #5.) Now build an ordered hash of start/end positions 
    my $intron_count = 0;
    foreach my $start (@starts){
      $intron_count++;
      my $end = shift @ends;
      $intronContent{$intron_count}{start} = $start;
      $intronContent{$intron_count}{end} = $end;
    }

    #Store the exonContent information hash for each gene
    $genes{$gene_id}{intron_content} = \%intronContent;
    $genes{$gene_id}{intron_count} = $intron_count;
  }

  return(\%genes)
}


=head2 junctionProbeCombinations()

=over 3

=item Function:

Get the number of theoretically possible valid exon-exon and intron-exon junction probes possible for a particular gene

Specify the max number of exons to consider in each exon skipping event (say 10 exons or less).

If '-exon_skip_limit' is not specified, no maximum will be considered and all valid combinations will be counted as theoretical probe combinations

=item Return:

Summary object (keyed on gene id) as a hash - Contains number of each type of probe

=item Args:

'-dbh' => database handle

'-gene_ids' => \@gene_ids

'-exon_skip_limit' => $exon_skip_threshold (say 10 exons)

=item Example(s):

my $probe_counts_ref = &junctionProbeCombinations ('-dbh'=>$dbh, '-gene_ids'=>\@gene_ids, '-exon_skip_limit'=>$exon_skip_threshold);

my $exon_exon = $probe_counts_ref->{$gene_id}->{exon_exon};

my $intron_exon = $probe_counts_ref->{$gene_id}->{intron_exon};

=back

=cut


############################
#junctionProbeCombinations #
############################
sub junctionProbeCombinations{

  my %args = @_;
  my $dbh = $args{'-dbh'};
  my @gene_ids = @{$args{'-gene_ids'}};
  my $silent = $args{'-silent'};

  unless($silent){
    $silent = "no";
  }
  unless ($silent eq "yes"){
    print "\nProcessing multiple Gene IDs supplied as a list\n";
  }

  my $gene_count = @gene_ids;

  #Make sure the user provided a $dbh
  unless ($dbh){
    print "\nCannot get genes without a valid -dbh\n\n";
    exit();
  }
  unless ($gene_count > 0){
    print "\n\nMust provide a valid ALEXA gene id or array of gene IDs for junctionProbeCombinations()\n";
    exit();
  }
  my $genes_processed = 0;
  my @sub_gene_list;

  my %gene_exons;

  #Process 100 genes at a time
  foreach my $gene (@gene_ids){
    $genes_processed++;
    push (@sub_gene_list, $gene);

    #If the sub gene list contains 100 genes or the last gene has been reached submit the query and reset variables
    my $current_gene_count = @sub_gene_list;

    if ($current_gene_count == 100 || $genes_processed == $gene_count){

      #Disable print buffer to get these print statements to display immediately instead of on return()
      unless ($silent eq "yes"){
	$| = 1;
	print ".";
	$| = 0;
      }
      my $id_string = join(",", @sub_gene_list);
      $id_string =~ s/^\,//;

      #Submit query.
      my $sql_gene = "SELECT Transcript.fk_Gene__id,Exon.id,ensembl_e_id,Exon.start,Exon.end FROM Exon,TranscriptExon,Transcript WHERE Transcript.id=TranscriptExon.fk_Transcript__id AND TranscriptExon.fk_Exon__id=Exon.id AND Transcript.fk_Gene__id IN ($id_string)";

      my $sth_gene = $dbh->prepare("$sql_gene");
      $sth_gene->execute();

      while (my ($gene_id,$exon_id,$ensembl_e_id,$exon_start,$exon_end) = $sth_gene->fetchrow_array()){

	#Populate hash
	if ($gene_exons{$gene_id}{exons}){
	  my $exons_ref = $gene_exons{$gene_id}{exons};
	  $exons_ref->{$exon_id}->{ensembl_e_id} = $ensembl_e_id;
	  $exons_ref->{$exon_id}->{exon_start} = $exon_start;
	  $exons_ref->{$exon_id}->{exon_end} = $exon_end;
	}else{
	  my %exons;
	  $exons{$exon_id}{ensembl_e_id} = $ensembl_e_id;
	  $exons{$exon_id}{exon_start} = $exon_start;
	  $exons{$exon_id}{exon_end} = $exon_end;
	  $gene_exons{$gene_id}{exons} = \%exons;
	}
      }
      $sth_gene->finish();

      #Reset variables
      @sub_gene_list = ();
    }
  }

  my $exon_skip_threshold = 0;
  if ($args{'-exon_skip_limit'}){
    $exon_skip_threshold = $args{'-exon_skip_limit'};
  }

  my %gene_probe_counts;

  #Now go through all of the exon for each gene and determine the number of possible junctions
  foreach my $gene_id (keys %gene_exons){

    #Get all exons for all transcripts for this particular gene id
    my $exons_ref = $gene_exons{$gene_id}{exons};

    #INTRON-EXON COMBINATIONS
    #Go through each exon and find the total number of unique start or end positions (junctions)
    my %intron_junctions;
    foreach my $exon_id (sort keys %{$exons_ref}){
    #Use a hash to compile a non-redundant list of start/end positions
    my $exon_start = $exons_ref->{$exon_id}->{exon_start};
    my $exon_end = $exons_ref->{$exon_id}->{exon_end};
    $intron_junctions{$exon_start}{tmp}='na';
    $intron_junctions{$exon_end}{tmp}='na';
  }

  #Total number of unique junctions (-2 for the first and last exon) will be the number of possible intron-exon junctions
  my $unique_intron_junctions = (keys %intron_junctions)-2;

  $gene_probe_counts{$gene_id}{intron_exon} = $unique_intron_junctions;

  #EXON-EXON COMBINATIONS
  #Go through each exon and find the total number of unique end positions
  #Also create a hash containing a non-redundant set of reference exons
  my %start_positions;
  my %end_positions;
  my %reference_exons;

  foreach my $exon_id (sort keys %{$exons_ref}){
    my $exon_start = $exons_ref->{$exon_id}->{exon_start};
    my $exon_end = $exons_ref->{$exon_id}->{exon_end};
    $start_positions{$exon_start}{tmp} = 'na';
    $end_positions{$exon_end}{tmp} = 'na';

    #Check each of the reference exons to see if one of these is the same, otherwise add it to the list
    my $redundant_exon = 0;
    foreach my $ref_exon_id (keys %reference_exons){
      if ($exon_start == $reference_exons{$ref_exon_id}{exon_start} && $exon_end == $reference_exons{$ref_exon_id}{exon_end}){
	$redundant_exon = 1;
      }
    }
    #Unless the current transcript exon was found to be redundant, add it to the list
    unless ($redundant_exon == 1){
      $reference_exons{$exon_id}{exon_start} = $exon_start;
      $reference_exons{$exon_id}{exon_end} = $exon_end;
    }
  }

  #Create arrays to represent a non-redundant set of reference exons
  my $ref_exon_count = keys %reference_exons;
  my @ref_exon_starts;
  my @ref_exon_ends;
  foreach my $ref_exon_id (sort {$reference_exons{$a}->{exon_start} <=> $reference_exons{$b}->{exon_start}} keys %reference_exons){
    push (@ref_exon_starts, $reference_exons{$ref_exon_id}{exon_start});
    push (@ref_exon_ends, $reference_exons{$ref_exon_id}{exon_end});
  }

  #Then find the number of valid connections between each of these end positions and all start positions
  my $valid_exon_connections = 0;
  foreach my $end_pos (sort {$a <=> $b} keys %end_positions){

    #For each end position, compare to all unique start positions of exons, any that are greater represent valid connections
    foreach my $start_pos (sort {$a <=> $b} keys %start_positions){
      if ($start_pos > $end_pos){

	#*********************ACCOUNT FOR THE EXON SKIPPING THRESHOLD SPECIFIED BY THE USER*********************
	#Determine how many exons are skipped by this exon1 ($end_pos) - exon2 ($start_pos) combination
	my $skipped_exons = 0;
	my $found_exon1 = 0;

	#First identify the last exon which corresponds to the exon1_end position
	for (my $i = 0; $i < $ref_exon_count-1; $i++){

	  #Look for the exon which corresponds to the exon2_start position
	  if ($start_pos >= $ref_exon_starts[$i] && $start_pos <= $ref_exon_ends[$i]){
	    last();
	  }

	  #Look for the exon which corresponds to the exon1_end position
	  if ($end_pos >= $ref_exon_starts[$i] && $end_pos <= $ref_exon_ends[$i]){
	    $found_exon1 = 1;
	    next();
	  }

	  #If exon1 has been found in the ordered reference set but exon2 has not been found yet, this is a skipped exon
	  #The boundaries of the current reference exon should also be completely between the two exons targetted
	  if ($found_exon1 == 1 && $end_pos < $ref_exon_starts[$i] && $start_pos > $ref_exon_starts[$i] && $end_pos < $ref_exon_ends[$i] && $start_pos > $ref_exon_ends[$i]){
	    $skipped_exons++;
	    next();
	  }
	}

	#This combination will be considered valid, but will not be counted if it involves skipping of more than the specified number of exons allowed
	if ($exon_skip_threshold == 0){
	  $valid_exon_connections++;
	}elsif($skipped_exons <= $exon_skip_threshold){
	  $valid_exon_connections++;
	}
      }
    }
  }
  $gene_probe_counts{$gene_id}{exon_exon} = $valid_exon_connections;
}
  return(\%gene_probe_counts);
}


=head2 displayGeneStats()

=over 3

=item Function:

Print out a variety of information about a gene and its transcripts, exons, etc.

=item Return:

Print results to standard output

=item Args:

'-dbh' => database handle

'-alexa_id' => $gene_id

'-ensembl_id' => $gene_id

=item Example(s):

&displayGeneStats ('-dbh'=>$dbh, '-alexa_id'=>$gene_id);

&displayGeneStats ('-dbh'=>$dbh, '-ensembl_id'=>$gene_id);

=back

=cut

############################
#displayGeneStats          #
############################
sub displayGeneStats{

  my %args = @_;
  my $dbh = $args{'-dbh'};
  my $alexa_id = $args{'-alexa_id'};
  my $ensembl_id = $args{'-ensembl_id'};

  #Make sure the user provided a $dbh
  unless ($dbh){
    print "\nCannot get genes without a valid -dbh\n\n";
    exit();
  }
  unless ($alexa_id || $ensembl_id){
    print "\n\nMust provide a valid ALEXA or Ensembl gene id for displayGeneStats()\n";
    exit();
  }

  my %gene;

  #1.) GET Gene Info
  if ($alexa_id){
    #Use the alexa gene id to get the related gene info
    my $sql_gene = "SELECT ensembl_g_id,species,gene_type,gene_name,evidence,chromosome,chr_start,chr_end FROM Gene WHERE id = '$alexa_id'";

    my $sth_gene = $dbh->prepare("$sql_gene");
    $sth_gene->execute();
    my ($ensembl_g_id,$species,$gene_type,$gene_name,$evidence,$chromosome,$chr_start,$chr_end) = $sth_gene->fetchrow_array();
    $sth_gene->finish();

    $gene{$alexa_id}{ensembl_g_id} = $ensembl_g_id;
    $gene{$alexa_id}{species} = $species;
    $gene{$alexa_id}{gene_type} = $gene_type;
    $gene{$alexa_id}{gene_name} = $gene_name;
    $gene{$alexa_id}{evidence} = $evidence;
    $gene{$alexa_id}{chromosome} = $chromosome;
    $gene{$alexa_id}{chr_start} = $chr_start;
    $gene{$alexa_id}{chr_end} = $chr_end;
  }else{
    #Use the ensembl gene id to get the related gene info
    my $sql_gene = "SELECT id,species,gene_type,gene_name,evidence,chromosome,chr_start,chr_end FROM Gene WHERE ensembl_g_id = '$ensembl_id'";
    my $sth_gene = $dbh->prepare("$sql_gene");
    $sth_gene->execute();
    my ($id,$species,$gene_type,$gene_name,$evidence,$chromosome,$chr_start,$chr_end) = $sth_gene->fetchrow_array();
    $sth_gene->finish();

    $alexa_id = $id;
    $gene{$alexa_id}{ensembl_g_id} = $ensembl_id;
    $gene{$alexa_id}{species} = $species;
    $gene{$alexa_id}{gene_type} = $gene_type;
    $gene{$alexa_id}{gene_name} = $gene_name;
    $gene{$alexa_id}{evidence} = $evidence;
    $gene{$alexa_id}{chromosome} = $chromosome;
    $gene{$alexa_id}{chr_start} = $chr_start;
    $gene{$alexa_id}{chr_end} = $chr_end;
  }

  #Print basic gene info
  print "\nALEXA: $alexa_id\tENSG: $gene{$alexa_id}{ensembl_g_id}\tType: $gene{$alexa_id}{gene_type}\tName: $gene{$alexa_id}{gene_name}\tEvidence: $gene{$alexa_id}{evidence}";
  #Print Gene size and position info
  my $gene_size = ($gene{$alexa_id}{chr_end} - $gene{$alexa_id}{chr_start})+1;
  print "\nChr: $gene{$alexa_id}{chromosome}\tPos: ($gene{$alexa_id}{chr_start} - $gene{$alexa_id}{chr_end})\tSize: $gene_size";

  #2.) GET External IDs
  #Get the external IDs of the type specified from ALEXA
  my $sql_entrez = "SELECT external_id FROM External_id WHERE fk_Gene__id = '$alexa_id' AND type = \"EntrezGene\"";
  my $sth_entrez = $dbh->prepare("$sql_entrez");
  $sth_entrez->execute();

  print "\nExternal IDs: ";

  my @entrez_ids;
  while (my ($external_id) = $sth_entrez->fetchrow_array()){
    push(@entrez_ids, $external_id);
  }
  $sth_entrez->finish();
  if (@entrez_ids){
    print "Entrez: @entrez_ids\t";
  }

  my $sql_refseq = "SELECT external_id FROM External_id WHERE fk_Gene__id = '$alexa_id' AND type = \"Refseq_dna\"";
  my $sth_refseq = $dbh->prepare("$sql_refseq");
  $sth_refseq->execute();

  my @refseq_ids;
  while (my ($external_id) = $sth_refseq->fetchrow_array()){
    push(@refseq_ids, $external_id);
  }
  $sth_refseq->finish();
  if (@refseq_ids){
    print "Refseq: @refseq_ids";
  }

  #3.) GET Transcript Info
  print "\nTranscripts: ";

  my $sql_transcript = "SELECT id,ensembl_t_id,start,end FROM Transcript WHERE fk_Gene__id = '$alexa_id'";
  my $sth_transcript = $dbh->prepare("$sql_transcript");
  $sth_transcript->execute();

  while (my ($t_id,$ensembl_t_id,$t_start,$t_end) = $sth_transcript->fetchrow_array()){
    #Print transcript info
    print "$ensembl_t_id ";

    #Get all the exons for this transcript and find the size of the transcript
    my $sql_transcript_exon = "SELECT Exon.id,Exon.start,Exon.end FROM Transcript,TranscriptExon,Exon WHERE Transcript.id=TranscriptExon.fk_Transcript__id AND TranscriptExon.fk_Exon__id=Exon.id AND Transcript.id = '$t_id'";

    my $sth_transcript_exon = $dbh->prepare("$sql_transcript_exon");
    $sth_transcript_exon->execute();

    my $transcript_size = 0;
    while (my ($e_id,$e_start,$e_end) = $sth_transcript_exon->fetchrow_array()){
      my $exon_size = ($e_end - $e_start)+1;
      $transcript_size += $exon_size;
    }
    $sth_transcript_exon->finish();
    print "($transcript_size bp)\t";
  }
  $sth_transcript->finish();



  #4.) GET Exon Info
  #Get all exons for all transcripts for this particular gene id
  my %exons;
  my $sql_exon = "SELECT Exon.id,ensembl_e_id,Exon.start,Exon.end FROM Exon,TranscriptExon,Transcript WHERE Transcript.fk_Gene__id = '$alexa_id' AND Transcript.id = TranscriptExon.fk_Transcript__id AND TranscriptExon.fk_Exon__id = Exon.id";
  my $sth_exon = $dbh->prepare("$sql_exon");
  $sth_exon->execute();

  while (my ($exon_id,$ensembl_e_id,$exon_start,$exon_end) = $sth_exon->fetchrow_array()){

    $exons{$exon_id}{ensembl_e_id} = $ensembl_e_id;
    $exons{$exon_id}{exon_start} = $exon_start;
    $exons{$exon_id}{exon_end} = $exon_end;

  }
  $sth_exon->finish();

  my $exon_count = keys %exons;

  print "\nTotal Exons: $exon_count\n";

  return(\%gene);
}


=head2 convertGeneCoordinates()

=over 3

=item Function:

Convert exon, transcript, or gene coordinates that are relative to the gene to chromosome coordinates

For a given gene ID, start position and end position, return the chromosome number, chromsome start and chromosome end

Note: Generally, all coordinates in ALEXA objects are relative to the gene (ie. 1 to length of genes genomic sequence)

=item Return:

Coordinate object (keyed on gene id) as a hash - Contains chromosome, chromosome start pos, and chromosome end pos

=item Args:

'-dbh' => database handle

'-gene_id' => $gene_id

'-start_pos' => $start

'-end_pos' => $end

=item Example(s):

my %coords = %{&convertGeneCoordinates ('-dbh'=>$dbh, '-gene_id'=>$gene_id, '-start_pos'=>$start, '-end_pos'=>$end)};

my $chr_start = $coords{$gene_id}{chr_start};

my $chr_end = $coords{$gene_id}{chr_end};

my $chromsome = $coords{$gene_id}{chr_name};

my $strand = $coords{$gene_id}{strand};

=back

=cut


############################
#convertGeneCoordinates    #
############################
sub convertGeneCoordinates{

  my %args = @_;
  my $dbh = $args{'-dbh'};
  my $gene_id = $args{'-gene_id'};
  my $start = $args{'-start_pos'};
  my $end = $args{'-end_pos'};

  my %coords;

  #Note: All gene coordinates stored in ALEXA are relative to the coding strand (i.e. start codon near beginning, end codon near the end)
  #This means that when coverting exon, or other coordinates back to the chromosome context, the original 'strand' of the gene on the chromosome
  #needs to be considered.
  my $sql_coords = "SELECT chromosome,chr_strand,chr_start,chr_end,gene_start,gene_end FROM Gene WHERE id = '$gene_id'";
  my $sth_coords = $dbh->prepare("$sql_coords");
  $sth_coords->execute();

  my ($chromosome,$chr_strand,$chr_start,$chr_end,$gene_start,$gene_end) = $sth_coords->fetchrow_array();

  $sth_coords->finish();

  #Make sure the supplied coordinates are actually within the specified gene
  unless ($start >= $gene_start-1 && $start <= $gene_end+1){
    print "\nStart coordinate ($start) provided to convertGeneCoordinates() does not appear valid for gene_id $gene_id\n\n";
    exit();
  }
  unless ($end >= $gene_start-1 && $end <= $gene_end+1){
    print "\nEnd coordinate ($end) provided to convertGeneCoordinates() does not appear valid for gene_id $gene_id\n\n";
    exit();
  }

  #Convert provided gene coordinates to coordinates relative to the chromosome
  if ($chr_strand == 1){
    my $query_chr_start = $chr_start + $start - 1;
    my $query_chr_end = $chr_start + $end - 1;

    $coords{$gene_id}{chr_start} = $query_chr_start;
    $coords{$gene_id}{chr_end} = $query_chr_end;
    $coords{$gene_id}{chr_name} = $chromosome;
    $coords{$gene_id}{strand} = "+";
  }elsif ($chr_strand == -1){

    my $query_chr_start = $chr_end - $end + 1;
    my $query_chr_end = $chr_end - $start + 1;

    $coords{$gene_id}{chr_start} = $query_chr_start;
    $coords{$gene_id}{chr_end} = $query_chr_end;
    $coords{$gene_id}{chr_name} = $chromosome;
    $coords{$gene_id}{strand} = "-";

  }else{
    print "\nStrand format: $chr_strand not understood by convertGeneCoordinates()!\n\n";
    exit();
  }

  return(\%coords);
}


=head2 convertGenomicCoordinates()

=over 3

=item Function:

Convert exon, transcript, or gene coordinates that are relative to the chromosome to gene coordinates

For a given gene ID, start position and end position, return the chromosome number, gene start and gene end

Note: Generally, all coordinates in ALEXA objects are already relative to the gene (ie. 1 to length of genes genomic sequence)

=item Return:

Coordinate object (keyed on gene id) as a hash - Contains chromosome, gene start pos, and gene end pos

=item Args:

'-dbh' => database handle

'-gene_id' => $gene_id

'-start_pos' => $start

'-end_pos' => $end

=item Example(s):

my %coords = %{&convertGenomicCoordinates ('-dbh'=>$dbh, '-gene_id'=>$gene_id, '-start_pos'=>$start, '-end_pos'=>$end)};

my $gene_start = $coords{$gene_id}{gene_start};

my $gene_end = $coords{$gene_id}{gene_end};

my $chromsome = $coords{$gene_id}{chr_name};

my $strand = $coords{$gene_id}{strand};

=back

=cut


############################
#convertGenomicCoordinates #
############################
sub convertGenomicCoordinates{

  my %args = @_;
  my $dbh = $args{'-dbh'};
  my $gene_id = $args{'-gene_id'};
  my $start = $args{'-start_pos'};
  my $end = $args{'-end_pos'};

  my %coords;

  #Note: All gene coordinates stored in ALEXA are relative to the coding strand (i.e. start codon near beginning, end codon near the end)
  #This means that when coverting exon, or other coordinates back to the chromosome context, the original 'strand' of the gene on the chromosome
  #needs to be considered.
  my $sql_coords = "SELECT id,chromosome,chr_strand,chr_start,chr_end,gene_start,gene_end FROM Gene WHERE id = '$gene_id'";

  #print "\nDEBUG SQL: $sql_coords";

  my $sth_coords = $dbh->prepare("$sql_coords");
  $sth_coords->execute();

  my ($id,$chromosome,$chr_strand,$chr_start,$chr_end,$gene_start,$gene_end) = $sth_coords->fetchrow_array();

  $sth_coords->finish();

  #Make sure the supplied coordinates are actually within the specified gene
  unless ($start >= $chr_start && $start <= $chr_end){
    print "\nStart coordinate ($start) provided to convertGenomicCoordinates() does not appear valid for gene_id $gene_id\n\n";
    exit();
  }
  unless ($end >= $chr_start && $end <= $chr_end){
    print "\nEnd coordinate ($end) provided to convertGenomicCoordinates() does not appear valid for gene_id $gene_id\n\n";
    exit();
  }

  #Chromosome coordinates should always be listed as start/end where start is smaller than end, regardless of strand
  unless ($start < $end){
    print "\nWhy is start_pos larger than end_pos in call to convertGenomicCoordinates()\n\n";
    exit();
  }

  #Convert provided genomic/chromosome coordinates to coordinates relative to the gene
  if ($chr_strand == 1){
    my $query_gene_start = $start - $chr_start + 1;
    my $query_gene_end = $end - $chr_start + 1;

    $coords{$gene_id}{gene_start} = $query_gene_start;
    $coords{$gene_id}{gene_end} = $query_gene_end;
    $coords{$gene_id}{chr_name} = $chromosome;
    $coords{$gene_id}{strand} = "+";
  }elsif ($chr_strand == -1){

    my $query_gene_start = $start - $chr_start + 1;
    my $query_gene_end = $end - $chr_start + 1;

    $coords{$gene_id}{gene_start} = $query_gene_start;
    $coords{$gene_id}{gene_end} = $query_gene_end;
    $coords{$gene_id}{chr_name} = $chromosome;
    $coords{$gene_id}{strand} = "-";

  }else{
    print "\nStrand format: $chr_strand not understood by convertGeneCoordinates()!\n\n";
    exit();
  }

  return(\%coords);
}


=head2 relativeGenePosition()

=over 3

=item Function:

Given a coordinate position for a particular ALEXA gene, determine the relative position of that coordinate within the gene

For example, for a coordinate representing the centre of a particular probe which represents a particular exon:

This function will determine where in the possible transcript of this gene this coordinate is found (distance from 5 prime end say)

Specifically it determines the amount of exon content (bases) upstream and downstream of the specified position given the known transcripts 

of the specified gene.

=item Return:

Coordinate object (keyed on gene id) as a hash - Contains chromosome, chromosome start pos, and chromosome end pos

=item Args:

'-dbh' => database handle

'-gene_id' => $gene_id

'-position' => $position

=item Example(s):

my %gene_pos = %{&relativeGenePosition ('-dbh'=>$dbh, '-gene_id'=>$gene_id, '-position'=>$position)};

my $distance_from_start = $gene_pos{$gene_id}{dist_5prime};

my $distance_from_end = $gene_pos{$gene_id}{dist_3prime};

my $total_expressed_bases = $gene_pos{$gene_id}{total_exon_bases};

=back

=cut


############################
#relativeGenePositions    #
############################
sub relativeGenePosition{

  my %args = @_;
  my $dbh = $args{'-dbh'};
  my $gene_id = $args{'-gene_id'};
  my $position = $args{'-position'};

  my %gene_pos;

  #First make sure the supplied coordinates are actually within the specified gene
  my $sql_coords = "SELECT id,chromosome,chr_strand,chr_start,chr_end,gene_start,gene_end FROM Gene WHERE id = '$gene_id'";
  print "\nDEBUG: $sql_coords";

  my $sth_coords = $dbh->prepare("$sql_coords");
  $sth_coords->execute();
  my ($id,$chromosome,$chr_strand,$chr_start,$chr_end,$gene_start,$gene_end) = $sth_coords->fetchrow_array();
  $sth_coords->finish();

  unless ($position >= $gene_start && $position <= $gene_end){
    print "\nCoordinate ($position) provided to relativeGenePosition() does not appear valid for gene_id $gene_id\n\n";
    exit();
  }

  #Next, get all exons for all transcripts for this particular gene id
  my %exons;
  my $sql_exon = "SELECT Exon.id,ensembl_e_id,Exon.start,Exon.end FROM Exon,TranscriptExon,Transcript WHERE Transcript.fk_Gene__id = '$gene_id' AND Transcript.id=TranscriptExon.fk_Transcript__id AND TranscriptExon.fk_Exon__id=Exon.id";
  my $sth_exon = $dbh->prepare("$sql_exon");
  $sth_exon->execute();

  while (my ($exon_id,$ensembl_e_id,$exon_start,$exon_end) = $sth_exon->fetchrow_array()){

    $exons{$exon_id}{ensembl_e_id} = $ensembl_e_id;
    $exons{$exon_id}{exon_start} = $exon_start;
    $exons{$exon_id}{exon_end} = $exon_end;

  }
  $sth_exon->finish();

  #Now build an array the same length as the number of bases in this gene (genomic sequence)
  #And initialize all positions in this array to '0'
  my $gene_length = ($gene_end - $gene_start)+1;
  my @positions;
  for (my $i = 1; $i <= $gene_length; $i++){
    push(@positions, 0);
  }

  #Now go through each exon, and mark all the positions in the gene that are expressed in any exon of any transcript
  foreach my $e_id (sort {$exons{$a}->{exon_start} <=> $exons{$b}{exon_start}} keys %exons){

    my $exon_length = ($exons{$e_id}{exon_end} - $exons{$e_id}{exon_start})+1;

    for (my $i=$exons{$e_id}{exon_start}; $i <= $exons{$e_id}{exon_end}; $i++){
      $positions[$i-1] = 1;
    }
  }

  #For the given input position, determine the number of exon bases upstream, downstream and the total number of exon bases
  my $total_exon_bases = 0;
  foreach my $base (@positions){
    $total_exon_bases += $base;
  }
  $gene_pos{$gene_id}{total_exon_bases} = $total_exon_bases;

  #Distance from 5' end of transcript
  my $distance_5prime = 0;
  for (my $i = 0; $i < $position; $i++){
    $distance_5prime += $positions[$i];
  }
  $gene_pos{$gene_id}{dist_5prime} = $distance_5prime;

  #Distance from 3' end of transcript
  my $distance_3prime = 0;
  for (my $i = $position-1; $i < $gene_length; $i++){
    $distance_3prime += $positions[$i];
  }
  $gene_pos{$gene_id}{dist_3prime} = $distance_3prime;

  return(\%gene_pos);
}


=head2 getMicroarrayProbes()

=over 3

=item Function:

Given a particular microarray design name (ex. 'Druggable_Genome'), return the probe IDs for all probes

=item Return:

Probe IDs or Probe Count IDs (which are not neccessarily the same - if in doubt use Probe Count IDs) and associated Gene IDs

Also a hash (keyed on probe_id) containing probe info such as (unit1 and unit2 start and end positions, probe type, etc.)

=item Args:

'-dbh' => database handle

'-name' => $microarray_name

'-id_type' => "probe_count_id" or "probe_id" - Determines what the hash will be keyed on

=item Example(s):

my %probe_ids = %{&getMicroarrayProbes ('-dbh'=>$dbh, '-name'=>"Druggable_Genome", '-id_type'=>"probe_count_id")};

=back

=cut


############################
#getMicroarrayProbes       #
############################
sub getMicroarrayProbes{

  my %args = @_;
  my $dbh = $args{'-dbh'};
  my $microarray_name = $args{'-name'};
  my $id_type = $args{'-id_type'};

  my %probe;

  #Make sure the user provided a $dbh
  unless ($dbh && $microarray_name){
    print "\nCannot getMicroarrayProbes() without a valid -dbh and microarray name\n\n";
    exit();
  }
  unless ($id_type eq "probe_count_id" || $id_type eq "probe_id"){
    print "\n\nMust choose to use 'probe_count_id' or 'probe_id' for getMicroarrayProbes()\n";
    exit();
  }

  #1.) Make sure a microarray with the specified name exists in the database
  my $sql_test = "SELECT id FROM Microarray WHERE name = '$microarray_name';";

  my $sth_test = $dbh->prepare("$sql_test");
  $sth_test->execute();
  my ($m_id) = $sth_test->fetchrow_array();
  $sth_test->finish();

  unless ($m_id){
    print "\n\ngetMicroarrayProbes() could not find a microarray with the name: $microarray_name!!\n\n";
    exit();
  }


  #2.) Get all the probes and neccessary info for the specified microarray
  my $sql = "SELECT Probe.id,Probe.probe_count,Probe.type,GeneProbe.fk_Gene__id,GeneProbe.unit1_start,GeneProbe.unit1_end,GeneProbe.unit2_start,GeneProbe.unit2_end FROM ArraySpot,Probe,Microarray,GeneProbe WHERE ArraySpot.fk_Probe__id=Probe.id AND ArraySpot.fk_Microarray__id=Microarray.id AND Microarray.name='Druggable_Genome' AND Probe.id=GeneProbe.fk_Probe__id;";



  my $sth = $dbh->prepare("$sql");
  $sth->execute();

  while (my ($id,$probe_count,$probe_type,$gene_id,$unit1_start,$unit1_end,$unit2_start,$unit2_end) = $sth->fetchrow_array()){

    if ($id_type eq "probe_count_id"){
      $probe{$probe_count}{id} = $id;
      $probe{$probe_count}{type} = $probe_type;
      $probe{$probe_count}{gene_id} = $gene_id;
      $probe{$probe_count}{unit1_start} = $unit1_start;
      $probe{$probe_count}{unit1_end} = $unit1_end;
      $probe{$probe_count}{unit2_start} = $unit2_start;
      $probe{$probe_count}{unit2_end} = $unit2_end;
    }else{
      $probe{$id}{probe_count_id} = $probe_count;
      $probe{$id}{type} = $probe_type;
      $probe{$id}{gene_id} = $gene_id;
      $probe{$id}{unit1_start} = $unit1_start;
      $probe{$id}{unit1_end} = $unit1_end;
      $probe{$id}{unit2_start} = $unit2_start;
      $probe{$id}{unit2_end} = $unit2_end;
    }
  }
  $sth->finish();

  return(\%probe);
}


=head2 selectFilteredSet()

=over 3

=item Function:

Determines what filtered sets of probes are available in an ALEXA instance and asks the user to select one

=item Return:

ALEXA filtered set ID (id from Probe_set table)

=item Args:

'-dbh' => database handle

=item Example(s):

my $filtered_set = &selectFilteredSet ('-dbh'=>$dbh);

=back

=cut


############################
#selectFilteredSet         #
############################
sub selectFilteredSet{
  my %args = @_;
  my $dbh = $args{'-dbh'};

  print BLUE, "\n\nSelecting a filtered probe set to use", RESET;

  my $filtered_set;
  my %sets;
  my $type = "Filtered_Probes";
  my $sql = "SELECT id,description FROM Probe_set WHERE type='$type';";

  my $sth = $dbh->prepare("$sql");
  $sth->execute();

  #Get all sets defined in the database of the type 'Filtered_Probe'
  while (my ($id,$description) = $sth->fetchrow_array()){
    $sets{$id}{description} = $description;
  }
  $sth->finish();

  my $sets_found = keys %sets;
  unless ($sets_found > 0){
    print RED, "\nCould not find any filtered sets of probes in the database!\n\n", RESET;
    exit;
  }

  #Display the possible choices to the user
  print BLUE, "\nFound the following filtered sets in the specified database:",RESET;
  my $first_set;

  my $last_set;
  foreach my $set (sort {$a <=> $b} keys %sets){
    print BLUE, "\n\t($set)\t$sets{$set}{description}", RESET;
    unless ($first_set){
      $first_set = $set;
    }
    $last_set = $set;
  }
  print "\n\n";

  if ($sets_found == 1){
    $filtered_set = $last_set;
    print YELLOW, "\tSelecting $filtered_set since it is the only option", RESET;
  }else{

    print YELLOW, "\tWhich set would you like to use ($first_set - $last_set): ", RESET;
    my $answer = <>;
    chomp($answer);

    unless ($sets{$answer}){
      print RED, "\nYour selection was not one of those listed!\n\n", RESET;
      exit();
    }
    $filtered_set = $answer;
  }

  print "\n\n";

  return($filtered_set);
}


1;


