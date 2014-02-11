#!/usr/bin/perl -w

#Written by Malachi Griffith

#NOTE: This entire analysis is based on build hg17

#1.) First get the list of ALEXA genes of interest from a design file containing info on all ALEXA probes
#    - At this time create a hash of ALEXA_PROBES(genes targeted -> probesets -> probes -> coordinates)
#    - Only grab Exon probes
#    - At this time also grab the strand and chromosome for the genes targeted by these probes from the appropriate ALEXA database
#    /home/malachig/AlternativeSplicing/Array_design/ALEXA_hs_35_35h/probe_dev/NimbleGen_Design_Submission/NimbleGenArray_mipVS5FUR_coords
#
#2.) The Affy exon array probesets which map to within an EnsEMBL gene are already known and defined in a mapfile
#    - For those genes identfied above, get all the relevant probesets
#    - Create a hash of AFFY_PROBES(genes -> probesets -> coordinates)
#    /home/malachig/AlternativeSplicing/perl_bin/Affy_analysis/mapfiles/hg17/COMPLETE_probesets_mapped_to_ensembl.txt
#
#3.) For all the Affy probesets identified in (2) get the individual probe IDs and coordinates from a liftOver file for this build
#    - Add these probes and coordinates to the AFFY_PROBES hash
#    /home/malachig/AlternativeSplicing/Array_analysis/Affymetrix_Exon_Arrays/Array_Design_and_Annotation_Files/HuEx-1_0-st_probe_tab-hg16/liftOverCoordinates/HuEx-1_0-st_probe_tab-hg17_liftOver.txt
#
#4.) Define Affy/ALEXA (AA) probesets
#    - For all genes identified in (1) get the coordinates for exons/exon regions (each of these will be a potential probeset)
#    - For each of these regions, create a probeset if at least one probe from either platform lies within it.
#    - Classify each of these regions as having probes from: AFFY, ALEXA, or BOTH
#    - Summarize the findings of these classes
#    - Count the number of AFFY and ALEXA probes for each exon region.
#    - Create a hash to store all this info.  AA_PROBES(genes -> exon_region/probeset -> coordinates, affy_probes, alexa_probes, class, probe_counts, etc.)
#
#5.) Print out a file summarizing each of the new AA probesets defined
#    - Print: AA_probeset_ID, chr, start, end, strand, Affy_probe_count, alexa_probe_count, class
#
#6.) Create datafiles
#    - For each platform, import a raw datafile containing triplicate MIP vs 5FUR data
#    - Print out each line to a new data file only if it belongs to an AA probeset with probes from BOTH platforms
#    - Keep the original probe ID but print out the new AA probeset ID for each probe data line
#    - The two resulting datafiles should contain the same number of probesets but will have different numbers of probe data entries

#Affy_data:
#/home/malachig/AlternativeSplicing/Array_analysis/Affymetrix_Exon_Arrays/MIP_vs_5FUR/MIP_vs_5FUR_8-26-2006/probelevel/HuEX-1_001-006_RawProbe.txt

#ALEXA_data:
#/home/malachig/AlternativeSplicing/Array_analysis/NimbleGen_Data/MIPvs5FUR_04Dec2006/AuxillaryData/PairData_formatted/Raw/All_hybes_withProbeInfo.txt

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);

use lib '/home/malachig/AlternativeSplicing/perl_bin/';
use utilities::utility qw(:all);
use utilities::ALEXA_DB qw(:all);

#Initialize command line options
my $database = '';
my $server = '';
my $user = '';
my $password = '';
my $design_file = '';
my $affy_ensembl_map_file = '';
my $affy_probe_tab_file = '';
my $affy_data_file = '';
my $alexa_data_file = '';
my $AA_summary_file = '';
my $AA_affy_data_file = '';
my $AA_alexa_data_file = '';

GetOptions ('database=s'=>\$database, 'server=s'=>\$server, 'user=s'=>\$user, 'password=s'=>\$password,
	    'design_file=s'=>\$design_file, 'affy_ensembl_map_file=s'=>\$affy_ensembl_map_file, 'affy_probe_tab_file=s'=>\$affy_probe_tab_file,
	    'affy_data_file=s'=>\$affy_data_file, 'alexa_data_file=s'=>\$alexa_data_file, 'AA_summary_file=s'=>\$AA_summary_file, 
	    'AA_affy_data_file=s'=>\$AA_affy_data_file, 'AA_alexa_data_file=s'=>\$AA_alexa_data_file);

#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\nExample: defineAffyNimbleGenProbesets.pl  --database=ALEXA_hs_35_35h  --server=server_name  --user=user_name  --password=pwd  --design_file=/home/malachig/AlternativeSplicing/Array_design/ALEXA_hs_35_35h/probe_dev/NimbleGen_Design_Submission/NimbleGenArray_mipVS5FUR_coords  --affy_ensembl_map_file=/home/malachig/AlternativeSplicing/perl_bin/Affy_analysis/mapfiles/hg17/COMPLETE_probesets_mapped_to_ensembl.txt  --affy_probe_tab_file=/home/malachig/AlternativeSplicing/Array_analysis/Affymetrix_Exon_Arrays/Array_Design_and_Annotation_Files/HuEx-1_0-st_probe_tab-hg16/liftOverCoordinates/HuEx-1_0-st_probe_tab-hg17_liftOver.txt  --affy_data_file=/home/malachig/AlternativeSplicing/Array_analysis/Affymetrix_Exon_Arrays/MIP_vs_5FUR/MIP_vs_5FUR_8-26-2006/probelevel/HuEX-1_001-006_RawProbe.txt  --alexa_data_file=/home/malachig/AlternativeSplicing/Array_analysis/NimbleGen_Data/MIPvs5FUR_04Dec2006/AuxillaryData/PairData_formatted/Raw/All_hybes_withProbeInfo.txt  --AA_summary_file=AffyAlexa_Probesets_Summary.txt  --AA_affy_data_file=AA_affy_data.txt  --AA_alexa_data_file=AA_alexa_data.txt\n\n", RESET;

#Make sure all options were specified
unless ($database && $server && $user && $password && $design_file && $affy_ensembl_map_file && $affy_probe_tab_file && $affy_data_file && $alexa_data_file && $AA_summary_file && $AA_affy_data_file && $AA_alexa_data_file){
  print RED, "\nOptions missing!\n\n", RESET;
  exit();
}

#1.) First get the list of ALEXA genes of interest from a design file containing info on all ALEXA probes
#    - At this time create a hash of ALEXA_PROBES(genes targeted -> probesets -> probes -> coordinates)
#    /home/malachig/AlternativeSplicing/Array_design/ALEXA_hs_35_35h/probe_dev/NimbleGen_Design_Submission/NimbleGenArray_mipVS5FUR_coords
my %alexa_probes;

#Establish connection with the Alternative Splicing Expression database
my $alexa_dbh = &connectDB('-database'=>$database, '-server'=>$server, '-user'=>$user, '-password'=>$password);

&getAlexaProbes('-design_file'=>$design_file, '-dbh'=>$alexa_dbh);

#print Dumper %alexa_probes;


#2.) The Affy exon array probesets which map to within an EnsEMBL gene are already known and defined in a mapfile
#    - For those genes identfied above, get all the relevant probesets
#    - Create a hash of AFFY_PROBES(genes -> probesets -> coordinates)
#    /home/malachig/AlternativeSplicing/perl_bin/Affy_analysis/mapfiles/hg17/COMPLETE_probesets_mapped_to_ensembl.txt
#
my %affy_probes;
my %affy_probeset_list;
&getAffyProbesets('-map_file'=>$affy_ensembl_map_file);

#3.) For all the Affy probesets identified in (2) get the individual probe IDs and coordinates from a liftOver file for this build
#    - Add these probes and coordinates to the AFFY_PROBES hash
#    /home/malachig/AlternativeSplicing/Array_analysis/Affymetrix_Exon_Arrays/Array_Design_and_Annotation_Files/HuEx-1_0-st_probe_tab-hg16/liftOverCoordinates/HuEx-1_0-st_probe_tab-hg17_liftOver.txt
&getAffyProbes('-probe_file'=>$affy_probe_tab_file);

#print Dumper %affy_probes;

#4.) Define Affy/ALEXA (AA) probesets
my %AA_probes;
my $AA_probeset_id = 0;
my $alexa_only_count = 0;
my $affy_only_count = 0;
my $both_count = 0;

&defineAffyAlexaProbesets;

#print Dumper %AA_probes;

#5.) Print out a file summarizing each of the new AA probesets defined
#    - Print: AA_probeset_ID, chr, start, end, strand, Affy_probe_count, alexa_probe_count, class
my %desired_alexa_probes;
my %desired_affy_probes;
&summarize_AA_probesets('-AA_summary_file'=>$AA_summary_file);

#6.) Create datafiles
#    - For each platform, import a raw datafile containing triplicate MIP vs 5FUR data
#    - Print out each line to a new data file only if it belongs to an AA probeset with probes from BOTH platforms
#    - Keep the original probe ID but print out the new AA probeset ID for each probe data line
#    - The two resulting datafiles should contain the same number of probesets but will have different numbers of probe data entries
&createAA_datafiles('-alexa_data_file'=>$alexa_data_file, '-affy_data_file'=>$affy_data_file,'-AA_alexa_data_file'=>$AA_alexa_data_file, '-AA_affy_data_file'=>$AA_affy_data_file);


#Close database connection
$alexa_dbh->disconnect();

exit();


###################################################################################################################
#1.) First get the list of ALEXA genes of interest from a design file containing info on all ALEXA probes         #
###################################################################################################################
sub getAlexaProbes{
  my %args = @_;
  my $design_file = $args{'-design_file'};
  my $dbh = $args{'-dbh'};

  my $exon_probe_count = 0;
  my $exon_probeset_count = 0;

  open (DESIGN, "$design_file") || die "\nCould not open design file: $design_file\n\n";

  my $first_line = 1;
  my %columns;

  while(<DESIGN>){
    chomp($_);
    my @line = split("\t", $_);

    #Watch for the header line
    if ($first_line == 1){
      $first_line = 0;
      my $column_count = 0;
      foreach my $column (@line){
	$columns{$column}{column_pos} = $column_count;
	$column_count++;
      }

      #Check for required columns
      unless($columns{'Probe_Count'} && $columns{'ProbeSet_ID'} && $columns{'Gene_ID'} && $columns{'Probe_Type'} && $columns{'Unit1_start_chr'} && $columns{'Unit1_end_chr'}){
	print RED, "\nExpected column missing\n\n", RESET;
	exit();
      }
      next();
    }

    my $probe_type = $line[$columns{'Probe_Type'}{column_pos}];
    unless ($probe_type eq "Exon"){
      next();
    }
    $exon_probe_count++;

    my $gene_id = $line[$columns{'Gene_ID'}{column_pos}];
    my $probe_id = $line[$columns{'Probe_Count'}{column_pos}];
    my $probeset_id = $line[$columns{'ProbeSet_ID'}{column_pos}];
    my $unit1_start = $line[$columns{'Unit1_start_chr'}{column_pos}];
    my $unit1_end = $line[$columns{'Unit1_end_chr'}{column_pos}];

    if ($alexa_probes{$gene_id}){
      #If this gene was previously observed 
      my $probesets_ref = $alexa_probes{$gene_id}{probesets};

      $alexa_probes{$gene_id}{gene_probe_count}++;  #Number of exon probes for this gene

      if ($probesets_ref->{$probeset_id}){
	#If this probeset was previously observed
	my $probes_ref = $probesets_ref->{$probeset_id}->{probes};
	$probes_ref->{$probe_id}->{start} = $unit1_start;
	$probes_ref->{$probe_id}->{end} = $unit1_end;

	$probesets_ref->{$probeset_id}->{probe_count}++;

      }else{
	#First time this probeset has been observed
	my %probes;
	$probes{$probe_id}{start} = $unit1_start;
	$probes{$probe_id}{end} = $unit1_end;

	$probesets_ref->{$probeset_id}->{probes} = \%probes;
	$probesets_ref->{$probeset_id}->{probe_type} = $probe_type;
	$probesets_ref->{$probeset_id}->{probe_count} = 1;
	$exon_probeset_count++;
      }
    }else{
      #First time gene was observed
      my %probes;
      my %probesets;

      $probes{$probe_id}{start} = $unit1_start;
      $probes{$probe_id}{end} = $unit1_end;

      $probesets{$probeset_id}{probes} = \%probes;
      $probesets{$probeset_id}{probe_type} = $probe_type;
      $probesets{$probeset_id}{probe_count} = 1;

      $exon_probeset_count++;

      $alexa_probes{$gene_id}{probesets} = \%probesets;
      $alexa_probes{$gene_id}{gene_probe_count} = 1;  #Number of probes that will be used to estimate gene expression
    }
  }

  close (DESIGN);

  #Now, for all genes gathered, get additional neccessary gene info from the alexa database
  my @gene_ids = keys %alexa_probes;

  my $gene_info_ref = &getGeneInfo ('-dbh'=>$dbh, '-gene_ids'=>\@gene_ids,'-sequence'=>"no");

  foreach my $gene_id (keys %alexa_probes){

    my $temp_chr = $gene_info_ref->{$gene_id}->{chromosome};
    $alexa_probes{$gene_id}{chromosome} = "chr"."$temp_chr";

    my $strand = $gene_info_ref->{$gene_id}->{chr_strand};
    if ($strand eq "-1"){
      $strand = "-";
    }else{
      $strand = "+";
    }

    $alexa_probes{$gene_id}{strand} = $strand;
  }

  my $gene_count = keys %alexa_probes;
  print BLUE "\nFound $gene_count genes, $exon_probeset_count exon-probesets and $exon_probe_count exon-probes in the ALEXA design\n\n", RESET;

  return();
}


###################################################################################################################
#2.) The Affy exon array probesets which map to within an EnsEMBL gene are already known and defined in a mapfile #
###################################################################################################################
sub getAffyProbesets{
  my %args = @_;
  my $map_file = $args{'-map_file'};

  my $potential_probesets = 0;

  open (MAP, "$map_file") || die "\nCould not open map file: $map_file\n\n";

  my $first_line = 1;
  my %columns;

  while(<MAP>){
    chomp($_);
    my @line = split("\t", $_);

    #Watch for the header line
    if ($first_line == 1){
      $first_line = 0;
      my $column_count = 0;
      foreach my $column (@line){
	$columns{$column}{column_pos} = $column_count;
	$column_count++;
      }

      #Check for required columns
      unless($columns{'probe_set_id'} && $columns{'alexa_gene_id'} && $columns{'chromosome'} && $columns{'strand'} && $columns{'start'} && $columns{'stop'}){
	print RED, "\nExpected column missing\n\n", RESET;
	exit();
      }
      next();
    }

    my $gene_id = $line[$columns{'alexa_gene_id'}{column_pos}];

    #Unless this probe belongs to one of the genes targeted on the ALEXA design, skip it.
    unless ($alexa_probes{$gene_id}){
      next();
    }
    $potential_probesets++;

    my $probeset_id = $line[$columns{'probe_set_id'}{column_pos}];
    my $chromosome = $line[$columns{'chromosome'}{column_pos}];
    my $strand = $line[$columns{'strand'}{column_pos}];
    my $start = $line[$columns{'start'}{column_pos}];
    my $end = $line[$columns{'stop'}{column_pos}];

    $affy_probeset_list{$probeset_id}{gene_id} = $gene_id;

    if ($affy_probes{$gene_id}){
      #If this gene was previously observed 
      my $probesets_ref = $affy_probes{$gene_id}{probesets};

      $affy_probes{$gene_id}{gene_probeset_count}++;  #Number of exon probes for this gene

      if ($probesets_ref->{$probeset_id}){
	#If this probeset was previously observed

      }else{
	#First time this probeset has been observed
	my %probes;

	$probesets_ref->{$probeset_id}->{probes} = \%probes;

      }
    }else{
      #First time gene was observed
      my %probes;
      my %probesets;

      $probesets{$probeset_id}{probes} = \%probes;


      $affy_probes{$gene_id}{probesets} = \%probesets;
      $affy_probes{$gene_id}{gene_probeset_count} = 1;  #Number of probes that will be used to estimate gene expression
    }
  }

  my $gene_count = keys %affy_probes;

  print BLUE "\nFound $gene_count correponding genes and $potential_probesets potential exon-probesets in the AFFY design\n\n", RESET;

  close (MAP);

  return();
}


####################################################################################################################################
#3.) For all the Affy probesets identified in (2) get the individual probe IDs and coordinates from a liftOver file for this build #
####################################################################################################################################
sub getAffyProbes{
  my %args = @_;
  my $probe_file = $args{'-probe_file'};

  my $potential_probes = 0;

  open (MAP, "$probe_file") || die "\nCould not open probe tab file: $probe_file\n\n";

  my $first_line = 1;
  my %columns;

  while(<MAP>){
    chomp($_);
    my @line = split("\t", $_);

    #Watch for the header line
    if ($first_line == 1){
      $first_line = 0;
      my $column_count = 0;
      foreach my $column (@line){
	$columns{$column}{column_pos} = $column_count;
	$column_count++;
      }

      #Check for required columns
      unless($columns{'seqname'} && $columns{'start'} && $columns{'stop'} && $columns{'Probe_ID'} && $columns{'Tm_Celsius'} && $columns{'strand'} && $columns{'Probe_Set_ID'}){
	print RED, "\nExpected column missing\n\n", RESET;
	exit();
      }
      next();
    }

    my $probeset_id = $line[$columns{'Probe_Set_ID'}{column_pos}];

    #Unless this probe belongs to one of the genes probesets identified in the previous step,  skip it.
    unless ($affy_probeset_list{$probeset_id}){
      next();
    }
    my $gene_id = $affy_probeset_list{$probeset_id}{gene_id};

    $potential_probes++;

    my $chromosome = $line[$columns{'seqname'}{column_pos}];
    my $start = $line[$columns{'start'}{column_pos}];
    my $end = $line[$columns{'stop'}{column_pos}];
    my $probe_id = $line[$columns{'Probe_ID'}{column_pos}];
    my $tm = $line[$columns{'Tm_Celsius'}{column_pos}];
    my $strand = $line[$columns{'strand'}{column_pos}];

    my $probesets_ref = $affy_probes{$gene_id}{probesets};
    my $probes_ref = $probesets_ref->{$probeset_id}->{probes};

    $probes_ref->{$probe_id}->{chromosome} = $chromosome;
    $probes_ref->{$probe_id}->{start} = $start;
    $probes_ref->{$probe_id}->{end} = $end;
    $probes_ref->{$probe_id}->{strand} = $strand;
    $probes_ref->{$probe_id}->{tm} = $tm;

  }

  print BLUE "\nFound $potential_probes potential exon-probes corresponding to the AFFY probesets imported above\n\n", RESET;

  close (MAP);

  return();
}


################################################################################################################################
#4.) Define Affy/ALEXA (AA) probesets                                                                                          #
################################################################################################################################
sub defineAffyAlexaProbesets{

  #Get a list of genes that have probes in both the ALEXA and AFFY designs
  my @gene_ids = sort {$a <=> $b} keys %affy_probes;

  print BLUE "\nCreating gene object to allow identification of exon regions\n\n", RESET;
  my $genes_ref = &createGeneObject ('-gene_ids'=> \@gene_ids);

  foreach my $gene_id (@gene_ids){

    print CYAN, "\n\tProcessing gene: $gene_id", RESET;

    my $exon_regions_found = 0;
    my $exon_regions_successfully_targeted = 0;

    #Go through each exon cluster - for clusters of a single exon get the probe sequences, for clusters of multiple exons, define regions
    my $exons_ref = $genes_ref->{$gene_id}->{exons};

    my $number_exons = keys %{$exons_ref};

    my $exon_clusters_ref = $genes_ref->{$gene_id}->{clusters};

    foreach my $cluster (sort {$a <=> $b} keys %{$exon_clusters_ref}){
      my $cluster_exons_aref = $exon_clusters_ref->{$cluster}->{exons};
      my $cluster_size = @{$cluster_exons_aref};

      #print BLUE, "\n\tProcessing cluster: $cluster consisting of $cluster_size exons", RESET;

      #1.) Single-exon clusters
      if ($cluster_size == 1){
	my $exon_id =@{$cluster_exons_aref}[0];
	my $region_start = $exons_ref->{$exon_id}->{exon_start};
	my $region_end = $exons_ref->{$exon_id}->{exon_end};
	#print BLUE, "\n\t\tSingle Exon Region: ($exons_ref->{$exon_id}->{exon_start} - $exons_ref->{$exon_id}->{exon_end}) covers exon ids: $exon_id", RESET;

	$exon_regions_found++;

	#Identify probes from each platform that fall within this region!!
	&identifyRegionProbes('-gene_id'=>$gene_id, '-region_start'=>$region_start, '-region_end'=>$region_end, '-exon_ids'=>$cluster_exons_aref);



      }else{
	#2.) Multi-exon clusters
	#Get the most informative regions from the cluster of overlapping exons (non-overlaping where possible)
	my $exon_regions_ref = &selectExonRegions('-exon_object'=>$exons_ref, '-exon_ids'=>$cluster_exons_aref);

	#For each region identified attempt to extract probe sequences as done for single-exon clusters
	foreach my $region (sort {$a <=> $b} keys %{$exon_regions_ref}){
	  my $region_start = $exon_regions_ref->{$region}->{region_start};
	  my $region_end = $exon_regions_ref->{$region}->{region_end};
	  my $region_exons_aref = $exon_regions_ref->{$region}->{exons};

	  #print BLUE, "\n\t\tExon Region: ($exon_regions_ref->{$region}->{region_start} - $exon_regions_ref->{$region}->{region_end}) covers exon ids: @{$region_exons_aref}", RESET;

	  $exon_regions_found++;

	  #Identify probes from each platform that fall within this region!!
	  &identifyRegionProbes('-gene_id'=>$gene_id, '-region_start'=>$region_start, '-region_end'=>$region_end, '-exon_ids'=>$region_exons_aref);

	}
      }
    }
  }

  my $AA_count = keys %AA_probes;
  print BLUE, "\n\n\nFound a total $AA_count exon regions with at least on Affy or ALEXA probe within it.  Breakdown as follows:", RESET;
  print BLUE, "\n\tALEXA only: $alexa_only_count\n\tAFFY only: $affy_only_count\n\tBOTH: $both_count\n\n", RESET;

  return();
}


################################################################################################
#Create gene object - return hash with gene sequence, masked sequence, exon info, etc.         #
#Identify overlapping exons and define exon clusters                                           #
################################################################################################
sub createGeneObject{
  my %args = @_;
  my @gene_ids = @{$args{'-gene_ids'}};

  my %gene_object;

  #Get the raw gene sequence
  print BLUE, "\nGetting gene info", RESET;

  my $gene_info_ref = &getGeneInfo ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids,'-sequence'=>"no");

  #Get all exons for this gene
  print BLUE, "\nGetting exon coordinate data", RESET;
  my $gene_exons_ref = &getExons ('-dbh'=>$alexa_dbh, '-gene_ids'=>\@gene_ids, '-sequence'=>"no");

  print BLUE, "\nConverting exon coordinates to chromosome system", RESET;

  #Go through all of these exons and convert their exon coordinates to chromosome coordinate system!!
  foreach my $gene_id (keys %{$gene_exons_ref}){

    my $exons_ref = $gene_exons_ref->{$gene_id}->{exons};

    foreach my $exon_id (keys %{$exons_ref}){
      my $gene_exon_start = $exons_ref->{$exon_id}->{exon_start};
      my $gene_exon_end = $exons_ref->{$exon_id}->{exon_end};

      my $coords_ref = &convertGeneCoordinates ('-dbh'=>$alexa_dbh, '-gene_id'=>$gene_id, '-start_pos'=>$gene_exon_start, '-end_pos'=>$gene_exon_end);

      $exons_ref->{$exon_id}->{exon_start} = $coords_ref->{$gene_id}->{chr_start};
      $exons_ref->{$exon_id}->{exon_end} = $coords_ref->{$gene_id}->{chr_end};

      #When chromosome coordinates are used, sometime the start position is a higher number than the end position (when original gene was on minus strand)
      #For simplicity make sure this is not the case
      #Not sure if this is neccessary ...
      my $temp_start = $exons_ref->{$exon_id}->{exon_start};
      my $temp_end = $exons_ref->{$exon_id}->{exon_end};

      if ($temp_start > $temp_end){
	$exons_ref->{$exon_id}->{exon_start} = $temp_end;
	$exons_ref->{$exon_id}->{exon_end} = $temp_start;
      }

    }
  }

  #Process each gene
  foreach my $gene_id (@gene_ids){

    my $exons_ref = $gene_exons_ref->{$gene_id}->{exons};

    my $exon_count = keys %{$exons_ref};

    #Go through each exon and create clusters of overlapping exons.  Any exon that shares at least one bp with another exon will be placed in the same cluster
    #An exon that is completely isolated from others will comprise a cluster of one
    my $cluster_count = 0;
    my %clusters;

  EXON: foreach my $exon_id (sort {$exons_ref->{$a}->{exon_start} <=> $exons_ref->{$b}->{exon_start}} keys %{$exons_ref}){

      #Get the start and end positions of the current exon
      my $exon_start = $exons_ref->{$exon_id}->{exon_start};
      my $exon_end = $exons_ref->{$exon_id}->{exon_end};

      #Go through each cluster and see if the current exon overlaps with one of the exons already present in one of these clusters
      foreach my $cluster_id (sort keys %clusters){
	my @tmp_exons = @{$clusters{$cluster_id}{exons}};
	foreach my $exon_id_test (@tmp_exons){

	  my $exon_start_test = $exons_ref->{$exon_id_test}->{exon_start};
	  my $exon_end_test = $exons_ref->{$exon_id_test}->{exon_end};

	  #See if the start of the current exon is within the range of the test exon
	  if ($exon_start >= $exon_start_test && $exon_start <= $exon_end_test){
	    push (@tmp_exons, $exon_id);
	    $clusters{$cluster_id}{exons} = \@tmp_exons;
	    next(EXON);
	  }
	  #See if the end of the current exon is within the range of the test exon
	  if ($exon_end >= $exon_start_test && $exon_end <= $exon_end_test){
	    push (@tmp_exons, $exon_id);
	    $clusters{$cluster_id}{exons} = \@tmp_exons;
	    next(EXON);
	  }
	  #See if the current exon completely flanks the test exon - if so it should be added to the cluster
	  #New condition to fix bug!!
	  if ($exon_start <= $exon_start_test && $exon_end >= $exon_end_test){
	    push (@tmp_exons, $exon_id);
	    $clusters{$cluster_id}{exons} = \@tmp_exons;
	    next(EXON);
	  }
	}
      }
      #If the current exon was not added to any of the current clusters - create a new cluster
      $cluster_count++;
      my @tmp_exons;
      push (@tmp_exons, $exon_id);
      $clusters{$cluster_count}{exons} = \@tmp_exons;
    }

    #Build the gene object from the info gathered above
    $gene_object{$gene_id}{ensembl_g_id} = $gene_info_ref->{$gene_id}->{ensembl_g_id};
    $gene_object{$gene_id}{chr_start} = $gene_info_ref->{$gene_id}->{chr_start};
    $gene_object{$gene_id}{chr_end} = $gene_info_ref->{$gene_id}->{chr_end};
    $gene_object{$gene_id}{number_exons} = $exon_count;
    $gene_object{$gene_id}{exons} = $exons_ref;
    $gene_object{$gene_id}{clusters} = \%clusters;
  }

  return(\%gene_object);
}


################################################################################################
#For exon clusters that have been derived from overlapping exons, identify informative regions #
#of these exons to use for probe design.                                                       #
#For each of these regions, make note of which of the exons from the cluster of overlapping    #
#exons are involved in each of the defined regions                                             #
################################################################################################
sub selectExonRegions{
  my %args = @_;

  my $exon_object_ref = $args{'-exon_object'};
  my $exon_ids_aref = $args{'-exon_ids'};

  #The exons in @exon_ids comprise a single cluster of overlapping exons
  #First try to identify regions of these exons that are unique to each exon

  #Compile a non-redundant list of start/end positions
  #Note that if two exons are in a cluster and have exactly the same start/end positions (which does happen in Ensembl!) then when
  #converted to a non-redundant list they will resolve to a single region which is good.
  my %junctions;
  my @junctions;
  foreach my $exon_id (@{$exon_ids_aref}){
    #Use a hash to compile a non-redundant list of start/end positions
    my $exon_start = $exon_object_ref->{$exon_id}->{exon_start};
    my $exon_end = $exon_object_ref->{$exon_id}->{exon_end};
    $junctions{$exon_start}{tmp}='na';
    $junctions{$exon_end}{tmp}='na';
    #print "\nEXON: $exon_id\tStart: $exon_start\tEnd: $exon_end";
 }
  #Create a sorted array of these junctions
  foreach my $junct (sort {$a <=> $b} keys %junctions){
    push (@junctions, $junct);
  }

  #Now consider the regions between the junctions and identify the exons associated with each
  my $number_junctions = @junctions;
  my %regions;
  my $region_count;

  for (my $i = 0; $i < $number_junctions-1; $i++){
    my $region_start = ($junctions[$i])+1;
    my $region_end = ($junctions[$i+1])-1;

    #print "\nCompare: $region_start - $region_end";

    #Confirm that the region selected is actually valid and larger than the required probe size
    if ($region_end <= $region_start){
      #print "\nExon junctions are too close to each other to allow probe design";
      next();
    }

    #Skip very small regions
    my $region_size = ($region_end - $region_start);
    if ($region_size <= 40){
      #print "\nExon junctions are too close to each other to allow probe design";
      next();
    }

    #Check each exon to see which overlap within this region at either end - or flank it completely
    my @exons_within_region;
    foreach my $exon_id (@{$exon_ids_aref}){
      my $exon_start = $exon_object_ref->{$exon_id}->{exon_start};
      my $exon_end = $exon_object_ref->{$exon_id}->{exon_end};

      #Is the start position of the selected region within this exon?
      if ($region_start >= $exon_start && $region_start <= $exon_end){
	push (@exons_within_region, $exon_id);
	next();
      }
      #Is the end position of the selected region within this exon?
      if ($region_end >= $exon_start && $region_end <= $exon_end){
	push (@exons_within_region, $exon_id);
	next();
      }
      #Does the selected region completely flank this exon?
      #New condition to fix bug!!
      if ($region_start <= $exon_start && $region_end >= $exon_end){
	push (@exons_within_region, $exon_id);
	next();
      }

    }
    $region_count++;
    $regions{$region_count}{region_start} = $region_start;
    $regions{$region_count}{region_end} = $region_end;
    $regions{$region_count}{exons} = \@exons_within_region;
    #print BLUE, "\n\t\tRegion: $region_count ($region_start - $region_end) covers exons: @exons_within_region", RESET;
  }

  #print Dumper %regions;

  #Keep track of the exons that are successfully covered by the regions defined
  #Check each region defined to see if it is completely within an exon.
  #For every exon, I want a region that is completely within the boundaries of the exon.  For those exons where this is not true,
  #define a region within the exon and note which other exons it covers - for this evaluation, if a probe covers any amount of sequence
  #of another exon it will be noted.  This is a conservative approach, if there is any ambiguity at all regarding which exons a probe covers,
  #it must be noted.  Nevertheless I want at least one probe that is completely within each exon.
  my %exons_covered;
  foreach my $region (sort {$a <=> $b} keys %regions){

    my $region_start = $regions{$region}{region_start};
    my $region_end = $regions{$region}{region_end};

    #Check this region against the initial list of exons for this exon cluster
    foreach my $exon_id (@{$exon_ids_aref}){
      my $exon_start = $exon_object_ref->{$exon_id}->{exon_start};
      my $exon_end = $exon_object_ref->{$exon_id}->{exon_end};

      #Which exons completely cover this region?
      if ($region_start >= $exon_start && $region_end <= $exon_end){
	$exons_covered{$exon_id}{tmp} = 'na';
      }
    }
  }

  return(\%regions);
}


######################################################################################################
#Given an exon-region start/end position, identify Affy and Alexa probes within this region          #
######################################################################################################
sub identifyRegionProbes{
  my %args = @_;
  my $gene_id = $args{'-gene_id'};
  my $region_start = $args{'-region_start'};
  my $region_end = $args{'-region_end'};
  my $exon_ids = $args{'-exon_ids'};

  my %temp_alexa_probes;
  my %temp_affy_probes;

  #1.) First go through the ALEXA probes
  my $alexa_probesets_ref = $alexa_probes{$gene_id}{probesets};

  foreach my $probeset_id (keys %{$alexa_probesets_ref}){

    my $probes_ref = $alexa_probesets_ref->{$probeset_id}->{probes};

    foreach my $probe_id (keys %{$probes_ref}){

      my $alexa_start = $probes_ref->{$probe_id}->{start};
      my $alexa_end = $probes_ref->{$probe_id}->{end};

      #test whether this probe is entirely within the current region
      if ($alexa_start >= $region_start && $alexa_start <= $region_end && $alexa_end >= $region_start && $alexa_end <= $region_end){
	$temp_alexa_probes{$probe_id}{start} = $probes_ref->{$probe_id}->{start};
	$temp_alexa_probes{$probe_id}{end} = $probes_ref->{$probe_id}->{end};
      }
    }
  }

  #2.) Now go through the Affy probes
  my $affy_probesets_ref = $affy_probes{$gene_id}{probesets};

  foreach my $probeset_id (keys %{$affy_probesets_ref}){

    my $probes_ref = $affy_probesets_ref->{$probeset_id}->{probes};

    foreach my $probe_id (keys %{$probes_ref}){

      my $affy_start = $probes_ref->{$probe_id}->{start};
      my $affy_end = $probes_ref->{$probe_id}->{end};

      #test whether this probe is entirely within the current region
      if ($affy_start >= $region_start && $affy_start <= $region_end && $affy_end >= $region_start && $affy_end <= $region_end){
	$temp_affy_probes{$probe_id}{start} = $probes_ref->{$probe_id}->{start};
	$temp_affy_probes{$probe_id}{end} = $probes_ref->{$probe_id}->{end};
      }
    }
  }

  #3.) If at least one probe was found on one of the platform, define an AA_probeset
  #4.) Classify the probeset as Affy-only, Alexa-only, or BOTH
  my $alexa_probes_found = keys %temp_alexa_probes;
  my $affy_probes_found = keys %temp_affy_probes;

  if ($alexa_probes_found > 0 || $affy_probes_found > 0){
    $AA_probeset_id++;

    if ($alexa_probes_found > 0 && $affy_probes_found > 0){
      $AA_probes{$AA_probeset_id}{gene_id} = $gene_id;
      $AA_probes{$AA_probeset_id}{region_start} = $region_start;
      $AA_probes{$AA_probeset_id}{region_end} = $region_end;
      $AA_probes{$AA_probeset_id}{alexa_probes} = \%temp_alexa_probes;
      $AA_probes{$AA_probeset_id}{affy_probes} = \%temp_affy_probes;
      $AA_probes{$AA_probeset_id}{alexa_probe_count} = $alexa_probes_found;
      $AA_probes{$AA_probeset_id}{affy_probe_count} = $affy_probes_found;
      $AA_probes{$AA_probeset_id}{exon_ids} = $exon_ids;
      $AA_probes{$AA_probeset_id}{class} = "BOTH";
      $both_count++;

    }elsif ($alexa_probes_found > 0){
      $AA_probes{$AA_probeset_id}{gene_id} = $gene_id;
      $AA_probes{$AA_probeset_id}{region_start} = $region_start;
      $AA_probes{$AA_probeset_id}{region_end} = $region_end;
      $AA_probes{$AA_probeset_id}{alexa_probes} = \%temp_alexa_probes;
      $AA_probes{$AA_probeset_id}{alexa_probe_count} = $alexa_probes_found;
      $AA_probes{$AA_probeset_id}{affy_probe_count} = 0;
      $AA_probes{$AA_probeset_id}{exon_ids} = $exon_ids;
      $AA_probes{$AA_probeset_id}{class} = "ALEXA-only";
      $alexa_only_count++;

    }elsif ($affy_probes_found > 0){
      $AA_probes{$AA_probeset_id}{gene_id} = $gene_id;
      $AA_probes{$AA_probeset_id}{region_start} = $region_start;
      $AA_probes{$AA_probeset_id}{region_end} = $region_end;
      $AA_probes{$AA_probeset_id}{affy_probes} = \%temp_affy_probes;
      $AA_probes{$AA_probeset_id}{affy_probe_count} = $affy_probes_found;
      $AA_probes{$AA_probeset_id}{alexa_probe_count} = 0;
      $AA_probes{$AA_probeset_id}{exon_ids} = $exon_ids;
      $AA_probes{$AA_probeset_id}{class} = "AFFY-only";
      $affy_only_count++;

    }
  }

  return();
}

#############################################################################################################
#5.) Print out a file summarizing each of the new AA probesets defined                                      #
#    - Print: AA_probeset_ID, chr, start, end, strand, Affy_probe_count, alexa_probe_count, class, exon_ids #
#############################################################################################################
sub summarize_AA_probesets{
  my %args = @_;
  my $AA_summary_file = $args{'-AA_summary_file'};

  open (SUMMARY, ">$AA_summary_file") || die "\nCould not open summary output file: $AA_summary_file\n\n";

  print SUMMARY "AA_probeset_ID\tgene_id\tchromosome\tregion_start\tregion_end\tstrand\tAffy_probe_count\tAlexa_probe_count\tclass\texon_ids\n";

  foreach my $aa_probeset_id (sort {$a <=> $b} keys %AA_probes){

    my $class = $AA_probes{$aa_probeset_id}{class};

    #Only summarize the common probesets which will be used for further analysis
    unless ($class eq "BOTH"){
      next();
    }

    my $gene_id = $AA_probes{$aa_probeset_id}{gene_id};
    my $chr = $alexa_probes{$gene_id}{chromosome};
    my $strand = $alexa_probes{$gene_id}{strand};

    print SUMMARY "$aa_probeset_id\t$gene_id\t$chr\t$AA_probes{$aa_probeset_id}{region_start}\t$AA_probes{$aa_probeset_id}{region_end}\t$strand\t$AA_probes{$aa_probeset_id}{affy_probe_count}\t$AA_probes{$aa_probeset_id}{alexa_probe_count}\t$AA_probes{$aa_probeset_id}{class}\t@{$AA_probes{$aa_probeset_id}{exon_ids}}\n";

    my $alexa_probes_ref = $AA_probes{$aa_probeset_id}{alexa_probes};
    my $affy_probes_ref = $AA_probes{$aa_probeset_id}{affy_probes};

    foreach my $alexa_probe_id (keys %{$alexa_probes_ref}){
      $desired_alexa_probes{$alexa_probe_id}{AA_probeset_id} = $aa_probeset_id;
    }
    foreach my $affy_probe_id (keys %{$affy_probes_ref}){
      $desired_affy_probes{$affy_probe_id}{AA_probeset_id} = $aa_probeset_id;
    }
  }

  close (SUMMARY);

  return();
}

##############################################################################################################################################
#6.) Create datafiles                                                                                                                        #
#    - For each platform, import a raw datafile containing triplicate MIP vs 5FUR data                                                       #
#    - Print out each line to a new data file only if it belongs to an AA probeset with probes from BOTH platforms                           #
#    - Keep the original probe ID but print out the new AA probeset ID for each probe data line                                              #
#    - The two resulting datafiles should contain the same number of probesets but will have different numbers of probe data entries         #
##############################################################################################################################################
sub createAA_datafiles{
  my %args = @_;
  my $alexa_data_file = $args{'-alexa_data_file'};
  my $affy_data_file = $args{'-affy_data_file'};
  my $AA_alexa_data_file = $args{'-AA_alexa_data_file'};
  my $AA_affy_data_file = $args{'-AA_affy_data_file'};

  #A.) Alexa data file
  print BLUE, "\nProcessing Alexa Data IN/OUT files\n\n", RESET;
  open (ALEXA_IN, "$alexa_data_file") || die "\nCould not open alexa_data_file: $alexa_data_file\n\n";
  open (ALEXA_OUT, ">$AA_alexa_data_file") || die "\nCould not open output alexa_data_file: $AA_alexa_data_file\n\n";

  my $alexa_header = 1;
  while (<ALEXA_IN>){
    chomp($_);

    if ($alexa_header == 1){
      print ALEXA_OUT "AA_probeset_id\t$_\n";
      $alexa_header = 0;
      next();
    }

    my @line = split ("\t", $_);
    my $probe_id = $line[1];

    if ($desired_alexa_probes{$probe_id}){
       print ALEXA_OUT "$desired_alexa_probes{$probe_id}{AA_probeset_id}\t$_\n";
    }
  }
  close (ALEXA_IN);
  close (ALEXA_OUT);

  #B.) Affy data file
  print BLUE, "\nProcessing Affy Data IN/OUT files\n\n", RESET;
  open (AFFY_IN, "$affy_data_file") || die "\nCould not open affy_data_file: $affy_data_file\n\n";
  open (AFFY_OUT, ">$AA_affy_data_file") || die "\nCould not open output affy_data_file: $AA_affy_data_file\n\n";

  my $affy_header = 1;
  while (<AFFY_IN>){
    chomp($_);

    #Skip comment lines
    if ($_ =~ /^\#/){
      next();
    }

    if ($affy_header == 1){
      print AFFY_OUT "AA_probeset_id\tAlexaGene_ID\t$_\n";
      $affy_header = 0;
      next();
    }

    my @line = split ("\t", $_);
    my $probe_id = $line[0];

    unless ($probe_id =~ /\d+/){
      next();
    }

    if ($desired_affy_probes{$probe_id}){
       my $aa_probeset_id = $desired_affy_probes{$probe_id}{AA_probeset_id};
       my $gene_id = $AA_probes{$aa_probeset_id}{gene_id};

       print AFFY_OUT "$aa_probeset_id\t$gene_id\t$_\n";

    }
  }

  close (AFFY_IN);
  close (AFFY_OUT);

  return();
}
