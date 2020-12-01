#!/usr/bin/perl -w

#Written by Malachi Griffith
#The purpose of this script is to select probe sequences within specified chromosomal regions
#The user will specify the desired target probe length (--target_length) and the degree of overlap/tiling (--overlap) (i.e. how closely to design probes to each other)
#If the value for --overlap is less than the --target_length, the designed probes will overlap each other
#At every design position, the length of the probe will be varied up to the --max_length_variance specified by the user to achieve the specified --target_tm
#If you want anisothermal probes of fixed length, simply set --max_variance_length=0

#In this way, many probes will be designed for each region.  From this pool, a selection of representative probes will be chosen

#In general the neccessary code can be summarized as follows:

# Import chromosomal regions of interest from EnsEMBL
# For each chromosomal region:
#     - Get the chromosomal sequence from EnsEMBL (both unmasked and masked)
#     - Extract probe sequences of length '--target_length' +/- '--max_length_variance', starting every '--overlap' bases as specified by the user
#     - Store the probe sequence, start position, and end position as a hash
#     - Eliminate probes that correspond to Genomic N's
#     - Eliminate probes that are comprised of too many repeat masked bases
#     - Eliminate probes that are outside the desired target Tm +/- the acceptable range
#     - Eliminate probes with simple repeats such as monopolymers and dipolymers

# Many probes will be designed across each region
#   - The end result of this script will be a file containing many probes for each region.  A seperate script will then evaluate these:
#   - Essentially we want the 'best' probe(s) for each region
#   - We want the probes with the highest specificity (ie. unique in the genome if possible)
#   - Also the 'best' probes are those that have a Tm that is as close as possible to the target Tm
#   - If there are multiple probes within x degrees of the Tm (where x is some reasonable range, say 5 degrees C) then chose the one closest to the centre
#   - Or alternatively come up with a scoring system to determine the probe that is close in Tm and close to the centre
#   - Store the Tm, coordinates, etc. and summarize for all chosen probe

use strict;
use Data::Dumper;
use Getopt::Long;
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
use utilities::Probes qw(:all);

#Initialize command line options
my $ensembl_api_version = ''; #Version of EnsEMBL to use
my $species = '';
my $ensembl_database = '';
my $ensembl_server = '';
my $ensembl_user = '';
my $ensembl_password = '';
my $region_file = '';
my $masked_bases_limit = '';
my $overlap = '';
my $target_tm = '';
my $tm_range = '';
my $target_length = '';
my $max_length_variance = '';
my $max_cycles = '';
my $verbose = '';
my $probe_dir = '';
my $outfile = '';
my $logfile = '';

GetOptions ('ensembl_api_version=s'=>\$ensembl_api_version, 'species=s'=>\$species,
	    'ensembl_database=s'=>\$ensembl_database, 'ensembl_server=s'=>\$ensembl_server,
	    'ensembl_user=s'=>\$ensembl_user, 'ensembl_password=s'=>\$ensembl_password,
	    'region_file=s'=>\$region_file, 'masked_bases_limit=s'=>\$masked_bases_limit,
	    'overlap=i'=>\$overlap, 'target_tm=f'=>\$target_tm, 'tm_range=f'=>\$tm_range, 'target_length=i'=>\$target_length,
	    'max_length_variance=i'=>\$max_length_variance, 'max_cycles=i'=>\$max_cycles,
	    'verbose=s'=>\$verbose, 'probe_dir=s'=>\$probe_dir, 'outfile=s'=>\$outfile,
	    'logfile=s'=>\$logfile);

#Provide instruction to the user
print GREEN, "\n\nNOTE: Before using this script, make sure the correct API version is hard coded!!\n\n", RESET;
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the correct EnsEMBL API version using: --ensembl_api_version (41, 42, etc.)", RESET;
print GREEN, "\n\tSpecify the SOURCE Ensembl Database, Server, User and Password using: --ensembl_database  --ensembl_server  --ensembl_user  and  --ensembl_password", RESET;
print GREEN, "\n\tThe species using: --species (e.g. --species=Human or --species='Homo sapiens')", RESET;
print GREEN, "\n\t\tMake sure the species you supply matches the EnsEMBL database you supply!!", RESET;
print GREEN, "\n\tSpecify a tab-delimited file containing chromosome coordinates of target regions using: --region_file", RESET;
print GREEN, "\n\t\tFormat of file. Should have a header line and each line should have the following columns first: 'unique_id, chromosome, start, end'", RESET;
print GREEN, "\n\tSpecify the full path to a directory containing fasta files for each chromosome in your input file", RESET;
print GREEN, "\n\t\tThis directory should contain one file for each masked chromosome and one for each unmasked chromosome", RESET;
print GREEN, "\n\tSpecify the maximum amount of repeat-masked bases allowed for a probe (as a percent of the length of the probe) using: --masked_bases_limit", RESET;
print GREEN, "\n\tSpecify the tiling distance or overlap between adjacent probe selection regions using: --overlap (eg. 5)", RESET;
print GREEN, "\n\tSpecify the target Tm for probe sequences using: --target_tm (eg. 67.0)", RESET;
print GREEN, "\n\tSpecify the Tm range allowed in case varying the length does not achieve the target Tm using: --tm_range (eg. 2.0)", RESET;
print GREEN, "\n\tSpecify the target probe length (must be a multiple of 2) using: --target_length (eg. 36)", RESET;
print GREEN, "\n\tSpecify the maximum this probe length will be allowed to vary using: --max_length_variance (eg. 5)", RESET;
print GREEN, "\n\tSpecify the maximum number of NimbleGen synthesis cycles allowed for each probe using: --max_cycles", RESET; 
print GREEN, "\n\t\tIf you want verbose output, use: --verbose=yes", RESET;
print GREEN, "\n\tSpecify the path to the directory of current probe files (used to determine starting probe ID value) using: --probe_dir", RESET;
print GREEN, "\n\tSpecify the name of the output probe file: --outfile", RESET;
print GREEN, "\n\tAlso specify the name of a log file to be used to store a record of the design process using: --logfile", RESET;
print GREEN, "\n\nExample: extractProbesFromEnsEMBL.pl  --ensembl_api_version=48  --species=Human  --ensembl_database=homo_sapiens_core_48_36j  --ensembl_server=ensembl01.bcgsc.ca  --ensembl_user=ensembl  --ensembl_password=pwd  --region_file=/projects/malachig/GTAC_Chip_Design/TargetRegions.txt  --masked_bases_limit=10.0  --overlap=1  --target_tm=67.0  --tm_range=10  --target_length=70  --max_length_variance=15  --max_cycles=180  --probe_dir=/projects/malachig/GTAC_Chip_Design/unfiltered_probes  --outfile=/projects/malachig/GTAC_Chip_Design/unfiltered_probes/regionProbes.txt  --logfile=/projects/malachig/GTAC_Chip_Design/logs/extractProbesFromEnsEMBL_LOG.txt\n\n", RESET;

#Make sure all options were specified
unless ($ensembl_api_version && $species && $ensembl_database && $ensembl_server && $ensembl_user && $ensembl_password && $region_file && $masked_bases_limit && $overlap && $target_tm && $tm_range && $target_length && ($max_length_variance || $max_length_variance eq "0") && ($max_cycles =~ /\d+/) && $probe_dir && $outfile && $logfile){
  print RED, "\nRequired input parameter(s) missing!\n\n", RESET;
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

##CONNECT TO ENSEMBL SERVER
my $registry = "Bio::EnsEMBL::Registry";
$registry->load_registry_from_db(-host =>$ensembl_server, -user =>$ensembl_user, -pass=>$ensembl_password, -db_version =>$ensembl_api_version);

#Get a connection to the local Ensembl CORE database
my $ensembl_core_api = Bio::EnsEMBL::Registry->get_DBAdaptor($species, "core");

my $slice_adaptor = $registry->get_adaptor($species, 'Core', 'Slice');

#Get the starting probe_id and probeset_id by examining the specified input probe file directory.  If it is empty start with 1
my @current_ids = &getCurrentProbeIds('-probe_dir'=>$probe_dir);
my $current_probe_id = $current_ids[0];          #Unique probe ID used for each successful probe
my $current_probeset_id = $current_ids[1];       #Unique count of exon-exon junctions with successful probes

open (LOG, ">$logfile") || die "\nCould not open logfile: $logfile\n\n";

#Print out the parameters supplied by the user to the logfile for future reference
print LOG "\nUser Specified the following options:\nensembl_api_version = $ensembl_api_version\nspecies = $species\nensembl_database = $ensembl_database\nensembl_server = $ensembl_server\nensembl_user = $ensembl_user\nensembl_password = pwd\nregion_file = $region_file\nmasked_bases_limit = $masked_bases_limit\noverlap = $overlap\ntarget_tm = $target_tm\ntm_range = $tm_range\ntarget_length = $target_length\nmax_length_variance = $max_length_variance\nverbose = $verbose\nprobe_dir = $probe_dir\noutfile = $outfile\nlogfile = $logfile\n\n";

#Import the regions of interest - perform basic sanity checks on the input file and organize regions accorinding their target sequence (chromosome)
my %regions;
my %sequence_names;
my $input_header;
&importTargetRegions('-region_file'=>$region_file);

#Test for region overlaps
&identifyRegionOverlaps('-region_object'=>\%regions);

#Keep track of the number of regions for which probes are designed.  This script will design many probes spanning each of these regions.
my $probe_number = 0 ;
my $total_regions = 0;
my $total_successful_regions = 0;

print YELLOW, "\n\nBegin extracting probe sequences", RESET;

print BLUE, "\nAll data will be written to $outfile\n\n", RESET;
print LOG "\nAll data will be written to $outfile\n\n";

#Open the probe output file
open (OUTFILE, ">$outfile") || die "\nCould not open $outfile";
print OUTFILE "Probe_Count\tProbeSet_ID\tRegion_ID\tSequence\tProbe_length\tProbe_Tm\tChromosome\tUnit1_start\tUnit1_end\tmasked_bases\tNimbleGenCycles\n";

#Counters for interests sake
my $total_extraction_attempts = 0;         #Total positions at which extraction of a probe was attempted (depends on number/size of regions, overlap, etc.)
my $ambigious_bases_failed_count = 0;      #Total positions where an ambiguous base was found (in at least one probe length attempted at that position)
my $repeat_bases_failed_count = 0;         #Total positions where too many repeat bases were found (in at least one probe length attempted at that position)
my $low_complexity_bases_failed_count = 0; #Total positions where a low complexity element was found (in at least one probe length attempted at that position)
my $cycles_failed_count = 0;               #Total positions where a probe that requires too many synthesis cycles was found (in at least one probe length ...)
my $tm_failed_count = 0;                   #Total positions which passed the other tests but failed to yield a probe with a Tm within the specified limit

my $one_percent_block_size; #Based on the size of the region and specified overlap - use to provide progress

#Go through each region, get corresponding source sequence (chromosome) and select probes for the desired region of that sequence
foreach my $chromosome (sort keys %sequence_names){

  print YELLOW, "\n\nExtracting probes for regions defined on chromosome: $chromosome\n", RESET;

  my $chr_regions_ref = $sequence_names{$chromosome}{regions};

  #Get a sequence object for this region (entire chromosome segment) from EnsEMBL
  my $sequences_ref = &getSeqsFromEnsembl('-chromosome'=>$chromosome);

  foreach my $region_id (@{$chr_regions_ref}){

    $total_regions++;

    my $region_size = $regions{$region_id}{region_size};

    my $temp = (($region_size/$overlap)/100);
    $one_percent_block_size = sprintf("%.0f", $temp);

    print BLUE, "\n\t$total_regions (REGION: $region_id)\tSize: $region_size\t... Extracting probe sequences\n", RESET;
    print LOG "\n\t$total_regions (REGION: $region_id)\tSize: $region_size\t... Extracting probe sequences\n";

    my $probes_ref;

    $probes_ref = &getRawProbes('-region_id'=>$region_id, '-chromosome'=>$chromosome, '-region_object'=>\%regions, '-sequence_object'=>$sequences_ref);
    my $probes_found = 0;

    if (%{$probes_ref}){
      $probes_found = keys %{$probes_ref};
    }
    print CYAN, "\n\t\tFound $probes_found probes for this region\n", RESET;
    print LOG "\n\t\tFound $probes_found probes for this region\n";

    if ($probes_found >= 1){
      $total_successful_regions++;
      &printProbeInfo('-probe_object'=>$probes_ref);
    }
  }
}

print YELLOW, "\n\nSUMMARY OF PROBE EXTRACTION ATTEMPTS", RESET;
print YELLOW, "\nTotal extraction attempts: $total_extraction_attempts", RESET;
print YELLOW, "\nTotal extraction attempts which encountered an ambiguous base: $ambigious_bases_failed_count", RESET;
print YELLOW, "\nTotal extraction attempts which encountered too many repeat masked bases: $repeat_bases_failed_count", RESET;
print YELLOW, "\nTotal extraction attempts which encountered a low complexity element: $low_complexity_bases_failed_count", RESET;
print YELLOW, "\nTotal extraction attempts which encountered probe with too many cycles: $cycles_failed_count", RESET;
print YELLOW, "\nTotal extraction attempts which passed these tests but failed to yield a probe with suitable Tm: $tm_failed_count", RESET;
print YELLOW, "\n\nSUMMARY OF PROBES PRINTED TO FILE",RESET;
print YELLOW, "\nTotal Regions Attempted: $total_regions\tTotal with at least one probe found: $total_successful_regions", RESET;
print YELLOW, "\nA total of $probe_number probes were printed to the output file\n\n", RESET;

print LOG "\n\nSUMMARY OF PROBE EXTRACTION ATTEMPTS";
print LOG "\nTotal extraction attempts: $total_extraction_attempts";
print LOG "\nTotal extraction attempts which encountered an ambiguous base: $ambigious_bases_failed_count";
print LOG "\nTotal extraction attempts which encountered too many repeat masked bases: $repeat_bases_failed_count";
print LOG "\nTotal extraction attempts which encountered a low complexity element: $low_complexity_bases_failed_count";
print LOG "\nTotal extraction attempts which encountered probe with too many cycles: $cycles_failed_count";
print LOG "\nTotal extraction attempts which passed these tests but failed to yield a probe with suitable Tm: $tm_failed_count";
print LOG "\n\nSUMMARY OF PROBES PRINTED TO FILE";
print LOG "\nTotal Regions Attempted: $total_regions\tTotal with at least one probe found: $total_successful_regions";
print LOG "\nA total of $probe_number probes were printed to the output file\n\n";


#Close the output file
close (OUTFILE);
close (LOG);
exit();


######################################################################################################################################
#Import the regions of interest                                                                                                      #
#- perform basic sanity checks on the input file and organize regions accorinding their target sequence (chromosome)                 #
######################################################################################################################################
sub importTargetRegions{
  my %args = @_;
  my $region_file = $args{'-region_file'};

  print YELLOW, "\n\nImporting target regions from: $region_file", RESET;
  print LOG "\n\nImporting target regions from: $region_file";

  open (REGIONS, "$region_file") || die "\nCould not open target region file: $region_file\n\n";

  my $first_line = 1;
  my $record_count = 0;

  while(<REGIONS>){
    chomp ($_);
    if ($first_line == 1){
      $first_line = 0;
      $input_header = $_;
      next();
    }
    $record_count++;
    my @line = split ("\t", $_);

    my $line = $_;
    my $region_id = $line[0];
    my $chromosome = $line[1];
    my $start = $line[2];
    my $end = $line[3];

    #Check data formats:
    unless ($region_id =~ /^\d+$/){print RED, "\nRegion ID format not correct: ($region_id) - check input file (line: $record_count + 1)!\n\n", RESET; exit();}
    unless ($chromosome =~ /^chr/){print RED, "\nChromosome name format not correct: ($chromosome) - check input file (line: $record_count + 1)!\n\n", RESET; exit();}
    unless ($start =~ /^\d+$/){print RED, "\nStart coordinate format not correct: ($start) - check input file (line: $record_count + 1)!\n\n", RESET; exit();}
    unless ($end =~ /^\d+$/){print RED, "\nEnd coordinate format not correct: ($end) - check input file (line: $record_count + 1)!\n\n", RESET; exit();}

    #Sanity checks
    if ($start >= $end){
      print RED, "\nStart coordinate ($start) must be smaller than end coordinate ($end) - check input file!\n\n", RESET;
      exit();
    }
    if ($regions{$region_id}){
      print RED, "\nRegion ID: $region_id appears to be a duplicate - check input file!\n\n", RESET;
      exit();
    }

    $regions{$region_id}{record_count} = $record_count;
    $regions{$region_id}{chromosome} = $chromosome;
    $regions{$region_id}{start} = $start;
    $regions{$region_id}{end} = $end;
    $regions{$region_id}{region_size} = ($end - $start)+1;
    $regions{$region_id}{line} = $line;

    if ($sequence_names{$chromosome}){
      $sequence_names{$chromosome}{count}++;
      push (@{$sequence_names{$chromosome}{regions}}, $region_id);
      my $region_list_ref = $sequence_names{$chromosome}{region_list};
      $region_list_ref->{$region_id}->{start} = $start;
    }else{
      my @regions;
      my %region_list;
      push (@regions, $region_id);
      $region_list{$region_id}{start} = $start;
      $sequence_names{$chromosome}{count} = 1;
      $sequence_names{$chromosome}{regions} = \@regions;
      $sequence_names{$chromosome}{region_list} = \%region_list;
    }

  }
  close (REGIONS);

  my $region_count = keys %regions;
  my $seq_count = keys %sequence_names;

  print BLUE, "\n\n\tFound $region_count regions corresponding to $seq_count sequences (chromosomes)\n\n", RESET;
  print LOG "\n\n\tFound $region_count regions corresponding to $seq_count sequences (chromosomes)\n\n";

  return();
}


######################################################################################################################
#Test for region overlaps                                                                                            #
######################################################################################################################
sub identifyRegionOverlaps{
  my %args = @_;
  my $regions_ref = $args{'-region_object'};
  my $regions_processed = 0;
  my $overlaping_regions = 0;

  print YELLOW, "\nChecking all regions against each other for overlaps\n", RESET;

  foreach my $chr (sort keys %sequence_names){

    print YELLOW, "\nProcessing regions on $chr", RESET;

    #Get the regions for this chromosome
    my $region_list_ref = $sequence_names{$chr}{region_list};

    #Define clusters of overlapping regions for this chromosome
    my %clusters;
    my $cluster_count = 0;

    #Sort the regions by start position and go through them one by one comparing coordinates to all previously considered regions
    foreach my $region_id (sort {$region_list_ref->{$a}->{start} <=> $region_list_ref->{$b}->{start}} keys %{$region_list_ref}){

      my $current_start = $regions_ref->{$region_id}->{start};
      my $current_end = $regions_ref->{$region_id}->{end};

      my $overlap_found = 0;
      my $overlaping_cluster = '';

      #Go through each cluster defined so far and check for overlaps to the current region
    CLUSTER:foreach my $cluster_id (sort {$a <=> $b} keys %clusters){

	my $ol_regions_ref = $clusters{$cluster_id}{ol_regions};

	foreach my $test_region_id (sort {$a <=> $b} keys %{$ol_regions_ref}){

	  my $test_start = $regions_ref->{$test_region_id}->{start};
	  my $test_end = $regions_ref->{$test_region_id}->{end};

	  if (($test_start > $current_start) && ($test_start < $current_end)){
	    $overlap_found = 1;
	  }
	  if (($test_end > $current_start) && ($test_end < $current_end)){
	    $overlap_found = 1;
	  }
	  if (($test_start <= $current_start) && ($test_end >= $current_end)){
	    $overlap_found = 1;
	  }
	  if ($overlap_found == 1){
	    $overlaping_cluster = $cluster_id;
	    last CLUSTER;
	  }
	}
      }

      if ($overlap_found == 1){
	#If this region does overlap one of the existing clusters, add it
	my $ol_regions_ref = $clusters{$overlaping_cluster}{ol_regions};
	$ol_regions_ref->{$region_id}->{start} = $regions_ref->{$region_id}->{start};
	$ol_regions_ref->{$region_id}->{end} = $regions_ref->{$region_id}->{end};

      }else{
	#If this region does not overlap any of the existing clusters, create a new cluster
	$cluster_count++;
	my %ol_regions;
	$ol_regions{$region_id}{start} = $regions_ref->{$region_id}->{start};
	$ol_regions{$region_id}{end} = $regions_ref->{$region_id}->{end};
	
	$clusters{$cluster_count}{ol_regions} = \%ol_regions;
	$clusters{$cluster_count}{chromosome} = $chr;
	
      }
    }

    #Go through each of the clusters defined for this chromosome and print a summary of it
    foreach my $cluster_id (sort {$a <=> $b} keys %clusters){

      #If a cluster consists of only a single region, simply keep the original record
      my $ol_regions_ref = $clusters{$cluster_id}{ol_regions};
      my $cluster_size = keys %{$ol_regions_ref};
      $regions_processed += $cluster_size;

      if ($cluster_size == 1){
	#Skip clusters consisting of only a single region
	next();
	
      }else{
	$overlaping_regions++;
	print YELLOW, "\n\n\tCluster: $overlaping_regions consists of the following overlaping regions", RESET;
	print LOG "\n\nCluster: $overlaping_regions consists of the following overlaping regions";

	foreach my $region_id (sort {$ol_regions_ref->{$a}->{start} <=> $ol_regions_ref->{$b}->{start}} keys %{$ol_regions_ref}){
	  my $size = ($ol_regions_ref->{$region_id}->{end} - $ol_regions_ref->{$region_id}->{start}) + 1;
	  print YELLOW, "\n\t\tRegion $region_id ($chr:$ol_regions_ref->{$region_id}->{start}-$ol_regions_ref->{$region_id}->{end} Size = $size)", RESET;
	  print LOG "\n\tRegion $region_id ($chr:$ol_regions_ref->{$region_id}->{start}-$ol_regions_ref->{$region_id}->{end} Size = $size)";

	}
      }
    }
  }

  print YELLOW, "\n\nProcessed $regions_processed regions and found a total of $overlaping_regions clusters of overlapping regions\n", RESET;
  print LOG "\n\nProcessed $regions_processed regions and found a total of $overlaping_regions clusters of overlapping regions\n";

  if ($overlaping_regions > 0){
    print YELLOW, "\nDo you wish to proceed (y/n)? ", RESET;
    my $answer = <>;
    chomp ($answer);
    unless ($answer =~ /^y$|^yes$/i){
      print RED, "\nAborting\n\n", RESET;
      exit();
    }
  }
  return();
}


##################################################################################################
#Get a sequence object for this region (entire chromosome segment) from EnsEMBL                  #
##################################################################################################
sub getSeqsFromEnsembl{
  my %args = @_;
  my $chromosome = $args{'-chromosome'};

  my %seqs;

  my $chr_format;
  if ($chromosome =~ /chr(.*)/){
    $chr_format = $1;
  }else{
    print RED, "\nChromosome format not understood: $chromosome\n\n", RESET;
    exit();
  }

  #Retrieve a slice representing the entire chromosome specified
  print BLUE, "\nFetching entire sequence for chromosome $chr_format from Ensembl: $ensembl_database\n\n", RESET;
  my $slice = $slice_adaptor->fetch_by_region('chromosome', $chr_format);

  #Get the unmasked sequence
  my $chr_seq = $slice->seq();
  my $upper_chr_seq = uc($chr_seq);

  #Get the entire sequence in hard-masked form (N's representing all bases identified by RepeatMasker, Tandem Repeat Finder, and Dust)
  my $hard_masked_chr_seq_object = $slice->get_repeatmasked_seq();
  my $hard_masked_chr_seq = $hard_masked_chr_seq_object->seq();
  my $upper_masked_chr_seq = uc($hard_masked_chr_seq);

  my $chr_length = length($upper_chr_seq);
  my $chr_masked_length = length($upper_masked_chr_seq);

  unless ($chr_length == $chr_masked_length){
    print RED, "\nChromosome seq and masked version do not have matching lengths!\n\n", RESET;
    exit();
  }
  my $masked_bases = ($upper_masked_chr_seq =~ tr/N/N/);


  my $masked_percent = ($masked_bases/$chr_length)*100;
  my $masked_percent_f = sprintf("%.2f", $masked_percent);
  print BLUE, "\nRetrieved $chromosome sequence.  Length: $chr_length\tPercent Masked: $masked_percent_f%\n", RESET;

  $seqs{$chromosome}{unmasked} = \$upper_chr_seq;
  $seqs{$chromosome}{masked} = \$upper_masked_chr_seq;

  return (\%seqs);
}


################################################################################################
#Given a sequence object and a sequence start and stop position, return valid probes.
#Disqualify probes that have genomic N's, too many repeat masked N's, and Tm outside
#the accepted range.  Return a hash of probes sequences with start/end coordinates, Tm, etc.
################################################################################################
sub getRawProbes{
  my %args = @_;
  my $region_id = $args{'-region_id'};
  my $chromosome = $args{'-chromosome'};
  my $region_object_ref = $args{'-region_object'};
  my $sequence_object_ref = $args{'-sequence_object'};

  my $seq_ref = $sequence_object_ref->{$chromosome}->{unmasked};
  my $masked_seq_ref = $sequence_object_ref->{$chromosome}->{masked};

  my $start_pos = $region_object_ref->{$region_id}->{start};
  my $end_pos = $region_object_ref->{$region_id}->{end};

  my $probe_count = 0;
  my %probe_object;

  #Probe characteristics defined globally ($target_length, $target_tm, $max_length_variance, $overlap)

  #If the desired probe length is longer than the sequence provided there is no point in continuing
  my $sequence_length = ($end_pos+1) - $start_pos;
  if (($target_length + $max_length_variance) > $sequence_length){
    return(\%probe_object);
  }


  #Define last position in the exon region that it is safe to select probes from without running past
  my $start_cutoff = ($end_pos - ($target_length + $max_length_variance + 1));

  my $counter = 0;

  #Parse probe sequences of $probe_length, starting every $overlap bases
 POS: for (my $i = $start_pos-1; $i < $start_cutoff; $i += $overlap){
    $total_extraction_attempts++;
    $counter++;

    if ($counter == $one_percent_block_size){
      $counter = 0;
      $| = 1;
      print BLUE, ".", RESET;
      $| = 0;
    }

    my %position_probes_object;
    my $position_probe_count = 0;  #Count of probes with varying length at a single position

    my $ambigious_bases_failed = 0;
    my $repeat_bases_failed = 0;
    my $low_complexity_bases_failed = 0;
    my $cycles_failed = 0;

  LENGTH: for (my $j = $target_length-$max_length_variance; $j <= $target_length+$max_length_variance; $j++){
      # Note: $i is the current probe start position and $j is the current probe length

      my $probe_seq = substr (${$seq_ref}, $i, $j);

      #Get the start and end position of the probe within the gene
      my $probe_seq_start = $i + 1;
      my $probe_seq_end = $i + $j;

      #Each probe must now pass through a number of filters to make it into the final list of candidates
      #If a probe fails any of the following basic checks - skip to the next probe

      #1.) Check for presence of genomic N's and other non-valid letters such as Ambiguiety codes
      #    - Watch for sequences with any single character other than A,T,C,G
      if($probe_seq =~ /[^ATCG]/){
	if ($verbose eq "yes"){
	  print YELLOW, "\n\t\tThis probe: $probe_seq contains an invalid bases such as an N or R", RESET;
	  print LOG "\n\t\tThis probe: $probe_seq contains an invalid bases such as an N or R";
	}
	$ambigious_bases_failed = 1;

	next(LENGTH); #Note if these are found, making the probe longer wont help, but we dont want to abandon any probes that were already successful
      }

      #2.) Check for presence of excessive repeatMasked bases
      #Specify the number of allowable RepeatMasked bases - Experiment with thresholds 
      my $allowed_masked_bases = ($j*($masked_bases_limit/100));  #If more than 1/4 of probe is masked, it will be rejected
      my $masked_seq = substr (${$masked_seq_ref}, $i, $j);
      my @n_count = $masked_seq =~ m/N/g;
      my $n_count = @n_count;
      if ($n_count > $allowed_masked_bases){
	if ($verbose eq "yes"){
	  print YELLOW, "\n\t\tToo many masked bases ($n_count)", RESET;
	  print LOG "\n\t\tToo many masked bases ($n_count)";
	}
	$repeat_bases_failed = 1;

	next(LENGTH); #Note making the probe longer might help in this case, but we dont want to abandon any probes that were already successful
      }

      #3.) Calculate the Tm of this probe
      my $temp_k = &tmCalc('-sequence'=>$probe_seq, '-silent'=>1);
      my $Tm_celsius = &tmConverter('-tm'=>$temp_k, '-scale'=>'Kelvin');

      #4.) Check for probes with simple repeats
      my $test_seq = $probe_seq;
      #Search for Monopolymeric repeats of 6 or longer
      if ($test_seq =~ /(A{6,}|T{6,}|C{6,}|G{6,})/){
	if ($verbose eq "yes"){
	  print YELLOW, "\n\t\tMonopolymer detected ($1): $probe_seq", RESET;
	  print LOG "\n\t\tMonopolymer detected ($1): $probe_seq";
	}
	$low_complexity_bases_failed = 1;

	next(LENGTH); #Note if these are found, making the probe longer wont help, but we dont want to abandon any probes that were already successful
      }
      #Search for dipolymeric repeats of 4 or longer
      my $match = "AT";
      my @matches = qw (AT AG AC TA TG TC CA CT CG GA GT GC);
      foreach my $match (@matches){
	if ($test_seq =~ /($match){4,}/g){
	  if ($verbose eq "yes"){
	    print YELLOW, "\n\t\tDipolymer detected ($1): $probe_seq", RESET;
	    print LOG "\n\t\tDipolymer detected ($1): $probe_seq";
	  }
	  $low_complexity_bases_failed = 1;

	  next(LENGTH); #Note if these are found, making the probe longer wont help, but we dont want to abandon any probes that were already successful
	}
      }

      #5.) Check for probes that require more synthesis cycles than the cutoff specified by the user
      my $cycles = &calculate_nimblegen_cycles('-sequence'=>$probe_seq);
      if ($cycles > $max_cycles){
	if ($verbose eq "yes"){
	  print YELLOW, "\n\t\tToo many cycles required to synthesize this probe ($cycles)", RESET;
	  print LOG "\n\t\tToo many cycles required to synthesize this probe ($cycles)";
	}

	$cycles_failed = 1;
	next(LENGTH); #Note making the probe longer might help in this case, but we dont want to abandon any probes that were already successful
      }


      #Store this probe that passed basic quality check temporarily until the best one for this probe/length for this spot can be chosen
      $position_probe_count++;
      $position_probes_object{$position_probe_count}{sequence} = $probe_seq;
      $position_probes_object{$position_probe_count}{start_pos} = $probe_seq_start;
      $position_probes_object{$position_probe_count}{end_pos} = $probe_seq_end;
      $position_probes_object{$position_probe_count}{Tm} = $Tm_celsius;
      $position_probes_object{$position_probe_count}{masked_bases} = $n_count;
      $position_probes_object{$position_probe_count}{length} = $j;
      $position_probes_object{$position_probe_count}{cycles} = $cycles;

    }#LENGTH variance loop - altering probe length to achieve target Tm

    #Make sure at least some probes were found for this position by varying the length
    my $position_probes_found = keys %position_probes_object;
    unless ($position_probes_found >= 1){
      if ($verbose eq "yes"){
	print YELLOW, "\n\t\tNo suitable probes found for start position: $i", RESET;
	print LOG "\n\t\tNo suitable probes found for start position: $i";
      }

      #No probes were found (of any length) that passed - Summarize the types of failures
      #First failed because of ambiguous bases, then failed because repeat-masked bases, finally failed because of low complexity elements
      if ($ambigious_bases_failed == 1){
	$ambigious_bases_failed_count++;
      }
      if ($repeat_bases_failed == 1){
	$repeat_bases_failed_count++;
      }
      if ($low_complexity_bases_failed == 1){
	$low_complexity_bases_failed_count++;
      }
      if ($cycles_failed == 1){
	$cycles_failed_count++;
      }
      next(POS);
    }

    #Go through each probe for this position, calculate absolute Tm difference from target, and note the 'best' probe
    my $best_probe_id;
    my $closest_tm = 1000;  #Arbitrarily large tm diff
    foreach my $pp_count (keys %position_probes_object){
      my $tm_diff = abs($position_probes_object{$pp_count}{Tm} - $target_tm);
      if ($tm_diff < $closest_tm){
	$closest_tm = $tm_diff;
	$best_probe_id = $pp_count;
      }
    }

    #Note: Even after selecting the best probe for each position, we still want to eliminate probes that are outside the desired Tm range
    #This will eliminate probes where even altering the length was not sufficient to get the target Tm.
    my $Tm_lower_limit = ($target_tm - $tm_range);
    my $Tm_upper_limit = ($target_tm + $tm_range);
    unless ($position_probes_object{$best_probe_id}{Tm} >= $Tm_lower_limit &&  $position_probes_object{$best_probe_id}{Tm} <= $Tm_upper_limit){
      if ($verbose eq "yes"){
	print YELLOW, "\n\t\tBest probe at this position has Tm: $position_probes_object{$best_probe_id}{Tm} outside of range ($Tm_lower_limit - $Tm_upper_limit)", RESET;
	print LOG "\n\t\tBest probe at this position has Tm: $position_probes_object{$best_probe_id}{Tm} outside of range ($Tm_lower_limit - $Tm_upper_limit)";

      }
      $tm_failed_count++;
      next(POS);
    }

    #If the probe made it this far, add it to the list of successful probes for this exon region and proceed to the next position within the region
    #Only select one probe for each tiling position - the one with the length that results in a Tm closest to the target_Tm
    $probe_count++;
    $probe_object{$probe_count}{region_id} = $region_id;
    $probe_object{$probe_count}{sequence} = $position_probes_object{$best_probe_id}{sequence};
    $probe_object{$probe_count}{chromosome} = $chromosome;
    $probe_object{$probe_count}{start_pos} = $position_probes_object{$best_probe_id}{start_pos};
    $probe_object{$probe_count}{end_pos} = $position_probes_object{$best_probe_id}{end_pos};
    $probe_object{$probe_count}{Tm} = $position_probes_object{$best_probe_id}{Tm};
    $probe_object{$probe_count}{masked_bases} = $position_probes_object{$best_probe_id}{masked_bases};
    $probe_object{$probe_count}{length} = $position_probes_object{$best_probe_id}{length};
    $probe_object{$probe_count}{cycles} = $position_probes_object{$best_probe_id}{cycles};

    #print "\nProgbeSEQ: $probe_seq";

  }#POS (Position) Tiling, Overlap loop - shifting position within the target exon region

  return(\%probe_object);
}



################################################################################################
#Given a probe object containing probe sequences for an exon region, print to a probe file.
################################################################################################
sub printProbeInfo{
  my %args = @_;
  my $probe_object_ref = $args{'-probe_object'};

  #Info that needs to be printed
  #Probe_Count\tProbeSet_ID\tRegion_ID\tSequence\tProbe_length\tProbe_Tm\tProbe_Type\tChromosome\tStart\tEnd\tmasked_bases\tNimbleGenCycles

  #Global file handle = OUTFILE

  $current_probeset_id++;

  foreach my $probe (sort {$a<=>$b} keys %{$probe_object_ref}){
    $probe_number++;
    $current_probe_id++;

    my $Tm = $probe_object_ref->{$probe}->{Tm};
    my $Tm_rounded = sprintf("%.3f", $Tm);

   print OUTFILE "$current_probe_id\t$current_probeset_id\t$probe_object_ref->{$probe}->{region_id}\t$probe_object_ref->{$probe}->{sequence}\t$probe_object_ref->{$probe}->{length}\t$Tm_rounded\t$probe_object_ref->{$probe}->{chromosome}\t$probe_object_ref->{$probe}->{start_pos}\t$probe_object_ref->{$probe}->{end_pos}\t$probe_object_ref->{$probe}->{masked_bases}\t$probe_object_ref->{$probe}->{cycles}\n";

  }
  return();
}


####################################################################################################
#The following subroutine was provided by NimbleGen                                                #
#Copyright 2001-2002, NimbleGen Systems Inc. All rights reserved.                                  #
####################################################################################################
sub calculate_nimblegen_cycles {
  my %args = @_;
  my $probe = $args{'-sequence'};

  my $synthesis_order = "ACGT";
  my @order = split("",$synthesis_order);

  # grab the probe sequence from the subroutine call and clean it up
  $probe = uc $probe;
  $probe =~ tr/AGCT//c;
  my @probe = split("",$probe);

  my $cycles = 0;
  # we know we have to do for every bp in the oligo
  for my $position (0 .. $#probe) {
    # an 'N' means skip 4 cycles.
    if($probe[-1] eq 'N') {
      pop @probe;
      $cycles+=4;
      if(scalar(@probe) == 0) {
	return $cycles;
      }
      next;
    }
    foreach my $base (@order) {
      $cycles++;
      # examine the 3' bp. If it equals the bp we're currently
      # making, pop it off the stack. When we get to the end
      # return the current cycle count
      if ($probe[-1] eq $base) { pop @probe };
      if (scalar(@probe) == 0) {
	return ($cycles);
      }
    }
  }
}
