#!/usr/bin/perl -w
#Written by Malachi Griffith

#The purpose of this script is to generate a design submission file for a custom microarray.
#This file is essentially a list of probe sequences with unique IDs.
#A custom array manufacturer (NimbleGen, Agilent, etc.) would use this to create a physical microarray representing the design file

#STEPS

#1.) Import target exon regions from an input file (region id, chromosome, start, end)

#2.) Import probe records from an input file (probe_id, region_id, start, end)

#3.) Import random sequence negative control probe records.  These will be used to fill space.  The minimum number will be specified by the user

#4.) Go through each exon region, find the probes that overlap the region.  Note that some probes will not overlap any region
#      - For each probe that does overlap a region, calculate the distance between the centre of the region and the centre of the probe
#      - Note that it is possible that a probe will overlap more than one region (if the regions themselves are close or overlapping)
#      - Make sure each probe is associated with the region it is the closest to the centre of
#      - For probes that do not overlap a region, determine the region they are closest to and the distance from the centre of this region

#5.) Region probe selection
#    - Select the 'best' probe (according to distance from the centre of the exon region) that has not already been selected
#    - Mark this probe as selected
#    - Once this has been done for every region, start back on the first region
#    - Continue this process until (array_capacity - min_nc_probes) has been achieved
#    - Display a brief summary at the end of each pass

#6.) Negative control probe selection
#    - Fill all remaining space on the array with negative control probes.  Chose them to uniformly represent the Tm and length of region probes

#7.) Go through the input probe files and print out a new probe files consisting only of those probes selected
#    - Do this for both the Region probes and negative control probes

#8.) Print out the design file (probe ID, region ID, sequence)
#    - This file should be named to indicate the capacity of the target array and the source genome build
#    - eg. GTAC_V1_385k_hg18
#    - Format of this file is as follows
#    - Probe ID, target region ID, Probe sequence (eg. 1, 10730, ATAAAAAGTGAGTATCCAAAGCATATAATGTACCCCAGGTAAACAGCTTGTGTCACTTGAAGCT)

#9.) Create a summary of the design for all regions.  Each record will contain the complete region record but also the following information
#    - Number of filtered probes found, number of probes selected, status (no probes, partial probeset, full probeset)
#    - This file should be named to indicate the number of probes which were attempted for each region and the capacity of the target array
#    - eg. GTAC_V1_385k_hg18

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);

#Initialize command line options
my $target_regions = '';
my $xmr_pcr_regions = '';
my $region_probes = '';
my $nc_probes = '';
my $array_capacity = '';
my $min_nc_probes = '';
my $selected_probes_dir = '';
my $design_name = '';
my $ucsc_dir = '';
my $web_path = '';

GetOptions ('target_regions=s'=>\$target_regions, 'xmr_pcr_regions=s'=>\$xmr_pcr_regions,
	    'region_probes=s'=>\$region_probes, 'nc_probes=s'=>\$nc_probes,
	    'array_capacity=i'=>\$array_capacity, 'min_nc_probes=i'=>\$min_nc_probes,
	    'selected_probes_dir=s'=>\$selected_probes_dir, 'design_name=s'=>\$design_name,
	    'ucsc_dir=s'=>\$ucsc_dir, 'web_path=s'=>\$web_path);

print GREEN, "\n\nThis script generates an array design submission file for all regions in a specified in an input file\n\n", RESET;

print GREEN, "\nParameters:", RESET;
print GREEN, "\n\tSpecify the path to the target exon region file used to select probes for the array using: --target_regions", RESET;
print GREEN, "\n\tSpecify the path to the a file describing the regions targeted by PCR amplification using: --xmr_pcr_regions", RESET;
print GREEN, "\n\tSpecify the path to the region probes file using: --region_probes", RESET;
print GREEN, "\n\tSpecify the path to the negative control probes file using: --nc_probes",RESET;
print GREEN, "\n\tSpecify the capacity of the microarray you are creating a design for using: --array_capacity (e.g. 385000)", RESET;
print GREEN, "\n\tSpecify the minimum number of random negative control probes you want on the final array using: --min_nc_probes (e.g 1000)", RESET;
print GREEN, "\n\tSpecify the directory for selected probe files using: --selected_probes_dir", RESET;
print GREEN, "\n\tSpecify the desired design name using: --design_name", RESET;
print GREEN, "\n\tSpecify the target UCSC directory for custom UCSC track files using: --ucsc_dir", RESET;
print GREEN, "\n\tSpecify the html web path to this directory using: --web_path", RESET;
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tcreateDesignSubmission.pl  --target_regions=/projects/malachig/GTAC_Chip_Design_A/TargetExonRegions.txt  --xmr_pcr_regions=/projects/malachig/GTAC_Chip_Design_A/XMR_Amplicons.bed  --region_probes=/projects/malachig/GTAC_Chip_Design_A/probes/filtered/regionProbes_filtered.txt  --nc_probes=/projects/malachig/GTAC_Chip_Design_A/probes/filtered/negativeControlProbes_filtered.txt  --array_capacity=385000  --min_nc_probes=1000  --selected_probes_dir=/projects/malachig/GTAC_Chip_Design_A/probes/selected/  --design_name=GTAC_V1_385k_hg18  --ucsc_dir=/home/malachig/www/public/htdocs/GTAC_Chip_Design/  --web_path=http://www.bcgsc.ca/people/malachig/htdocs/GTAC_Chip_Design/\n\n", RESET;

unless ($target_regions && $xmr_pcr_regions && $region_probes && $nc_probes && ($array_capacity =~ /\d+/) && ($min_nc_probes =~ /\d+/) && $selected_probes_dir && $design_name && $ucsc_dir && $web_path){
  print RED, "\nRequired input parameter(s) missing or incorrect format!\n\n", RESET;
  exit();
}

my $log_file = "createDesignSubmission_"."$design_name"."_LOG.txt";
open (LOG, ">$log_file") || die "\nCould not open log file: $log_file!\n\n";

print LOG "\nUsed the following parameters for design submission:\ntarget_regions = $target_regions\nregion_probes = $region_probes\nnc_probes = $nc_probes\narray_capacity = $array_capacity\nmin_nc_probes = $min_nc_probes\nselected_probes_dir = $selected_probes_dir\ndesign_name = $design_name\n\n";


#1.) Import target regions from an input file (region id, chromosome, start, end)
my $input_header;
my $chromosomes_ref = &importTargetRegions('-region_file'=>$target_regions);

my $pcr_regions_ref = &importPcrRegions('-xmr_pcr_regions'=>$xmr_pcr_regions);

#2.) Import region probe records from an input file (probe_id, region_id, start, end)
#    - Which exon region each probe maps to will also be imported at this time
my %region_columns;
my $total_overlaping_probes = 0;
my $probes_ref = &importRegionProbes('-input_file'=>$region_probes);

#3.) Import random sequence negative control probe records.  These will be used to fill space.  The minimum number will be specified by the user
my %nc_columns;
my $nc_probes_ref = &importNCProbes('-input_file'=>$nc_probes);


#5.) Region probe selection
#    - Select the 'best' probe (according to distance from the centre of the exon region) that has not already been selected
#    - Mark this probe as selected
#    - Once this has been done for every region, start back on the first region
#    - Continue this process until (array_capacity - min_nc_probes) has been achieved
#    - Display a brief summary at the end of each pass
&selectRegionProbes('-regions_object'=>$chromosomes_ref);


#6.) Negative control probe selection
#    - Fill all remaining space on the array with negative control probes.  Chose them to uniformly represent the Tm and length of region probes
selectNegativeControlProbes('-nc_probe_count'=>$min_nc_probes, '-nc_probes_ref'=>$nc_probes_ref);

#7.) Go through the input probe files and print out a new probe file consisting only of those probes selected
#    - Do this for both the Region probes and negative control probes
#    - At the same time retrieve the probe sequences for only those probes selected for the array
unless ($selected_probes_dir =~ /\/$/){
  $selected_probes_dir = "$selected_probes_dir"."/";
}
my $selected_region_probes_file = "$selected_probes_dir"."$design_name"."_RegionProbes.txt";
my $selected_nc_probes_file = "$selected_probes_dir"."$design_name"."_randomControlProbes.txt";
my %selected_region_probes;
&printSelectedProbeFiles('-region_probes_in'=>$region_probes, '-nc_probes_in'=>$nc_probes,
			 '-region_probes_out'=>$selected_region_probes_file, '-nc_probes_out'=>$selected_nc_probes_file);


#8.) Print out the design file (probe ID, region ID, sequence)
#    - This file should be named to indicate the capacity of the target array and the source genome build
#    - eg. GTAC_V1_385k_hg18
#    - Format of this file is as follows
#    - Probe ID, target region ID, Probe sequence (eg. 1, 10730, ATAAAAAGTGAGTATCCAAAGCATATAATGTACCCCAGGTAAACAGCTTGTGTCACTTGAAGCT)

my $design_file = "$design_name"."_design.txt";
&printDesignFile('-design_file'=>$design_file);


#9.) Create a summary of the design for all regions.  Each record will contain the complete region record but also the following information
#    - Number of filtered probes found, number of probes selected, status (no probes, partial probeset, full probeset)
#    - This file should be named to indicate the number of probes which were attempted for each region and the capacity of the target array
#    - eg. GTAC_V1_385k_hg18
#    - At this time, also create a UCSC custom track for each chromosome and create a link for each exon region in the summary file

my $summary_file = "$design_name"."_summary.txt";
&printSummaryFile('-summary_file'=>$summary_file, '-ucsc_dir'=>$ucsc_dir, '-web_path'=>$web_path);

print "\n\n";

exit();



######################################################################################################################################
#1.) Import target regions from an input file (region id, chromosome, start, end)                                                    #
#- perform basic sanity checks on the input file and organize regions accorinding their target sequence (chromosome)                 #
######################################################################################################################################
sub importTargetRegions{
  my %args = @_;
  my $region_file = $args{'-region_file'};

  my %chromosomes;
  my %columns;

  print YELLOW, "\n\nImporting target regions from: $region_file", RESET;
  print LOG "\n\nImporting target regions from: $region_file";

  open (REGIONS, "$region_file") || die "\nCould not open target region file: $region_file\n\n";

  my $first_line = 1;
  my $record_count = 0;
  my $region_count = 0;

  while(<REGIONS>){
    chomp ($_);
    my @line = split ("\t", $_);

    if ($first_line == 1){
      $first_line = 0;
      $input_header = $_;

      my $col_count = 0;
      foreach my $col_name (@line){
	$columns{$col_name}{position} = $col_count;
	$col_count++;
      }

      next();
    }
    $record_count++;

    my $line = $_;
    my $region_id = $line[$columns{'ExonRegion_ID'}{position}];
    my $chromosome = $line[$columns{'Chromosome'}{position}];
    $chromosome = "chr"."$chromosome";
    my $start = $line[$columns{'Start'}{position}];
    my $end = $line[$columns{'End'}{position}];

    #Check data formats:
    unless ($region_id =~ /^\d+$/){print RED, "\nRegion ID format not correct: ($region_id) - check input file (line: $record_count + 1)!\n\n", RESET; exit();}
    unless ($chromosome =~ /^chr\d+$|^chr\w+$/){print RED, "\nChromosome name format not correct: ($chromosome) - check input file (line: $record_count + 1)!\n\n", RESET; exit();}
    unless ($start =~ /^\d+$/){print RED, "\nStart coordinate format not correct: ($start) - check input file (line: $record_count + 1)!\n\n", RESET; exit();}
    unless ($end =~ /^\d+$/){print RED, "\nEnd coordinate format not correct: ($end) - check input file (line: $record_count + 1)!\n\n", RESET; exit();}

    my $temp = $start + (($end - $start)/2);
    my $region_centre = sprintf("%.0f", $temp);

    #Sanity checks
    if ($start >= $end){
      print RED, "\nStart coordinate ($start) must be smaller than end coordinate ($end) - check input file!\n\n", RESET;
      exit();
    }

    #Organize the regions according to the chromosome they come from
    if ($chromosomes{$chromosome}){

      my $regions_ref = $chromosomes{$chromosome}{regions};

      if ($regions_ref->{$region_id}){
	print RED, "\nRegion ID: $region_id appears to be a duplicate - check input file!\n\n", RESET;
	exit();
      }
      $region_count++;

      $regions_ref->{$region_id}->{record_count} = $record_count;
      $regions_ref->{$region_id}->{chromosome} = $chromosome;
      $regions_ref->{$region_id}->{start} = $start;
      $regions_ref->{$region_id}->{end} = $end;
      $regions_ref->{$region_id}->{region_centre} = $region_centre;
      $regions_ref->{$region_id}->{region_size} = ($end - $start)+1;
      $regions_ref->{$region_id}->{line} = $line;
      $regions_ref->{$region_id}->{selected_probe_count} = 0;

      my %mapped_probe_list;
      $regions_ref->{$region_id}->{mapped_probe_list} = \%mapped_probe_list; #List of probe that have been mapped to this region.  Probes that are closest to or overlap this region

    }else{
      my %regions;
      $region_count++;
      $regions{$region_id}{record_count} = $record_count;
      $regions{$region_id}{chromosome} = $chromosome;
      $regions{$region_id}{start} = $start;
      $regions{$region_id}{end} = $end;
      $regions{$region_id}{region_centre} = $region_centre;
      $regions{$region_id}{region_size} = ($end - $start)+1;
      $regions{$region_id}{line} = $line;
      $regions{$region_id}{selected_probe_count} = 0;

      my %mapped_probe_list;
      $regions{$region_id}{mapped_probe_list} = \%mapped_probe_list; #List of probe that have been mapped to this region.  Probes that are closest to or overlap this region

      $chromosomes{$chromosome}{regions} = \%regions;
    }
  }
  close (REGIONS);



  print BLUE, "\n\n\tFound $region_count exon regions\n\n", RESET;
  print LOG "\n\n\tFound $region_count exon regions\n\n";

  return(\%chromosomes);
}


##########################################################################################################################################
#1b.) Import regions targeted by PCR amplification                                                                                       #
##########################################################################################################################################
sub importPcrRegions{
  my %args = @_;
  my $file = $args{'-xmr_pcr_regions'};

  my %pcr_regions;

  open (REGIONS, "$file") || die "\nCould not open PCR regions file: $file\n\n";

  my $region_count = 0;
  my $first_line = 1;

  while(<REGIONS>){
    chomp ($_);
    my @line = split ("\t", $_);

    if ($first_line == 1){
      $first_line = 0;
      next();
    }
    $region_count++;

    my $chromosome = $line[0];
    my $start = $line[1];
    my $end = $line[2];
    my $region_name = $line[3];

    if ($pcr_regions{$chromosome}){
      my $regions_ref = $pcr_regions{$chromosome}{regions};

      $regions_ref->{$region_count}->{start} = $start;
      $regions_ref->{$region_count}->{end} = $end;
      $regions_ref->{$region_count}->{name} = $region_name;

    }else{
      my %regions;
      $regions{$region_count}{start} = $start;
      $regions{$region_count}{end} = $end;
      $regions{$region_count}{name} = $region_name;

      $pcr_regions{$chromosome}{regions} = \%regions;
    }

  }

  close(REGIONS);

  return(\%pcr_regions);
}


##########################################################################################################################################
#2.) Import region probe records from an input file (probe_id, region_id, start, end)                                                    #
##########################################################################################################################################
sub importRegionProbes{
  my %args = @_;
  my $probe_file = $args{'-input_file'};

  my $progress_count = 0;
  my $blocks_imported = 0;
  my $probe_count = 0;

  my %probes;

  #Open the probe file and read the neccessary probe data into a hash keyed on probe ID
  open (PROBES, "$probe_file") || die "\nCould not open input probe file: $probe_file\n\n";


  my $first_line = 1;

  print BLUE, "\nImporting Region probe records from: $probe_file\n", RESET;
  print LOG "\nImporting Region probe records from: $probe_file\n";

  while (<PROBES>){
    $progress_count++;
    if ($progress_count == 10000){
      $blocks_imported++;
      print BLUE, "\n\tImported $blocks_imported blocks of 10,000 probes", RESET;
      print LOG "\n\tImported $blocks_imported blocks of 10,000 probes";
      $progress_count = 0;
    }
    chomp($_);
    my @line = split("\t", $_);

    #Get the header line and identify column names and their positions
    if ($first_line == 1){

      my $col_count = 0;
      foreach my $column (@line){
	$region_columns{$column}{column_pos} = $col_count;
	$col_count++;
      }
      $first_line = 0;

      #Check for critical columns and their names
      unless ($region_columns{Probe_Count} && $region_columns{Chromosome} && $region_columns{Unit1_start} && $region_columns{Unit1_end} && $region_columns{TargetExonRegion_ID} && $region_columns{DistanceFromExonRegion} && $region_columns{OverlapsExonRegion}){
	print RED, "\nCritical column missing or named incorrectly, check input file", RESET;
	exit();
      }
      next();
    }

    #Get the values of interest from each line (probe record)
    my $probe_id = $line[$region_columns{Probe_Count}{column_pos}];

    unless ($probe_id =~ /^\d+/){
      print RED, "\nInvalid probe or probeset ID\n\nLINE: $_\n\n", RESET;
      exit();
    }

    $probe_count++;

    #Required column headings:
    #Probe_Count, ProbeSet_ID, Region_ID, Probe_length, Probe_Tm, Unit1_start, Unit1_end, TargetExonRegion_ID, DistanceFromExonRegion, OverlapsExonRegion

    my $probeset_id = $line[$region_columns{ProbeSet_ID}{column_pos}];
    my $chromosome = $line[$region_columns{Chromosome}{column_pos}];
    my $unit1_start = $line[$region_columns{Unit1_start}{column_pos}];
    my $unit1_end = $line[$region_columns{Unit1_end}{column_pos}];
    my $exon_region_id = $line[$region_columns{TargetExonRegion_ID}{column_pos}];
    my $distance_from_exon = $line[$region_columns{DistanceFromExonRegion}{column_pos}];
    my $overlaps_exon = $line[$region_columns{OverlapsExonRegion}{column_pos}];

    if ($unit1_end < $unit1_start){
      print RED, "\nEnd coordinate of probe $probe_id is smaller than start coordinate!\n\n", RESET;
      exit();
    }

    $probes{$probe_id}{probeset_id} = $probeset_id;
    $probes{$probe_id}{chromosome} = $chromosome;
    $probes{$probe_id}{unit1_start} = $unit1_start;
    $probes{$probe_id}{unit1_end} = $unit1_end;
    $probes{$probe_id}{selected} = "no";
    $probes{$probe_id}{exon_region_id} = $exon_region_id;
    $probes{$probe_id}{distance_from_exon} = $distance_from_exon;
    $probes{$probe_id}{overlaps_exon} = $overlaps_exon;

    if ($overlaps_exon == 1){
      $total_overlaping_probes++;
    }

    #Add this probe to the list of for its target exon
    my $regions_ref = $chromosomes_ref->{$chromosome}->{regions};
    my $mapped_probe_list_ref = $regions_ref->{$exon_region_id}->{mapped_probe_list};
    $mapped_probe_list_ref->{$probe_id}->{distance_to_closest_region} = $distance_from_exon;

  }
  return(\%probes);
}


##################################################################################################################################################
#3.) Import random sequence negative control probe records.  These will be used to fill space.  The minimum number will be specified by the user #
##################################################################################################################################################
sub importNCProbes{
  my %args = @_;
  my $probe_file = $args{'-input_file'};

  my %nc_probes;
  my $progress_count = 0;
  my $blocks_imported = 0;
  my $probe_count = 0;

  #Open the probe file and read the neccessary probe data into a hash keyed on probe ID
  open (PROBES, "$probe_file") || die "\nCould not open input probe file: $probe_file\n\n";

  my $first_line = 1;

  print BLUE, "\n\nImporting NC probe records from: $probe_file\n", RESET;
  print LOG "\n\nImporting NC probe records from: $probe_file\n";

  while (<PROBES>){
    $progress_count++;
    if ($progress_count == 10000){
      $blocks_imported++;
      print BLUE, "\n\tImported $blocks_imported blocks of 10,000 probes", RESET;
      print LOG "\n\tImported $blocks_imported blocks of 10,000 probes";
      $progress_count = 0;
    }
    chomp($_);
    my @line = split("\t", $_);

    #Get the header line and identify column names and their positions
    if ($first_line == 1){

      my $col_count = 0;
      foreach my $column (@line){
	$nc_columns{$column}{column_pos} = $col_count;
	$col_count++;
      }
      $first_line = 0;

      #Check for critical columns and their names
      unless ($nc_columns{Probe_Count}){
	print RED, "\nCritical column missing or named incorrectly, check input file", RESET;
	exit();
      }
      next();
    }

    #Get the values of interest from each line (probe record)
    my $probe_id = $line[$nc_columns{Probe_Count}{column_pos}];

    unless ($probe_id =~ /^\d+/){
      print RED, "\nInvalid probe or probeset ID\n\nLINE: $_\n\n", RESET;
      exit();
    }
    $probe_count++;

    #Associate the required probe info with the target region involved

    #Required column headings:
    #Probe_Count, ProbeSet_ID, Probe_length, Probe_Tm
    my $probeset_id = $line[$nc_columns{ProbeSet_ID}{column_pos}];
    my $probe_length = $line[$nc_columns{Probe_length}{column_pos}];
    my $probe_tm = $line[$nc_columns{Probe_Tm}{column_pos}];

    $nc_probes{$probe_id}{probeset_id} = $probeset_id;
    $nc_probes{$probe_id}{probe_length} = $probe_length;
    $nc_probes{$probe_id}{probe_tm} = $probe_tm;
    $nc_probes{$probe_id}{selected} = "no";

  }
  return(\%nc_probes);
}


#######################################################################################################
#5.) Region probe selection                                                                           #
#######################################################################################################
sub selectRegionProbes{
  my %args = @_;
  my $regions_ref = $args{'-regions_object'};

  my $probes_needed = $array_capacity - $min_nc_probes;

  print BLUE, "\n\nAttempting to select ($array_capacity - $min_nc_probes) = $probes_needed region probes\n", RESET;
  print LOG "\n\nAttempting to select ($array_capacity - $min_nc_probes) = $probes_needed region probes\n";

  print BLUE, "\n\nFirst selecting all probes which actually overlap a target exon region ($total_overlaping_probes)\n", RESET;
  print LOG "\n\nFirst selecting all probes which actually overlap a target exon region ($total_overlaping_probes)\n";

  my $probes_selected = 0;
  my $pass_count = 0;

  if ($total_overlaping_probes > $probes_needed){
    print RED, "\nFound $total_overlaping_probes and only $probes_needed probe spots are available on the array!!\n\n", RESET;
    exit();
  }


  #First go through all probes and immediately select those that overlap the target region
  foreach my $probe_id (sort keys %{$probes_ref}){

    if ($probes_ref->{$probe_id}->{overlaps_exon} == 1){
      $probes_ref->{$probe_id}->{selected} = "yes";
      my $chr = $probes_ref->{$probe_id}->{chromosome};
      my $region_id = $probes_ref->{$probe_id}->{exon_region_id};

      my $regions_ref = $chromosomes_ref->{$chr}->{regions};
      $regions_ref->{$region_id}->{selected_probe_count}++;
      $probes_selected++;
    }
  }

  print BLUE, "\n\nSelected $probes_selected probes which actually overlap the target region\n", RESET;
  print BLUE, "\tNow proceeding to select additional probes to fill the array by selecting additional flanking probes for each region\n", RESET;

  print LOG "\n\nSelected $probes_selected probes which actually overlap the target region\n";
  print LOG "\tNow proceeding to select additional probes to fill the array by selecting additional flanking probes for each region\n";


  #Go through all target exon regions an arbitrarily large number of times (if there were only one target region, then $probes_needed passes would be required
 PASS:for (my $i = 1; $i < $probes_needed; $i++){

    #Go through each region and select the best probe per region
    #Make sure it has not already been selected and that it does not overlap any probes selected (unless overlaps were allowed by the user)

    my $pass_count++;

    print YELLOW, "\n\tSelecting region probes for pass $i (probes already selected: $probes_selected)", RESET;
    print LOG "\n\tSelecting region probes for pass $i (probes already selected: $probes_selected)";


    foreach my $chromosome (sort keys %{$chromosomes_ref}){

      my $regions_ref = $chromosomes_ref->{$chromosome}->{regions};


    REGION:foreach my $region_id (sort {$regions_ref->{$b}->{region_centre} <=> $regions_ref->{$a}->{region_centre}} keys %{$regions_ref}){

	#Skip regions with 0 probes
	unless ($regions_ref->{$region_id}->{mapped_probe_list}){
	  next REGION;
	}

	my $mapped_probe_list_ref = $regions_ref->{$region_id}->{mapped_probe_list};

	foreach my $probe_id (sort {$mapped_probe_list_ref->{$a}->{distance_to_closest_region} <=> $mapped_probe_list_ref->{$b}->{distance_to_closest_region}} keys %{$mapped_probe_list_ref}){

	  #If the required number of probes has already been achieved, stop looking
	  if ($probes_selected == $probes_needed){
	    last PASS;
	  }

	  #Make sure this probe has not already been selected
	  if ($probes_ref->{$probe_id}->{selected} eq "yes"){
	    next();
	  }else{

	    #Assign the status of this probe to selected in the master probe hash
	    $probes_ref->{$probe_id}->{selected} = "yes";

	    $regions_ref->{$region_id}->{selected_probe_count}++;

	    $probes_selected++;
	    next REGION;
	  }
	}
      }
    }
  }

  print BLUE, "\n\nSelected $probes_selected probes\n", RESET;
  print LOG "\n\nSelected $probes_selected probes\n";

  return();
}


#######################################################################################################
#6.) Negative control probe selection
#######################################################################################################
sub selectNegativeControlProbes{
  my %args = @_;
  my $target_nc_probe_count = $args{'-nc_probe_count'};
  my $nc_probes_ref = $args{'-nc_probes_ref'};

  print BLUE, "\n\nSelecting the desired number ($target_nc_probe_count) of negative control probes, drawing evenly from passing probes of each Tm/Length bin\n", RESET;
  print LOG "\n\nSelecting the desired number ($target_nc_probe_count) of negative control probes, drawing evenly from passing probes of each Tm/Length bin\n";

  #Organize all probes according to their probe length and Tm (to 0.1 decimal places)
  my %probe_bins;

  my $potential_candidates = 0;
  foreach my $probe_id (keys %{$nc_probes_ref}){

    $potential_candidates++;
    my $length = $nc_probes_ref->{$probe_id}->{probe_length};
    my $probe_tm = $nc_probes_ref->{$probe_id}->{probe_tm};
    my $tm = sprintf("%.1f", $probe_tm);

    if ($probe_bins{$length}){
      #Length has been seen before
      my $probe_length_bin_ref = $probe_bins{$length}{probe_length_bins};

      if ($probe_length_bin_ref->{$tm}){
	#Tm has been observed before
	my $probe_tm_bin_ref = $probe_length_bin_ref->{$tm}->{tm_probes};
	$probe_tm_bin_ref->{$probe_id}->{tm} = $probe_tm;

      }else{
	#First probe of this Tm
	my %probe_tm_bin;
	$probe_tm_bin{$probe_id}{tm} = $probe_tm;
	$probe_length_bin_ref->{$tm}->{tm_probes} = \%probe_tm_bin;
      }
    }else{
      #First probe of this length and Tm
      my %probe_tm_bin;
      $probe_tm_bin{$probe_id}{tm} = $probe_tm;

      my %probe_length_bin;
      $probe_length_bin{$tm}{tm_probes} = \%probe_tm_bin;

      $probe_bins{$length}{probe_length_bins} = \%probe_length_bin;
    }
  }

  #Make sure enough passing NC probes were found to accomodate the desired target number
  unless ($potential_candidates > $target_nc_probe_count){
    print RED, "\nNot enough passing negative control probes ($potential_candidates) were found to accomodate the desired target: $target_nc_probe_count", RESET;
    exit();
  }

  #Now go through each length-bin -> tm-bin -> probe_list  and select one probe from each list.
  #Repeat entire process until desired number are chosen
  my %final_nc_probes;
  my $final_nc_probe_count = 0;
  my $nc_probes_selected = 0;

  my $iter = 0;
  #Until desired probe count is reached
  ITER:while ($final_nc_probe_count < $target_nc_probe_count){
      $iter++;

      my $length_bin_count = keys %probe_bins;
      print YELLOW, "\n\tSelecting negative control probes for pass $iter (probes already selected: $nc_probes_selected)", RESET;
      print LOG "\n\tSelecting negative control probes for pass $iter (probes already selected: $nc_probes_selected)";

      #Go through each length bin
      LENGTH:foreach my $length (sort {$a <=> $b} keys %probe_bins){

	  my $tm_bin_ref = $probe_bins{$length}{probe_length_bins};
	  my $bin_count = keys %{$tm_bin_ref};

	  #Go through each Tm bin of this length
	TM:foreach my $tm (sort {$a <=> $b} keys %{$tm_bin_ref}){

	    my $tm_probes_ref = $tm_bin_ref->{$tm}->{tm_probes};

	    #Select 1 probe from this bin (make sure it hasnt already been added)
	    foreach my $probe_id (sort {$tm_probes_ref->{$a}->{tm} <=> $tm_probes_ref->{$b}->{tm}} keys %{$tm_probes_ref}){

	      unless($final_nc_probes{$probe_id}){
		$final_nc_probe_count++;
		$final_nc_probes{$probe_id}{tm_celsius} = $tm_probes_ref->{$probe_id}->{tm};
		$final_nc_probes{$probe_id}{nc_probe_count} = $final_nc_probe_count;

		#Mark this probe as selected
		$nc_probes_ref->{$probe_id}->{selected} = "yes";
		$nc_probes_selected++;

		#See if the desired number of probes has been reached
		if ($final_nc_probe_count == $target_nc_probe_count){
		  last LENGTH;
		}else{
		  #otherwise just go to the next Tm bin of probes
		  next TM;
		}
	      }
	    }
	  }
	}
    }

  my $success_count = keys %final_nc_probes;
  print BLUE, "\n\nSelected $success_count probes\n", RESET;
  print LOG "\n\nSelected $success_count probes\n";

  return();
}


################################################################################################################
#7.) Go through the input probe files and print out a new probe files consisting only of those probes selected #
#    - Do this for both the Region probes and negative control probes                                          #
################################################################################################################
sub printSelectedProbeFiles{
  my %args = @_;
  my $region_in = $args{'-region_probes_in'};
  my $nc_in = $args{'-nc_probes_in'};
  my $region_out = $args{'-region_probes_out'};
  my $nc_out = $args{'-nc_probes_out'};

  print BLUE, "\n\nPrinting selected probe files for region and negative control probes", RESET;
  print BLUE, "\n\t$region_out\n\t$nc_out\n", RESET;
  print LOG "\n\nPrinting selected probe files for region and negative control probes";
  print LOG "\n\t$region_out\n\t$nc_out\n";

  #First go through the regions object and get a final list of selected probes
  foreach my $probe_id (sort keys %{$probes_ref}){
    if ($probes_ref->{$probe_id}->{selected} eq "yes"){
      $selected_region_probes{$probe_id}{printed} = "no";
    }
  }

  #Open the region probe input and output files - check each probe line in the infile, if the probe was selected print it to the outfile
  open (REG_IN, "$region_in") || die "\nCould not open region probe input file: $region_in\n\n";
  open (REG_OUT, ">$region_out") || die "\nCould not open region probe output file: $region_out\n\n";

  while (<REG_IN>){
    chomp($_);
    my $line = $_;
    my @line = split("\t", $line);

    if ($line[0] =~ /\d+/){
      my $probe_id = $line[0];
      if ($selected_region_probes{$probe_id}){
	my $sequence = $line[$region_columns{Sequence}{column_pos}];
	$selected_region_probes{$probe_id}{sequence} = $sequence;
	$selected_region_probes{$probe_id}{printed} = "yes";

	print REG_OUT "$line\n";
      }
    }else{
      print REG_OUT "$line\n";
    }
  }

  close (REG_IN);
  close (REG_OUT);

  #Open the NC probe input and output files - check each probe line in the infile, if the probe was selected print it to the outfile
  open (NEG_IN, "$nc_in") || die "\nCould not open NC input file: $nc_in\n\n";
  open (NEG_OUT, ">$nc_out") || die "\nCould not open NC output file: $nc_out\n\n";

  while (<NEG_IN>){
    chomp($_);
    my $line = $_;
    my @line = split("\t", $line);

    if ($line[0] =~ /\d+/){
      my $probe_id = $line[0];
      if ($nc_probes_ref->{$probe_id}->{selected} eq "yes"){
	my $sequence = $line[$nc_columns{Sequence}{column_pos}];
	$nc_probes_ref->{$probe_id}->{sequence} = $sequence;
	print NEG_OUT "$line\n";
      }
    }else{
      print NEG_OUT "$line\n";
    }
  }

  close (NEG_IN);
  close (NEG_OUT);

  return();
}


################################################################################################################
#8.) Print out the design file (probe ID, region ID, sequence)                                                 #
################################################################################################################
sub printDesignFile{
  my %args = @_;
  my $design_file = $args{'-design_file'};

  print BLUE, "\n\nPrinting design file: $design_file\n", RESET;
  print LOG "\n\nPrinting design file: $design_file\n";

  open (DESIGN, ">$design_file") || die "\nCould not open design file: $design_file\n\n", RESET;

  print DESIGN "Probe_ID\tRegion_ID\tSequence\n";

  #First print out the selected region probes
  foreach my $probe_id (sort {$a <=> $b} keys %{$probes_ref}){
    if ($selected_region_probes{$probe_id}){
      print DESIGN "$probe_id\t$probes_ref->{$probe_id}->{exon_region_id}\t$selected_region_probes{$probe_id}{sequence}\n";
    }
  }

  #Next print out the selected NC probes
  foreach my $probe_id (sort {$a <=> $b} keys %{$nc_probes_ref}){

    if ($nc_probes_ref->{$probe_id}->{selected} eq "yes"){
      print DESIGN "$probe_id\tna\t$nc_probes_ref->{$probe_id}->{sequence}\n";
    }
  }
  close (DESIGN);
  return();
}

##############################################################################################################################################
#9.) Create a summary of the design for all regions.  Each record will contain the complete region record but also the following information #
#    - Number of filtered probes found, number of probes selected, status (no probes, partial probeset, full probeset)                       #
##############################################################################################################################################
sub printSummaryFile{
  my %args = @_;
  my $summary_file = $args{'-summary_file'};
  my $ucsc_dir = $args{'-ucsc_dir'};
  my $web_path = $args{'-web_path'};

  unless ($ucsc_dir =~ /.*\/$/){
    $ucsc_dir = "$ucsc_dir"."/";
  }
  unless ($web_path =~ /.*\/$/){
    $web_path = "$web_path"."/";
  }


  print BLUE, "\n\nPrinting summary file: $summary_file\n", RESET;
  print LOG "\n\nPrinting summary file: $summary_file\n";

  open (SUMMARY, ">$summary_file") || die "\nCould not open summary file: $summary_file\n\n", RESET;

  print SUMMARY "$input_header\tProbe_Count\n";

  foreach my $chromosome (sort keys %{$chromosomes_ref}){

    my $ucsc_file = "$ucsc_dir"."$chromosome".".txt";

    print BLUE, "\nPrinting UCSC file for $chromosome: $ucsc_file", RESET;

    #Create ucsc track records for all PCR targeted regions on this chromosomes
    my $xmr_ref = $pcr_regions_ref->{$chromosome}->{regions};
    my @xmr_records;

    foreach my $pcr_region (sort {$xmr_ref->{$a}->{start} <=> $xmr_ref->{$b}->{start}} keys %{$xmr_ref}){
      my $xmr_record = "\n$chromosome\tDesign\tPcrRegion\t$xmr_ref->{$pcr_region}->{start}\t$xmr_ref->{$pcr_region}->{end}\t.\t.\t.\t$xmr_ref->{$pcr_region}->{name}";
      push (@xmr_records, $xmr_record);
    }


    #Create ucsc track records for all target regions and probes
    my @region_records;
    my @probe_records;
    my $regions_ref = $chromosomes_ref->{$chromosome}->{regions};

    foreach my $region_id (sort {$regions_ref->{$a}->{region_centre} <=> $regions_ref->{$b}->{region_centre}} keys %{$regions_ref}){

      my $region_record = "\n$chromosome\tDesign\tTargetRegion\t$regions_ref->{$region_id}->{start}\t$regions_ref->{$region_id}->{end}\t.\t.\t.\t$region_id";
      push(@region_records, $region_record);

      #Go through each probe associated with this region and print it to the UCSC file
      my $regions_ref = $chromosomes_ref->{$chromosome}->{regions};
      my $mapped_probe_list_ref = $regions_ref->{$region_id}->{mapped_probe_list};

      foreach my $probe_id (sort {$mapped_probe_list_ref->{$a}->{distance_to_closest_region} <=> $mapped_probe_list_ref->{$b}->{distance_to_closest_region}} keys %{$mapped_probe_list_ref}){

	if ($selected_region_probes{$probe_id}){
	  my $probe_record = "\n$chromosome\tDesign\tProbe\t$probes_ref->{$probe_id}->{unit1_start}\t$probes_ref->{$probe_id}->{unit1_end}\t.\t.\t.\t$probe_id";
	  push (@probe_records, $probe_record);
	}
      }

    #Example link
    #http://genome.ucsc.edu/cgi-bin/hgTracks?db=ce4&position=chrV:10384148-10385484&hgt.customText=http://www.bcgsc.ca/people/malachig/htdocs/alexa_tracks/ALEXA_ce_45_170b/chrV_25.txt.gz&ctfile_ce4=

      my $link = "http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg18&position=$chromosome:$regions_ref->{$region_id}->{start}-$regions_ref->{$region_id}->{end}&hgt.customText=$web_path"."$chromosome".".txt.gz"."&ctfile_hg18=";

      #Print a summary line for this region including a link to the appropriate track file
      print SUMMARY "$regions_ref->{$region_id}->{line}\t$regions_ref->{$region_id}->{selected_probe_count}\t$link\n";


    }

    #Print out all the UCSC track records gathered above for this chromosome
    open (UCSC, ">$ucsc_file") || die "\nCould not open ucsc file: $ucsc_file\n\n";

    #Browser line
    print UCSC "#Browser line";
    print UCSC "\nbrowser hide all";
    print UCSC "\nbrowser full knownGene";
    print UCSC "\nbrowser pack multiz28way";

    #1.) Track line for PCR targeted regions
    print UCSC "\n\n#XMR amplicon regions line";
    my $xmr_track_name = "XMR_Amplicons";
    #RED 204,0,0
    print UCSC "\ntrack name=$xmr_track_name description=\"XMR_Amplicon bounds for Exon pipeline\" color=204,0,0 useScore=0 visibility=3";
    print UCSC "\n\n#Begin DATA";

    foreach my $xmr_record (@xmr_records){
      print UCSC "$xmr_record";
    }

    #2.) Track line for regions track
    print UCSC "\n\n#Track line";
    my $region_track_name = "$design_name"."_Regions";
    #GREEN 51,153,0
    print UCSC "\ntrack name=$region_track_name description=\"$design_name Target Regions\" color=51,153,0 useScore=0 visibility=3";
    print UCSC "\n\n#Begin DATA";

    foreach my $region_record (@region_records){
      print UCSC "$region_record";
    }

    #3.) Track line for probes track
    print UCSC "\n\n#Track line";
    my $probe_track_name = "$design_name"."_Probes";
    #BLUE: 0,0,255
    print UCSC "\ntrack name=$design_name description=\"$design_name Probe Positions\" color=0,0,255 useScore=0 visibility=4";
    print UCSC "\n\n#Begin DATA";

    foreach my $probe_record (@probe_records){
      print UCSC "$probe_record";
    }

    close(UCSC);

    #Now gzip this file
    my $cmd = "/usr/bin/gzip $ucsc_file";
    print BLUE, "\nExecuting: $cmd", RESET;
    system($cmd);

  }
  close(SUMMARY);

  return();
}


