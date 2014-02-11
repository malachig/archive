#!/usr/bin/perl -w
#Written by Malachi Griffith

#The purpose of this script is to generate a design submission file for a custom microarray.
#This file is essentially a list of probe sequences with unique IDs.
#A custom array manufacturer (NimbleGen, Agilent, etc.) would use this to create a physical microarray representing the design file

#STEPS

#1.) Import target regions from an input file (region id, chromosome, start, end)

#2.) Import region probe records from an input file (probe_id, region_id, start, end)
#    - Create a hash for all of the probes associated with each region ID and associate this hash with the corresponding region record

#3.) Import random sequence negative control probe records.  These will be used to fill space.  The minimum number will be specified by the user

#4.) Region probe selection
#    - Go through each region, get the probes for that region
#    - Select the 'best' probe (according to Tm & length) that hasnt been selected already and has NO OVERLAP with previously selected probes
#    - Mark this probe as selected
#    - Once this has been done for every region, start back on the first region
#    - Continue this process until (array_capacity - min_nc_probes) has been achieved
#    - Display a brief summary at the end of each pass

#5.) Negative control probe selection
#    - Fill all remaining space on the array with negative control probes.  Chose them to uniformly represent the Tm and length of region probes

#6.) Go through the input probe files and print out a new probe files consisting only of those probes selected
#    - Do this for both the Region probes and negative control probes

#7.) Print out the design file (probe ID, region ID, sequence)
#    - This file should be named to indicate the capacity of the target array and the source genome build
#    - eg. MR_CNV_V1_72k_hg18
#    - Format of this file is as follows
#    - Probe ID, target region ID, Probe sequence (eg. 1, 10730, ATAAAAAGTGAGTATCCAAAGCATATAATGTACCCCAGGTAAACAGCTTGTGTCACTTGAAGCT)

#8.) Create a summary of the design for all regions.  Each record will contain the complete region record but also the following information
#    - Number of filtered probes found, number of probes selected, status (no probes, partial probeset, full probeset)
#    - This file should be named to indicate the number of probes which were attempted for each region and the capacity of the target array
#    - eg. MR_CNV_V1_72k_hg18

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);

#Initialize command line options
my $target_regions = '';
my $region_probes = '';
my $nc_probes = '';
my $target_tm = '';
my $target_length = '';
my $allow_overlap = '';
my $array_capacity = '';
my $min_nc_probes = '';
my $region_probeset_threshold = '';
my $region_size_threshold = '';
my $selected_probes_dir = '';
my $design_name = '';

GetOptions ('target_regions=s'=>\$target_regions, 'region_probes=s'=>\$region_probes, 'nc_probes=s'=>\$nc_probes,
	    'target_tm=f'=>\$target_tm, 'target_length=i'=>\$target_length, 'allow_overlap=s'=>\$allow_overlap,
	    'array_capacity=i'=>\$array_capacity, 'min_nc_probes=i'=>\$min_nc_probes,
	    'region_probeset_threshold=i'=>\$region_probeset_threshold, 'region_size_threshold=i'=>\$region_size_threshold,
	    'selected_probes_dir=s'=>\$selected_probes_dir, 'design_name=s'=>\$design_name);

print GREEN, "\n\nThis script generates an array design submission file for all regions in a specified in an input file\n\n", RESET;

print GREEN, "\nParameters:", RESET;
print GREEN, "\n\tSpecify the path to the target region file using: --target_regions", RESET;
print GREEN, "\n\tSpecify the path to the region probes file using: --region_probes", RESET;
print GREEN, "\n\tSpecify the path to the negative control probes file using: --nc_probes",RESET;
print GREEN, "\n\tSpecify the target probe Tm used for this design using: --target_tm", RESET;
print GREEN, "\n\tSpecify the target probe length used for this design using: --target_length", RESET;
print GREEN, "\n\tSpecify whether you want to allow probes to overlap using: --allow_overlap (yes/no)", RESET;
print GREEN, "\n\tSpecify the capacity of the microarray you are creating a design for using: --array_capacity (e.g. 385000)", RESET;
print GREEN, "\n\tSpecify the minimum number of negative control probes you want on the final array using: --min_nc_probes (e.g 3850)", RESET;
print GREEN, "\n\tSpecify the desired number of probes per region using: --region_probeset_threshold and --region_size_threshold", RESET;
print GREEN, "\n\t\tOnce a region has --region_probeset_threshold probes, more probe will only be added if the region is larger than --region_size_threshold", RESET;
print GREEN, "\n\t\tSelection of probes within regions SMALLER than --region_size_threshold will be based primarily of TM", RESET;
print GREEN, "\n\t\tSelection of probes within regions LARGER than --region_size_threshold will also try to spread them evenly across the region", RESET;

print GREEN, "\n\tSpecify the directory for selected probe files using: --selected_probes_dir", RESET;
print GREEN, "\n\tSpecify the desired design name using: --design_name", RESET;
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tcreateDesignSubmission.pl  --target_regions=MasterRegionList-hg18-Malachi_FinalAfterMerging_V8.txt  --region_probes=RegionProbes_filtered.txt  --nc_probes=randomControlProbes_filtered.txt  --target_tm=76.0  --target_length=54  --allow_overlap=no  --array_capacity=72000  --min_nc_probes=720  --region_probeset_threshold=6  --region_size_threshold=1300  --selected_probes_dir=/home/malachig/MR_Chip_Design/probes/selected_probes  --design_name=MR_CNV_V1_72k_hg18\n\n", RESET;

unless ($target_regions && $region_probes && $nc_probes && ($target_tm =~ /\d+/) && ($target_length =~ /\d+/) && ($allow_overlap =~ /^yes$|^no$/i) && ($array_capacity =~ /\d+/) && ($min_nc_probes =~ /\d+/) && ($region_probeset_threshold =~ /\d+/) && ($region_size_threshold =~ /\d+/) && $selected_probes_dir && $design_name){
  print RED, "\nRequired input parameter(s) missing or incorrect format!\n\n", RESET;
  exit();
}

my $log_file = "createDesignSubmission_"."$design_name"."_LOG.txt";
open (LOG, ">$log_file") || die "\nCould not open log file: $log_file!\n\n";

print LOG "\nUsed the following parameters for design submission:\ntarget_regions = $target_regions\nregion_probes = $region_probes\nnc_probes = $nc_probes\ntarget_tm = $target_tm\ntarget_length = $target_length\nallow_overlap = $allow_overlap\narray_capacity = $array_capacity\nmin_nc_probes = $min_nc_probes\nregion_probeset_threshold = $region_probeset_threshold\nregion_size_threshold = $region_size_threshold\nselected_probes_dir = $selected_probes_dir\ndesign_name = $design_name\n\n";


#1.) Import target regions from an input file (region id, chromosome, start, end)
my $input_header;
my $regions_ref = &importTargetRegions('-region_file'=>$target_regions);

#2.) Import region probe records from an input file (probe_id, region_id, start, end)
#    - Create a hash for all of the probes associated with each region ID and associate this hash with the corresponding region record
my %region_columns;
&importRegionProbes('-input_file'=>$region_probes, '-regions_object'=>$regions_ref);

#Deal with special case regions by importing a second set of probes
my $extra_probes_file = "/home/malachig/MR_Chip_Design/probes/filtered_probes_poorSpecificity/RegionProbes_filtered.txt";
&importRegionProbes('-input_file'=>$extra_probes_file, '-regions_object'=>$regions_ref);

#3.) Import random sequence negative control probe records.  These will be used to fill space.  The minimum number will be specified by the user
my %nc_columns;
my $nc_probes_ref = &importNCProbes('-input_file'=>$nc_probes);

#4.) Region probe selection
#    - Go through each region, get the probes for that region
#    - Select the 'best' probe (according to Tm & length) that hasnt been selected already and has NO OVERLAP with previously selected probes
#    - Mark this probe as selected
#    - Once this has been done for every region, start back on the first region
#    - Continue this process until (array_capacity - min_nc_probes) has been achieved
#    - Display a brief summary at the end of each pass
&selectRegionProbes('-regions_object'=>$regions_ref);


#5.) Negative control probe selection
#    - Fill all remaining space on the array with negative control probes.  Chose them to uniformly represent the Tm and length of region probes
selectNegativeControlProbes('-nc_probe_count'=>$min_nc_probes, '-nc_probes_ref'=>$nc_probes_ref);

#6.) Go through the input probe files and print out a new probe file consisting only of those probes selected
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


#7.) Print out the design file (probe ID, region ID, sequence)
#    - This file should be named to indicate the capacity of the target array and the source genome build
#    - eg. MR_CNV_V1_72k_hg18
#    - Format of this file is as follows
#    - Probe ID, target region ID, Probe sequence (eg. 1, 10730, ATAAAAAGTGAGTATCCAAAGCATATAATGTACCCCAGGTAAACAGCTTGTGTCACTTGAAGCT)
my $design_file = "$design_name"."_design.txt";
&printDesignFile('-design_file'=>$design_file);


#8.) Create a summary of the design for all regions.  Each record will contain the complete region record but also the following information
#    - Number of filtered probes found, number of probes selected, status (no probes, partial probeset, full probeset)
#    - This file should be named to indicate the number of probes which were attempted for each region and the capacity of the target array
#    - eg. MR_CNV_V1_72k_hg18
my $summary_file = "$design_name"."_summary.txt";
&printSummaryFile('-summary_file'=>$summary_file);


print "\n\n";

exit();



######################################################################################################################################
#1.) Import target regions from an input file (region id, chromosome, start, end)                                                    #
#- perform basic sanity checks on the input file and organize regions accorinding their target sequence (chromosome)                 #
######################################################################################################################################
sub importTargetRegions{
  my %args = @_;
  my $region_file = $args{'-region_file'};

  my %regions;
  my %columns;

  print YELLOW, "\n\nImporting target regions from: $region_file", RESET;
  print LOG "\n\nImporting target regions from: $region_file";

  open (REGIONS, "$region_file") || die "\nCould not open target region file: $region_file\n\n";

  my $first_line = 1;
  my $record_count = 0;

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
    my $region_id = $line[$columns{'Region ID'}{position}];
    my $chromosome = $line[$columns{'Chromosome'}{position}];
    my $start = $line[$columns{'Start position (with flank)'}{position}];
    my $end = $line[$columns{'End position (with flank)'}{position}];

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
    $regions{$region_id}{selected_probe_count} = 0;

  }
  close (REGIONS);

  my $region_count = keys %regions;

  print BLUE, "\n\n\tFound $region_count regions\n\n", RESET;
  print LOG "\n\n\tFound $region_count regions\n\n";

  return(\%regions);
}


##########################################################################################################################################
#2.) Import region probe records from an input file (probe_id, region_id, start, end)                                                    #
#    - Create a hash for all of the probes associated with each region ID and associate this hash with the corresponding region record   #
##########################################################################################################################################
sub importRegionProbes{
  my %args = @_;
  my $probe_file = $args{'-input_file'};
  my $regions_ref = $args{'-regions_object'};

  my $progress_count = 0;
  my $blocks_imported = 0;
  my $probe_count = 0;

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
      unless ($region_columns{Probe_Count} && $region_columns{Region_ID}){
	print RED, "\nCritical column missing or named incorrectly, check input file", RESET;
	exit();
      }
      next();
    }

    #Get the values of interest from each line (probe record)
    my $probe_id = $line[$region_columns{Probe_Count}{column_pos}];
    my $region_id = $line[$region_columns{Region_ID}{column_pos}];

    unless ($probe_id =~ /^\d+/ && $region_id =~ /^\d+/){
      print RED, "\nInvalid probe or probeset ID\n\nLINE: $_\n\n", RESET;
      exit();
    }

    #Make sure the region ID in the probe file was found in the target region file!
    unless ($regions_ref->{$region_id}){
      print RED, "\nRegion ID in region probe file was not found in the target region file!\n\n", RESET;
      exit();
    }

    $probe_count++;

    #Associate the required probe info with the target region involved

    #Required column headings:
    #Probe_Count, ProbeSet_ID, Region_ID, Probe_length, Probe_Tm, Unit1_start, Unit1_end

    my $probeset_id = $line[$region_columns{ProbeSet_ID}{column_pos}];
    my $probe_length = $line[$region_columns{Probe_length}{column_pos}];
    my $probe_tm = $line[$region_columns{Probe_Tm}{column_pos}];
    my $unit1_start = $line[$region_columns{Unit1_start}{column_pos}];
    my $unit1_end = $line[$region_columns{Unit1_end}{column_pos}];

    if ($unit1_end < $unit1_start){
      print RED, "\nEnd coordinate of probe $probe_id is smaller than start coordinate!\n\n", RESET;
      exit();
    }

    #Calculate a quality score as follows so that Tm gives up to a value of 6 times the Tm range of probes and length is taken at par
    #This means that if region probes were allowed to vary by +/- 4.5 degrees and +/- 10 bp in length the total score possible is 27 + 10 = 37
    #A perfect score of 0 would occur for probes with the target length and tm.  Slightly more weight is given to the influence of Tm than length
    #This should be okay considering that the Tm of region probes is only allowed to vary by 4.5 degrees anyway.
    my $tm_diff = abs($probe_tm - $target_tm);
    my $length_diff = abs($probe_length - $target_length);
    my $quality_score = ($tm_diff * 6) + $length_diff;


    if ($regions_ref->{$region_id}->{probes}){
      my $probes_ref = $regions_ref->{$region_id}->{probes};
      $probes_ref->{$probe_id}->{probeset_id} = $probeset_id;
      $probes_ref->{$probe_id}->{probe_length} = $probe_length;
      $probes_ref->{$probe_id}->{probe_tm} = $probe_tm;
      $probes_ref->{$probe_id}->{unit1_start} = $unit1_start;
      $probes_ref->{$probe_id}->{unit1_end} = $unit1_end;
      $probes_ref->{$probe_id}->{quality_score} = $quality_score;
      $probes_ref->{$probe_id}->{selected} = "no";

    }else{
      my %probes;
      $probes{$probe_id}{probeset_id} = $probeset_id;
      $probes{$probe_id}{probe_length} = $probe_length;
      $probes{$probe_id}{probe_tm} = $probe_tm;
      $probes{$probe_id}{unit1_start} = $unit1_start;
      $probes{$probe_id}{unit1_end} = $unit1_end;
      $probes{$probe_id}{quality_score} = $quality_score;
      $probes{$probe_id}{selected} = "no";
      $regions_ref->{$region_id}->{probes} = \%probes;
    }
  }
  return();
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
#4.) Region probe selection                                                                           #
#######################################################################################################
sub selectRegionProbes{
  my %args = @_;
  my $regions_ref = $args{'-regions_object'};

  my $probes_needed = $array_capacity - $min_nc_probes;

  print BLUE, "\n\nAttempting to select ($array_capacity - $min_nc_probes) = $probes_needed region probes\n", RESET;
  print LOG "\n\nAttempting to select ($array_capacity - $min_nc_probes) = $probes_needed region probes\n";

  my $probes_selected = 0;
  my $pass_count = 0;

 PASS:for (my $i = 1; $i < 1000; $i++){

    #Go through each region and select the best probe per region
    #Make sure it has not already been selected and that it does not overlap any probes selected (unless overlaps were allowed by the user)

    my $pass_count++;

    print YELLOW, "\n\tSelecting region probes for pass $i (probes already selected: $probes_selected)", RESET;
    print LOG "\n\tSelecting region probes for pass $i (probes already selected: $probes_selected)";

  REGION:foreach my $region_id (sort {$regions_ref->{$b}->{region_size} <=> $regions_ref->{$a}->{region_size}} keys %{$regions_ref}){

      #Skip regions with 0 probes
      unless ($regions_ref->{$region_id}->{probes}){
	next REGION;
      }

      #If a region already has the desired number of probes, skip it unless it is larger than the specified cutoff
      if (($regions_ref->{$region_id}->{selected_probe_count} >= $region_probeset_threshold) && ($regions_ref->{$region_id}->{region_size} < $region_size_threshold)){
	next REGION;
      }

      my $probes_ref = $regions_ref->{$region_id}->{probes};

      #A.)Go through each probe and calculate an adjusted quality score which takes into account it's distance from already selected probes
      #   - If no probes have been selected yet, then the adjusted quality score is the same as the original quality score (based on Tm/length only)
      if ($regions_ref->{$region_id}->{selected_probes}){

	my $selected_probes_ref = $regions_ref->{$region_id}->{selected_probes};
	my $selected_probe_count = keys %{$selected_probes_ref};

      PROBE:foreach my $probe_id (keys %{$probes_ref}){
	  #Recalculate the quality score of this probe taking into account the minimum percent distance to any probe selected so far

	  my $min_distance = 1000000000; #start with arbitrarily large distance
	  my $region_size = $regions_ref->{$region_id}->{region_size};
	  my $probe_start = $probes_ref->{$probe_id}->{unit1_start};

	  foreach my $selected_probe_id (keys %{$selected_probes_ref}){

	    #Avoid self-self comparisons - Also avoid cases where the current probe is the ONLY selected probe
	    if ($selected_probe_id == $probe_id && $selected_probe_count == 1){
	      $probes_ref->{$probe_id}->{adjusted_quality_score} = $probes_ref->{$probe_id}->{quality_score};
	      next PROBE;
	    }

	    if ($selected_probe_id == $probe_id){
	      next();
	    }

	    my $selected_probe_start = $probes_ref->{$selected_probe_id}->{unit1_start};
	    my $distance = abs($probe_start - $selected_probe_start);
	    if ($distance < $min_distance){
	      $min_distance = $distance;
	    }
	  }

	  #Use the min distance between this probe and selected probes to calculate an adjusted score
	  if ($min_distance > $region_size){
	    print RED, "\nRegion: $region_id\tProbe: $probe_id\tMin distance ($min_distance) should not be larger than region size ($region_size) ...\n\n", RESET;
	    exit();
	  }

	  #Calculate adjusted score:
	  #Current score should have values 27 (Tm) + 10 (length) = 37 (where 0 would be a perfect score)
	  #Percent min distance will give a value of 0 - 100 (The higher the number the better!)
	  #Make sure negative or zero values are avoided here
	  #A.) Divide this number by 4 to reduce the weight of distance to other probes in the final adjusted score (values 0 - 25) - for large regions
	  #B.) Divide this number by 20 to reduce the weight of distance to other probes in the final adjusted score (values 0 - 5) - for small regions
	  #C.) Divide this number by 2 to increase the weigth of distance to other probes in the final adjusted score (values 0 - 50) - when overlap is allowed

	  #Note that if overlap is allowed by the user, the position will be given even more weight to encourage a better spread of probes across the region

	  #Final score range: [0 - 27 (tm)] + [0 - 10 (length)] - [0 - 50 (distance)] = (-?? to 37) where lowest score is best - when overlap is allowed
	  #Final score range: [0 - 27 (tm)] + [0 - 10 (length)] - [0 - 25 (distance)] = (-25 to 37) where lowest score is best - for large regions
 	  #Final score range: [0 - 27 (tm)] + [0 - 10 (length)] - [0 - 10 (distance)] = (-5 to 37) where lowest score is best - for small regions

	  if ($min_distance <= 0){$min_distance = 1;} #Avoid division by 0
	  my $min_distance_ratio = ($min_distance / $region_size);
	  my $min_distance_percent = ($min_distance_ratio * 100);

	  my $weighted_min_distance_percent;

	  if ($allow_overlap =~ /^yes$/i){
	    #Give most weight to probe position
	    #Compensate for the fact that the region size available is getting smaller each time when calculating the min_distance_percent
	    my $effective_region_size = ($region_size/($selected_probe_count));
	    my $effective_min_distance_ratio = ($min_distance/$effective_region_size);
	    my $effective_min_distance_percent = ($effective_min_distance_ratio * 100);
	    $weighted_min_distance_percent = ($effective_min_distance_percent / 2);
	  }elsif ($region_size > $region_size_threshold){
	    #Give more weight to probe position
	    $weighted_min_distance_percent = ($min_distance_percent / 4);
	  }else{
	    #Give less weight to probe position
	    $weighted_min_distance_percent = ($min_distance_percent / 10);
	  }

	  $probes_ref->{$probe_id}->{position_score} = -$weighted_min_distance_percent;
	  $probes_ref->{$probe_id}->{adjusted_quality_score} = $probes_ref->{$probe_id}->{quality_score} - $weighted_min_distance_percent;
	}

      }else{
	#First probe - No need to adjust score to account for position - simply pick the best one based on Tm/length
	foreach my $probe_id (keys %{$probes_ref}){
	  $probes_ref->{$probe_id}->{adjusted_quality_score} = $probes_ref->{$probe_id}->{quality_score};
	  $probes_ref->{$probe_id}->{position_score} = 0;
	}
      }

      #B.) Now go through probes for this region again (sorted by adjusted quality score)
    PROBE:foreach my $probe_id (sort {$probes_ref->{$a}->{adjusted_quality_score} <=> $probes_ref->{$b}->{adjusted_quality_score}} keys %{$probes_ref}){

	#Make sure the desired number of probes has not been reached yet - If so break out
	if ($probes_selected == $probes_needed){
	  last PASS;
	}

	#If some probes have already been selected for this region
	if ($regions_ref->{$region_id}->{selected_probes}){
	
	  my $selected_probes_ref = $regions_ref->{$region_id}->{selected_probes};

	  #If this probe has already been selected for this region, skip to the next probe
	  if ($selected_probes_ref->{$probe_id}){
	    next PROBE;
	  }

	  #Check if the current probe overlaps any of the selected probes
	  my $overlap = 0;
	  my $probe_start = $probes_ref->{$probe_id}->{unit1_start};
	  my $probe_end = $probes_ref->{$probe_id}->{unit1_end};

	  #Only test for overlap if the user does not want to allow overlaps
	  if ($allow_overlap =~ /^no$/i){
	    foreach my $selected_probe_id (keys %{$selected_probes_ref}){
	      my $test_start = $probes_ref->{$selected_probe_id}->{unit1_start};
	      my $test_end = $probes_ref->{$selected_probe_id}->{unit1_end};

	      if ($probe_start >= $test_start && $probe_start <= $test_end){$overlap = 1;}
	      if ($probe_end >= $test_start && $probe_end <= $test_end){$overlap = 1;}
	      if ($probe_start <= $test_start && $probe_end >= $test_end){$overlap = 1;}
	    }
	  }

	  #If overlap was found, proceed to the next probe - otherwise add it to the list of selected probes
	  if ($overlap == 1){
	    next PROBE;
	  }else{
	    $probes_selected++;
	    $probes_ref->{$probe_id}->{selected} = "yes";
	    $selected_probes_ref->{$probe_id}->{tmp} = '';
	    $regions_ref->{$region_id}->{selected_probe_count}++;

	    print CYAN, "\n\t\tQualityScore: $probes_ref->{$probe_id}->{quality_score}\tPostionScore: $probes_ref->{$probe_id}->{position_score}\tAdjustedQualityScore: $probes_ref->{$probe_id}->{adjusted_quality_score}", RESET;
	    next REGION;
	  }

	}else{
	  my %selected_probes;
	  $probes_selected++;
	  $probes_ref->{$probe_id}->{selected} = "yes";
	  $selected_probes{$probe_id}{tmp} = '';
	  $regions_ref->{$region_id}->{selected_probes} = \%selected_probes;
	  $regions_ref->{$region_id}->{selected_probe_count} = 1;

	  next REGION;
	}
      }
    }
  }

  print BLUE, "\n\nSelected $probes_selected probes\n", RESET;
  print LOG "\n\nSelected $probes_selected probes\n";

  return();
}


#######################################################################################################
#5.) Negative control probe selection
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
#6.) Go through the input probe files and print out a new probe files consisting only of those probes selected #
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
  foreach my $region_id (sort {$a <=> $b} keys %{$regions_ref}){

    unless ($regions_ref->{$region_id}->{probes}){
      next;
    }
    my $probes_ref = $regions_ref->{$region_id}->{probes};

    #Go through probes for this region
    foreach my $probe_id (sort keys %{$probes_ref}){
      if ($probes_ref->{$probe_id}->{selected} eq "yes"){
	$selected_region_probes{$probe_id}{printed} = "no";
      }
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

  #Open the alternate probes file and print these out as well if they were selected
  if ($extra_probes_file){
    open (REG_IN_2, "$extra_probes_file") || die "\nCould not open extra region probe input file: $extra_probes_file\n\n";

    while(<REG_IN_2>){
      chomp($_);
      my $line = $_;
      my @line = split("\t", $line);

      if ($line[0] =~ /\d+/){
	my $probe_id = $line[0];
	if ($selected_region_probes{$probe_id}){

	  #Unless it was already printed in the main file:
	  unless ($selected_region_probes{$probe_id}{printed} eq "yes"){
	    my $sequence = $line[$region_columns{Sequence}{column_pos}];
	    $selected_region_probes{$probe_id}{sequence} = $sequence;
	    $selected_region_probes{$probe_id}{printed} = "yes";

	    print REG_OUT "$line\n";
	  }
	}
      }
    }
    close (REG_IN_2);
  }

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
#7.) Print out the design file (probe ID, region ID, sequence)                                                 #
################################################################################################################
sub printDesignFile{
  my %args = @_;
  my $design_file = $args{'-design_file'};

  print BLUE, "\n\nPrinting design file: $design_file\n", RESET;
  print LOG "\n\nPrinting design file: $design_file\n";

  open (DESIGN, ">$design_file") || die "\nCould not open design file: $design_file\n\n", RESET;

  print DESIGN "Probe_ID\tRegion_ID\tSequence\n";

  #First print out the selected region probes
  foreach my $region_id (sort {$regions_ref->{$a}->{record_count} <=> $regions_ref->{$b}->{record_count}} keys %{$regions_ref}){

    unless ($regions_ref->{$region_id}->{probes}){
      next;
    }
    my $probes_ref = $regions_ref->{$region_id}->{probes};

    #Go through probes for this region
    foreach my $probe_id (sort {$probes_ref->{$a}->{unit1_start} <=> $probes_ref->{$b}->{unit1_start}} keys %{$probes_ref}){
      if ($selected_region_probes{$probe_id}){
	print DESIGN "$probe_id\t$region_id\t$selected_region_probes{$probe_id}{sequence}\n";
      }
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
#8.) Create a summary of the design for all regions.  Each record will contain the complete region record but also the following information #
#    - Number of filtered probes found, number of probes selected, status (no probes, partial probeset, full probeset)                       #
##############################################################################################################################################
sub printSummaryFile{
  my %args = @_;
  my $summary_file = $args{'-summary_file'};

  print BLUE, "\n\nPrinting summary file: $summary_file\n", RESET;
  print LOG "\n\nPrinting summary file: $summary_file\n";

  open (SUMMARY, ">$summary_file") || die "\nCould not open summary file: $summary_file\n\n", RESET;

  print SUMMARY "$input_header\tProbe_Count\tFinal_Status\n";

  my $no_probes_count = 0;
  my $partial_count = 0;
  my $full_count = 0;

  foreach my $region_id (sort {$regions_ref->{$a}->{record_count} <=> $regions_ref->{$b}->{record_count}} keys %{$regions_ref}){

    if ($regions_ref->{$region_id}->{selected_probe_count} == 0){
      $no_probes_count++;
      print SUMMARY "$regions_ref->{$region_id}->{line}\t$regions_ref->{$region_id}->{selected_probe_count}\tNo Probes\n";

    }elsif ($regions_ref->{$region_id}->{selected_probe_count} < $region_probeset_threshold){
      $partial_count++;
      print SUMMARY "$regions_ref->{$region_id}->{line}\t$regions_ref->{$region_id}->{selected_probe_count}\tPartial Probeset\n";

    }else{
      $full_count++;
      print SUMMARY "$regions_ref->{$region_id}->{line}\t$regions_ref->{$region_id}->{selected_probe_count}\tFull Probeset\n";
    }
  }

  print BLUE, "\n\tDesired probeset size for each region was: $region_probeset_threshold", RESET;
  print BLUE, "\n\tFull probesets: $full_count\tPartial probesets: $partial_count\tNo probes: $no_probes_count\n\n", RESET;
  print LOG "\n\tDesired probeset size for each region was: $region_probeset_threshold";
  print LOG "\n\tFull probesets: $full_count\tPartial probesets: $partial_count\tNo probes: $no_probes_count\n\n";

  close(SUMMARY);

  return();
}


