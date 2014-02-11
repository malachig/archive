#!/usr/bin/perl -w
#Written by Malachi Griffith

#The purpose of this script is to map probes to target exon regions
#STEPS

#1.) Import target exon regions from an input file (region id, chromosome, start, end)

#2.) Import probe records from an input file (probe_id, region_id, start, end)

#3.) Go through each exon region, find the probes that overlap the region.  Note that some probes will not overlap any region
#      - For each probe that does overlap a region, calculate the distance between the centre of the region and the centre of the probe
#      - Note that it is possible that a probe will overlap more than one region (if the regions themselves are close or overlapping)
#      - Make sure each probe is associated with the region it is the closest to the centre of
#      - For probes that do not overlap a region, determine the region they are closest to and the distance from the centre of this region

#4.) Update the region probe file with the exon region ID it maps to (closest one), the distance from the centre of this region and whether is actually overlaps
#    - Append: TargetExonRegion_ID, DistanceFromRegion, Overlaps

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);

#Initialize command line options
my $target_regions = '';
my $region_probes = '';
my $out_file = '';
my $log_file = '';

GetOptions ('target_regions=s'=>\$target_regions, 'region_probes=s'=>\$region_probes, 'out_file=s'=>\$out_file, 'log_file=s'=>\$log_file);

print GREEN, "\n\nThis script generates an array design submission file for all regions in a specified in an input file\n\n", RESET;

print GREEN, "\nParameters:", RESET;
print GREEN, "\n\tSpecify the path to the target region file using: --target_regions", RESET;
print GREEN, "\n\tSpecify the path to the region probes file using: --region_probes", RESET;
print GREEN, "\n\tSpecify the path to the new appended region probes file using: --out_file", RESET;

print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tmapProbesToTargetExonRegions.pl  --target_regions=/projects/malachig/GTAC_Chip_Design/TargetExonRegions.txt  --region_probes=/projects/malachig/GTAC_Chip_Design/probes/filtered/regionProbes_filtered.txt  --out_file=/projects/malachig/GTAC_Chip_Design/probes/filtered/regionProbes_filtered_mapped.txt  --log_file=/projects/malachig/GTAC_Chip_Design/logs/mapProbesToTargetExonRegions_LOG.txt\n\n", RESET;

unless ($target_regions && $region_probes && $out_file && $log_file){
  print RED, "\nRequired input parameter(s) missing or incorrect format!\n\n", RESET;
  exit();
}

open (LOG, ">$log_file") || die "\nCould not open log file: $log_file!\n\n";


#1.) Import target regions from an input file (region id, chromosome, start, end)
my $input_header;
my $chromosomes_ref = &importTargetRegions('-region_file'=>$target_regions);

#2.) Import region probe records from an input file (probe_id, region_id, start, end)
my %region_columns;
my $probes_ref = &importRegionProbes('-input_file'=>$region_probes);


#3.) Go through each exon region, find the probes that overlap the region.  Note that some probes will not overlap any region
#      - For each probe that does overlap a region, calculate the distance between the centre of the region and the centre of the probe
#      - Note that it is possible that a probe will overlap more than one region (if the regions themselves are close or overlapping)
#      - Make sure each probe is associated with the region it is the closest to the centre of
#      - For probes that do not overlap a region, determine the region they are closest to and the distance from the centre of this region
&mapProbesToExonRegions();

#4.) Update the probe file
&updateProbeFile('-probe_file'=>$region_probes, '-out_file'=>$out_file);




print BLUE, "\n\n", RESET;

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
      $regions_ref->{$region_id}->{mapped_probe_list} = \%mapped_probe_list; #List of probe that will be assigned to this region.  Probes that are closest to or overlap this region

      my %selected_probe_list;
      $regions_ref->{$region_id}->{selected_probe_list} = \%selected_probe_list; #List of probe that will from this region that are actually selected for the array

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
      $regions{$region_id}{mapped_probe_list} = \%mapped_probe_list; #List of probe that will be assigned to this region.  Probes that are closest to or overlap this region

      my %selected_probe_list;
      $regions{$region_id}{selected_probe_list} = \%selected_probe_list; #List of probe that will from this region that are actually selected for the array

      $chromosomes{$chromosome}{regions} = \%regions;
    }
  }
  close (REGIONS);



  print BLUE, "\n\n\tFound $region_count exon regions\n\n", RESET;
  print LOG "\n\n\tFound $region_count exon regions\n\n";

  return(\%chromosomes);
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
      unless ($region_columns{Probe_Count}){
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
    #Probe_Count, ProbeSet_ID, Region_ID, Probe_length, Probe_Tm, Unit1_start, Unit1_end

    my $probeset_id = $line[$region_columns{ProbeSet_ID}{column_pos}];
    my $chromosome = $line[$region_columns{Chromosome}{column_pos}];
    my $unit1_start = $line[$region_columns{Unit1_start}{column_pos}];
    my $unit1_end = $line[$region_columns{Unit1_end}{column_pos}];

    if ($unit1_end < $unit1_start){
      print RED, "\nEnd coordinate of probe $probe_id is smaller than start coordinate!\n\n", RESET;
      exit();
    }

    $probes{$probe_id}{probeset_id} = $probeset_id;
    $probes{$probe_id}{chromosome} = $chromosome;
    $probes{$probe_id}{unit1_start} = $unit1_start;
    $probes{$probe_id}{unit1_end} = $unit1_end;

  }
  return(\%probes);
}


#######################################################################################################
#4.) Determine the distance for each probe from its closest exon region                               #
#######################################################################################################
sub mapProbesToExonRegions{

  print BLUE, "\n\nMapping all probes to their closest region in the exon region list.  Noting actual overlaps as well\n", RESET;
  print LOG "\n\nMapping all probes to their closest region in the exon region list.  Noting actual overlaps as well\n";

  #- Go through each exon region, find the probes that overlap the region.  Note that some probes will not overlap any region
  #- For each probe that does overlap a region, calculate the distance between the centre of the region and the centre of the probe
  #- Note that it is possible that a probe will overlap more than one region (if the regions themselves are close or overlapping)
  #- Make sure each probe is associated with the region it is the closest to the centre of
  # - For probes that do not overlap a region, determine the region they are closest to and the distance from the centre of this region

  my $overlaping_probes_found = 0;
  my $probes_mapped = 0;
  my $counter = 0;

  #Go through each probe.  Find the distance from the centre of this probe to all regions on the same chromosome
  foreach my $probe_id (sort {$a <=> $b} keys %{$probes_ref}){
    $counter++;

    if ($counter == 10000){
      $| = 1;
      print ".";
      $| = 0;
      $counter = 0;
    }

    my $probe_chr = $probes_ref->{$probe_id}->{chromosome};
    my $probe_start = $probes_ref->{$probe_id}->{unit1_start};
    my $probe_end = $probes_ref->{$probe_id}->{unit1_end};

    #Sanity check
    if ($probe_start > $probe_end){
      print RED, "\nProbe start is larger than end!\n\n", RESET;
      exit();
    }

    unless ($chromosomes_ref->{$probe_chr}){
      print RED, "\nCould not find any target regions for the chromosome of probe: $probe_id ($probe_chr)\n\n", RESET;
      exit();
    }

    my $regions_ref = $chromosomes_ref->{$probe_chr}->{regions};

    my $temp = $probe_start + (($probe_end - $probe_start)/2);
    my $probe_centre = sprintf("%.0f", $temp);

    my $smallest_diff = 1000000000000000000000000;  #arbitrarily large value
    my $best_region_id;  #ID of region that the probe is closest to
    my $overlap_found = 0;  #Whether the probe actually overlaps the best region

    foreach my $region_id (sort {$regions_ref->{$a}->{region_centre} <=> $regions_ref->{$b}->{region_centre}} keys %{$regions_ref}){

      my $region_centre = $regions_ref->{$region_id}->{region_centre};

      my $centre_diff = abs($probe_centre - $region_centre);
      if ($centre_diff <= $smallest_diff){
	$smallest_diff = $centre_diff;
	$best_region_id = $region_id;
      }else{
	last();
      }
    }

    if ($best_region_id){
      $probes_mapped++;
    }else{
      print RED, "\nCould not find a closest region for probe: $probe_id\t($probe_chr:$probe_start-$probe_end)\n\n", RESET;
      exit();
    }
    #Store the best region found for this probe (one the probe is closest to)
    $probes_ref->{$probe_id}->{closest_region_id} = $best_region_id;
    $probes_ref->{$probe_id}->{distance_to_closest_region} = $smallest_diff;

    #Check if this probe actually overlap the closest region
    my $region_start = $regions_ref->{$best_region_id}->{start};
    my $region_end = $regions_ref->{$best_region_id}->{end};

    if ($probe_start >= $region_start && $probe_start <= $region_end){$overlap_found = 1;}
    if ($probe_end >= $region_start && $probe_end <= $region_end){$overlap_found = 1;}
    if ($probe_start <= $region_start && $probe_end >= $region_end){$overlap_found = 1;}

    if ($overlap_found == 1){
      $probes_ref->{$probe_id}->{overlaps_closest_region} = "1";
      $overlaping_probes_found++;
    }else{
      $probes_ref->{$probe_id}->{overlaps_closest_region} = "0";
    }

    #print YELLOW, "\nProbe: $probe_id is closest to region $best_region_id (distance = $smallest_diff bp) (overlap = $overlap_found)", RESET;
  }


  print BLUE, "\tMapped $probes_mapped probes and found $overlaping_probes_found that actually overlap an exon region\n", RESET;
  print LOG "\tMapped $probes_mapped probes and found $overlaping_probes_found that actually overlap an exon region\n";

  return();
}


##########################################################################################################################################
#4.) Update the probe file
##########################################################################################################################################
sub updateProbeFile{
  my %args = @_;
  my $probe_file = $args{'-probe_file'};
  my $out_file = $args{'-out_file'};

  #Open the probe file and read the neccessary probe data into a hash keyed on probe ID
  open (PROBES, "$probe_file") || die "\nCould not open input probe file: $probe_file\n\n";
  open (OUT, ">$out_file") || die "\nCould not open output probe file: $out_file\n\n";

  my $first_line = 1;

  print BLUE, "\nUpdating Region probe records from: $probe_file to $out_file\n", RESET;
  print LOG "\nUpdating Region probe records from: $probe_file to $out_file\n";

  my $blocks_imported = 0;
  my $progress_count = 0;

  while (<PROBES>){
    $progress_count++;
    if ($progress_count == 10000){
      $blocks_imported++;
      print BLUE, "\n\tImported $blocks_imported blocks of 10,000 probes", RESET;
      print LOG "\n\tImported $blocks_imported blocks of 10,000 probes";
      $progress_count = 0;
    }
    chomp($_);
    my $line = $_;
    my @line = split("\t", $line);

    #Get the header line and identify column names and their positions
    if ($first_line == 1){

      print OUT "$line\tTargetExonRegion_ID\tDistanceFromExonRegion\tOverlapsExonRegion\n";

      my $col_count = 0;
      foreach my $column (@line){
	$region_columns{$column}{column_pos} = $col_count;
	$col_count++;
      }
      $first_line = 0;

      #Check for critical columns and their names
      unless ($region_columns{Probe_Count}){
	print RED, "\nCritical column missing or named incorrectly, check input file", RESET;
	exit();
      }
      next();
    }

    #Get the values of interest from each line (probe record)
    my $probe_id = $line[$region_columns{Probe_Count}{column_pos}];

    print OUT "$line\t$probes_ref->{$probe_id}->{closest_region_id}\t$probes_ref->{$probe_id}->{distance_to_closest_region}\t$probes_ref->{$probe_id}->{overlaps_closest_region}\n";

  }

  close(PROBES);
  close(OUT);


  return();
}
