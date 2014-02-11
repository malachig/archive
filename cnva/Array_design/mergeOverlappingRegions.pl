#!/usr/bin/perl -w

#Written by Malachi Griffith
#The purpose of this script is to merge overlapping chromosomal regions supplied in an input file
#A new ID will be assigned to the merged region
#The information associated with the first record encountered will be kept for the resulting cluster of merged regions
#The regions which were merged to create the new region will also be noted
#The user will be given the option of merging only those regions below a certain size 
# - This is done in case you want to deal with larger regions in a different way

#A new output file will be created - consisting of the merged regions

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;

BEGIN {
  my $script_dir = &File::Basename::dirname($0);
  push (@INC, $script_dir);
}
use utilities::Descriptive;

#Initialize command line options
my $region_file = '';
my $size_limit = '';
my $outfile = '';

GetOptions ('region_file=s'=>\$region_file, 'size_limit=i'=>\$size_limit, 'outfile=s'=>\$outfile);

#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify a tab-delimited file containing chromosome coordinates of target regions using: --region_file", RESET;
print GREEN, "\n\t\tFormat of file. Should have a header line and each line should have the following columns first: 'unique_id, chromosome, start, end'", RESET;
print GREEN, "\n\tSpecify the max size of regions to be considered for merging using:  --size_limit", RESET;
print GREEN, "\n\tSpecify the name of the output merged region file: --outfile", RESET;
print GREEN, "\n\nExample: mergeOverlappingRegions.pl  --region_file=MasterRegionList-hg18-Malachi_FinalBeforeMerging.txt  --size_limit=602  --outfile=MasterRegionList-hg18-Malachi_FinalAfterMerging_V1.txt\n\n", RESET;

unless ($region_file && $size_limit && $outfile){
  print RED, "\nRequired input parameter(s) missing!\n\n", RESET;
  exit();
}

#Import the regions of interest - perform basic sanity checks on the input file and organize regions accorinding their target sequence (chromosome)
my %regions;
my %sequence_names;
my $input_header;
my $current_max_region_id = 0;

&importTargetRegions('-region_file'=>$region_file);

open (OUT, ">$outfile") || die "\nCould not open output file: $outfile\n\n";
print OUT "$input_header\tMergedRegionId\tMerged_Regions\tMerged_Start\tMerged_End\n";

my %merged_regions;
&mergeTargetRegions('-region_object'=>\%regions);







close (OUT);


exit();


######################################################################################################################################
#Import the regions of interest                                                                                                      #
#- perform basic sanity checks on the input file and organize regions accorinding their target sequence (chromosome)                 #
######################################################################################################################################
sub importTargetRegions{
  my %args = @_;
  my $region_file = $args{'-region_file'};

  print YELLOW, "\n\nImporting target regions from: $region_file", RESET;

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
    my $chromosome = $line[2];
    my $start = $line[3];
    my $end = $line[4];

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

    if ($region_id > $current_max_region_id){
      $current_max_region_id = $region_id;
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

  return();
}


######################################################################################################################################
#Merge target regions                                                                                                                #
######################################################################################################################################
sub mergeTargetRegions{
  my %args = @_;
  my $regions_ref = $args{'-region_object'};
  my $regions_processed = 0;
  my $merged_regions = 0;

  print BLUE, "\nChecking for overlaping regions and merging those (unless they are larger than $size_limit)\n\n", RESET;

  foreach my $chr (sort keys %sequence_names){

    print BLUE, "\n\tProcessing regions on $chr", RESET;

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

      #Watch for regions larger than the max size to be considered - mark these special cases and skip to the next region
      my $size = ($current_end - $current_start) + 1;
      if ($size > $size_limit){
	$cluster_count++;
	my %ol_regions;
	$ol_regions{$region_id}{start} = $regions_ref->{$region_id}->{start};
	$ol_regions{$region_id}{end} = $regions_ref->{$region_id}->{end};
	
	$clusters{$cluster_count}{ol_regions} = \%ol_regions;
	$clusters{$cluster_count}{chromosome} = $chr;
	$clusters{$cluster_count}{large_region} = "yes";
	next();
      }

      #Go through each cluster defined so far and check for overlaps to the current region
    CLUSTER:foreach my $cluster_id (sort {$a <=> $b} keys %clusters){

	#Do not consider overlap to regions above the allowed size limit
	if ($clusters{$cluster_id}{large_region} eq "yes"){
	  next();
	}

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
	$clusters{$cluster_count}{large_region} = "no";
	
      }

    }

    #DEBUG: Print dumper of clusters for this chromosome
    #print Dumper %clusters;

    #Go through each of the clusters defined for this chromosome and print new merged regions
    foreach my $cluster_id (sort {$a <=> $b} keys %clusters){

      #If a cluster consists of only a single region, simply keep the original record
      my $ol_regions_ref = $clusters{$cluster_id}{ol_regions};
      my $cluster_size = keys %{$ol_regions_ref};
      $regions_processed += $cluster_size;

      if ($cluster_size == 1){
	#Process cluster of only a single region
	foreach my $region_id (sort keys %{$ol_regions_ref}){
	  my $line = $regions_ref->{$region_id}->{line};
	  my $merged_start = $regions_ref->{$region_id}->{start};
	  my $merged_end = $regions_ref->{$region_id}->{end};

	  print OUT "$line\t$region_id\t$region_id\t$merged_start\t$merged_end\n";
	  $merged_regions++;
	}

      }else{
	#Process actual cluster of overlapping regions
	#First get the min region ID of the cluster as well as the new start/end coordinates of the merged region
	$current_max_region_id++;

	my @ids;
	my @starts;
	my @ends;

	foreach my $region_id (sort keys %{$ol_regions_ref}){
	  push (@ids, $region_id);
	  push (@starts, $regions_ref->{$region_id}->{start});
	  push (@ends, $regions_ref->{$region_id}->{end});
	}
	my $stat;

	$stat = Statistics::Descriptive::Full->new();
	$stat->add_data(@ids);
	my $min_region_id = $stat->min();

	$stat = Statistics::Descriptive::Full->new();
	$stat->add_data(@starts);
	my $merged_start = $stat->min();

	$stat = Statistics::Descriptive::Full->new();
	$stat->add_data(@ends);
	my $merged_end = $stat->max();

	my $line = $regions_ref->{$min_region_id}->{line};
	my @merged_regions = keys %{$ol_regions_ref};

	my @sorted_merged_regions = sort {$a <=> $b} (@merged_regions);

	print OUT "$line\t$current_max_region_id\t@sorted_merged_regions\t$merged_start\t$merged_end\n";
	$merged_regions++;
      }
    }
  }

  print BLUE, "\n\nProcessed a total of $regions_processed regions to identify overlaps", RESET;
  print BLUE, "\nPrinted a total of $merged_regions merged regions to the output file", RESET;
  return();
}
