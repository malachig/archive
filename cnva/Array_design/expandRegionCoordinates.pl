#!/usr/bin/perl -w
#Written by Malachi Griffith

#This script will expand the size of genomic regions provided in an input file.
#It assumes the regions have already been tested to ensure they do not overlap!
#It will attempt to increase the size of each region in the input file by adding the specified amount of flank
#It will check the space available on each side of the region and add the flank only if enough space is available
#When determining this it also assumes that the same flank will be added to the next region
#If there is not enough space, it will add as much as possible without causing an overlap with adjacent regions 
#    - again assuming that the same flank is added to these

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);

#Initialize command line options
my $target_regions = '';
my $sequence_dir = '';
my $flank = '';
my $out_file = '';

GetOptions ('target_regions=s'=>\$target_regions, 'sequence_dir=s'=>\$sequence_dir, 'flank=i'=>\$flank, 'out_file=s'=>\$out_file);

print GREEN, "\n\nThis script adds flank to a list of NON-OVERLAPING chromosomal regions\n\n", RESET;
print GREEN, "\nParameters:", RESET;
print GREEN, "\n\tSpecify the path to the target region file using: --target_regions", RESET;
print GREEN, "\n\tSpecify the path to a directory of chromosome fasta files using: --sequence_dir", RESET;
print GREEN, "\n\tSpecify the desired amount of flank to add to each side of every region using: --flank", RESET;
print GREEN, "\n\tSpecify the path to the updated target region file using: --outfile", RESET;
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\texpandRegionCoordinates.pl  --target_regions=MasterRegionList-hg18-Malachi_FinalAfterMerging_V8.txt  --sequence_dir=/home/malachig/MR_Chip_Design/hg18_genome_fasta  --flank=300  --out_file=MasterRegionList_plusFlank300.txt\n\n", RESET;

unless ($target_regions && $sequence_dir && $flank && $out_file){
  print RED, "\nRequired input parameter(s) missing or incorrect format!\n\n", RESET;
  exit();
}

my %sequence_names;
my $input_header;
my $regions_ref = &importTargetRegions('-region_file'=>$target_regions);

#Test for region overlaps
&identifyRegionOverlaps('-region_object'=>$regions_ref);

#Within each sequence (chromosome) assign an index to each region based on their order
foreach my $chr (sort keys %sequence_names){
  #Get the regions for this chromosome
  my $region_list_ref = $sequence_names{$chr}{region_list};
  my $index = 0;

  my @regions_index;

  #Sort the regions by start position and go through them one by one comparing coordinates to all previously considered regions
  foreach my $region_id (sort {$region_list_ref->{$a}->{start} <=> $region_list_ref->{$b}->{start}} keys %{$region_list_ref}){
    $region_list_ref->{$region_id}->{position_index} = $index;
    push (@regions_index, $region_id);
    $index++;
  }
  $sequence_names{$chr}{regions_index} = \@regions_index;
}



#Get the length of each chromosome involved in the regions imported
chomp($sequence_dir);
unless ($sequence_dir =~ /\/$/){
  $sequence_dir = "$sequence_dir"."/";
}
unless (-e $sequence_dir && -d $sequence_dir){
  print RED, "\nSpecified sequence directory: $sequence_dir does not appear valid!\n\n", RESET;
  exit();
}
&getSequenceLengths('-sequence_dir'=>$sequence_dir, '-sequence_names'=>\%sequence_names);

print BLUE, "\n\nProcessing each chromosome and expanding the regions on each where possible\n", RESET;

#Working one chromosome at a time, go through each region and attempt to add flank to it
foreach my $chr (sort keys %sequence_names){
  print YELLOW, "\n\tExpanding regions on $chr", RESET;
  my $regions_expanded = 0;

  #Get the regions for this chromosome
  my $region_list_ref = $sequence_names{$chr}{region_list};

  my $seq_region_count = keys %{$region_list_ref};
  my $region_count = 0;

  my @regions_index = @{$sequence_names{$chr}{regions_index}};
  my $chr_start = $sequence_names{$chr}{chr_start};
  my $chr_end = $sequence_names{$chr}{chr_end};

  #Sort the regions by start position and go through them one by one comparing coordinates to all previously considered regions
  foreach my $region_id (sort {$region_list_ref->{$a}->{start} <=> $region_list_ref->{$b}->{start}} keys %{$region_list_ref}){
    $region_count++;

    my $start = $regions_ref->{$region_id}->{start};
    my $end = $regions_ref->{$region_id}->{end};

    my $new_start;
    my $new_end;

    #Sanity check - make sure end coord > start coord
    unless ($end > $start){
      print RED, "\nEnd coordinate is smaller than start coordinate! - not allowed\n\n", RESET;
      exit();
    }

    my $current_index = $region_list_ref->{$region_id}->{position_index};

    #A.) Deal with the left side of this region
    if (($current_index - 1) >= 0){
      my $left_region_id = $regions_index[$current_index - 1];

      #Calculate the distance to the region on the left
      my $dist_to_left = $start - $regions_ref->{$left_region_id}->{end};

      if($dist_to_left < 5){
	#Regions are already beside each other
	$new_start = $start;

      }elsif ($dist_to_left < $flank){
	#Regions have some space between them but not enough for full flank distance
	$new_start = ($start - ($dist_to_left - 2));

      }else{
	#Regions have enough space between them to add the complete flank distance
	$new_start = ($start - $flank);
      }

    }else{
      #No region to the left at all
      $new_start = ($start - $flank);
      print CYAN, "\n\tNo region to left - adding flank", RESET;
    }


    #B.) Deal with the right side of this region
    if (exists $regions_index[$current_index + 1]){
      my $right_region_id = $regions_index[$current_index + 1];

      #Calculate the distance to the region on the left
      my $dist_to_right = $regions_ref->{$right_region_id}->{start} - $end;
      my $test_to_right = $regions_ref->{$right_region_id}->{start} - (($flank * 2) + 2);

      if($dist_to_right < 5){
	#Regions are already beside each other
	$new_end = $end;

      }elsif ($end > $test_to_right){
	#Region to the right is within (2 x flank) of the current region
	#Keep in mind that the region to the right will also be expanded (next iteration)
	my $x = $dist_to_right/2;
	my $allowable_flank;
	if ($x =~ /^(\d+)/){
	  $allowable_flank = $1;
	}else{
	  print RED, "\nUnrecognized number format!\n\n", RESET;
	  exit();
	}
	$new_end = $end + ($allowable_flank - 2);

      }else{
	#Regions have enough space between them to add the complete flank distance
	$new_end = ($end + $flank);
      }

    }else{
      #No region to the right at all
      $new_end = ($end + $flank);
      print CYAN, "\n\tNo region to right - adding flank", RESET;
    }

    #Sanity test, new end must be larger than new start
    unless ($new_end > $new_start){
      print RED, "\nEnd coordinate is smaller than start coordinate! - not allowed\n\n", RESET;
      exit();
    }

    #Final test - Make sure that expansion of this region did not exceed the boundaries of the chromosome!
    #If this is okay, change the coordinates of this region to the new coords with flank
    if (($new_start-1 < $chr_start) || ($new_end+1 > $chr_end)){
      print "\n\tNew coordinates out of chromosome bounds - leaving untouched\n", RESET;
    }else{

      unless (($new_start == $start) && ($new_end == $end)){
	$regions_expanded++;
      }

      $regions_ref->{$region_id}->{start} = $new_start;
      $regions_ref->{$region_id}->{end} = $new_end;

    }
  }
  print YELLOW, "\n\t\tFound $region_count regions and expanded $regions_expanded", RESET;
}

#Print out input file with new coordinates
print BLUE, "\n\nPrinting new output file: $out_file with new coordinates appended\n\n", RESET;
open (OUT, ">$out_file") || die "\nCould not open output file: $out_file\n\n";
print OUT "$input_header\tNewStart\tNewEnd\n";
foreach my $region_id (sort {$regions_ref->{$a}->{record_count} <=> $regions_ref->{$b}->{record_count}} keys %{$regions_ref}){
  print OUT "$regions_ref->{$region_id}->{line}\t$regions_ref->{$region_id}->{start}\t$regions_ref->{$region_id}->{end}\n";
}
close(OUT);

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

  print BLUE, "\n\nImporting target regions from: $region_file", RESET;

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

  print YELLOW, "\n\n\tFound $region_count regions\n\n", RESET;

  return(\%regions);
}


######################################################################################################################
#Test for region overlaps                                                                                            #
######################################################################################################################
sub identifyRegionOverlaps{
  my %args = @_;
  my $regions_ref = $args{'-region_object'};
  my $regions_processed = 0;
  my $overlaping_regions = 0;

  print BLUE, "\nChecking all regions against each other for overlaps\n", RESET;

  foreach my $chr (sort keys %sequence_names){

    print YELLOW, "\n\tProcessing regions on $chr", RESET;

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

	foreach my $region_id (sort {$ol_regions_ref->{$a}->{start} <=> $ol_regions_ref->{$b}->{start}} keys %{$ol_regions_ref}){
	  my $size = ($ol_regions_ref->{$region_id}->{end} - $ol_regions_ref->{$region_id}->{start}) + 1;
	  print YELLOW, "\n\t\tRegion $region_id ($chr:$ol_regions_ref->{$region_id}->{start}-$ol_regions_ref->{$region_id}->{end} Size = $size)", RESET;
	}
      }
    }
  }

  print BLUE, "\n\nProcessed $regions_processed regions and found a total of $overlaping_regions clusters of overlapping regions\n", RESET;

  if ($overlaping_regions > 0){
    print RED, "\nFound overlaping regions - resolve these before proceeding\n\n", RESET;
    exit();
  }
  return();
}


######################################################################################################################
#Check sequence files in the specified directory                                                                     #
######################################################################################################################
sub getSequenceLengths{
  my %args = @_;
  my $sequence_dir = $args{'-sequence_dir'};
  my $sequence_names_ref = $args{'-sequence_names'};

  print BLUE, "\nGetting the length of each chromosome\n", RESET;

  #Go through each sequence name and make sure a corresponding masked and unmasked sequence file exists

  foreach my $sequence_name (sort keys %{$sequence_names_ref}){
    print YELLOW, "\n\tChecking fasta file for $sequence_name", RESET;

    my $unmasked_file = "$sequence_dir"."$sequence_name".".fa";

    unless (-e $unmasked_file){
      print RED, "\nCould not find the unmasked file: $unmasked_file for the sequence: $sequence_name!\n\n", RESET;
      exit();
    }

    $sequence_names_ref->{$sequence_name}->{unmasked_file} = $unmasked_file;

    #Open the two files and get the sequence data
    my $unmasked_seq_ref = &parseFastaFile('-file'=>$unmasked_file, '-sequence_name'=>$sequence_name);

    my $unmasked_seqs_found = keys %{$unmasked_seq_ref};

    if ($unmasked_seqs_found > 1){
      print RED, "\nFound more than one sequence in the file: $unmasked_file - Only one chromosome expected per file!\n\n", RESET;
      exit();
    }

    #Compare length of masked and unmasked sequences
    my $unmasked_length = length(${$unmasked_seq_ref->{$sequence_name}->{sequence}});

    $sequence_names_ref->{$sequence_name}->{chr_start} = 1;
    $sequence_names_ref->{$sequence_name}->{chr_end} = $unmasked_length - 1;

    print YELLOW, "\tCoords = $sequence_names_ref->{$sequence_name}->{chr_start} - $sequence_names_ref->{$sequence_name}->{chr_end}", RESET;

    $unmasked_seq_ref->{$sequence_name}->{sequence} = ();
  }
  return();
}


################################################################################################################
#Parse a fasta file with a single sequence record within in it (and the name is specified)
################################################################################################################
sub parseFastaFile{
  my %args = @_;
  my $file = $args{'-file'};
  my $sequence_name = $args{'-sequence_name'};

  my %seqs;

  my $input_seperator = $/;
  $/ = ">";

  open (FASTA, "$file") || die "\nCould not open file: $file\n\n";
  while(<FASTA>){

    my $name;
    my $seq;
    if ($_ =~ /^(chr.*?)\n(.*)/ms){
      $name = $1;
      $seq = $2;
      $seq =~ s/\n//g;

    }elsif($_ =~ /^\>$/){
      next();
    }else{
      print RED, "\nRecord not understood: $file", RESET;
      exit();
    }

    unless ($name eq $sequence_name){
      print RED, "\nThe sequence name within the sequence file ($file) does not match the expected sequence name!\n\n", RESET;
      exit();
    }

    $seq = uc($seq);
    $seqs{$sequence_name}{sequence} = \$seq;
  }
  close (FASTA);

  $/ = $input_seperator;

  return(\%seqs);
}

