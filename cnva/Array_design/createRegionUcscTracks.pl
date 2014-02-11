#!/usr/bin/perl -w

#Written by Malachi Griffith
#The purpose of this script is to create custom UCSC track to show the position of target regions
#One file will be created per chromosome

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;

#Initialize command line options
my $region_file = '';
my $region_probe_file = '';
my $ucsc_dir = '';
my $ucsc_build = '';
my $web_path = '';
my $gzip_path = '';

GetOptions ('region_file=s'=>\$region_file, 'region_probe_file=s'=>\$region_probe_file, 'ucsc_dir=s'=>\$ucsc_dir, 'ucsc_build=s'=>\$ucsc_build,
	    'web_path=s'=>\$web_path, 'gzip_path=s'=>\$gzip_path);

print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\nExample: createRegionUcscTracks.pl  --region_file=/home/malachig/MR_Chip_Design/MasterRegionList-hg18-Malachi_V3.txt  --region_probe_file=/home/malachig/MR_Chip_Design/probes/selected_probes/MR_CNV_V1_72k_hg18_RegionProbes.txt(OPTIONAL)  --ucsc_dir=/home/malachig/www/public/htdocs/MR_Chip_Design/TargetRegions_hg18_V3/  --web_path='http://www.bcgsc.ca/people/malachig/htdocs/MR_Chip_Design/TargetRegions_hg18_V3/'  --ucsc_build=hg18  --gzip_path=/usr/bin/gzip\n\n", RESET;

#Make sure all options were specified
unless ($region_file && $ucsc_dir && $ucsc_build && $web_path && $gzip_path){
  print RED, "\nRequired input parameter(s) missing!\n\n", RESET;
  exit();
}

#Import the regions of interest - perform basic sanity checks on the input file and organize regions accorinding their target sequence (chromosome)
my %regions;
my %sequence_names;
my $input_header;
print BLUE, "\nImporting target regions from input file\n\n", RESET;
&importTargetRegions('-region_file'=>$region_file);


my %region_columns;
if ($region_probe_file){
  print BLUE, "\nImporting target region probes from input file\n\n", RESET;
  &importRegionProbes('-probe_file'=>$region_probe_file, '-regions_object'=>\%regions);
}

unless ($ucsc_dir =~ /.*\/$/){
  $ucsc_dir = "$ucsc_dir"."/";
}
unless ($web_path =~ /.*\/$/){
  $web_path = "$web_path"."/";
}

#Make sure the ucsc_dir provided exists and is a directory
unless (-e $ucsc_dir && -d $ucsc_dir){
  print RED, "\nSpecified directory: $ucsc_dir does not appear valid!\n\n", RESET;
  exit();
}

#Go through each chromosome and print a track file for the regions on it
foreach my $chr (sort keys %sequence_names){

  my @regions = @{$sequence_names{$chr}{regions}};

  my $track_file = "$ucsc_dir"."$chr".".txt";
  open (UCSC, ">$track_file") || die "\nCould not open trackfile: $track_file\n\n";

  #Browser line
  print UCSC "#Browser line";
  print UCSC "\nbrowser hide all";
  print UCSC "\nbrowser full knownGene";
  print UCSC "\nbrowser pack multiz28way";


  #Track line for regions track
  print UCSC "\n\n#Track line";
  print UCSC "\ntrack name=target_regions description=\"Target Regions - MR_CNV_V1_385k_hg18 - Malachi\" color=0,0,255 useScore=0 visibility=3";
  print UCSC "\n\n#Begin DATA";

  foreach my $region_id (@regions){
    my $start = $regions{$region_id}{start};
    my $end = $regions{$region_id}{end};
    my $region_name = "Region_$region_id";
    print UCSC "\n$chr\tMG\tTargetRegion\t$start\t$end\t.\t.\t.\t$region_name";
  }

  if ($region_probe_file){
    #Track line for probes track
    print UCSC "\n\n#Track line";
    print UCSC "\ntrack name=selected_probes description=\"Microarray Probes - MR_CNV_V1_385k_hg18 - Malachi\" color=0,100,0 useScore=0 visibility=3";
    print UCSC "\n\n#Begin DATA";

    foreach my $region_id (@regions){
      my $probes_ref = $regions{$region_id}{probes};

      foreach my $probe_id (sort {$probes_ref->{$a}->{unit1_start} <=> $probes_ref->{$b}->{unit1_start}} keys %{$probes_ref}){
	my $start = $probes_ref->{$probe_id}->{unit1_start};
	my $end = $probes_ref->{$probe_id}->{unit1_end};
	my $region_name = "Probe_$probe_id";
	print UCSC "\n$chr\tMG\tProbe\t$start\t$end\t.\t+\t.\t$region_name";
      }
    }
  }

  close (UCSC);
}


#Go to the target dir and compress the track files
my $command = "$gzip_path". " $ucsc_dir"."*";
print BLUE, "\nCompressing custom track files with command: $command\n\n", RESET;
system ("$command");

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

    $regions{$region_id}{record_count} = $record_count;
    $regions{$region_id}{chromosome} = $chromosome;
    $regions{$region_id}{start} = $start;
    $regions{$region_id}{end} = $end;
    $regions{$region_id}{region_size} = ($end - $start)+1;
    $regions{$region_id}{line} = $line;

    if ($sequence_names{$chromosome}){
      $sequence_names{$chromosome}{count}++;
      push (@{$sequence_names{$chromosome}{regions}}, $region_id);
    }else{
      my @regions;
      push (@regions, $region_id);
      $sequence_names{$chromosome}{count} = 1;
      $sequence_names{$chromosome}{regions} = \@regions;
    }

  }
  close (REGIONS);

  my $region_count = keys %regions;
  my $seq_count = keys %sequence_names;

  print BLUE, "\n\n\tFound $region_count regions corresponding to $seq_count sequences (chromosomes)\n\n", RESET;

  return();
}

######################################################################################################################################
#Import the probes selected for these regions                                                                                        #
######################################################################################################################################
sub importRegionProbes{
  my %args = @_;
  my $probe_file = $args{'-probe_file'};
  my $regions_ref = $args{'-regions_object'};

  my $progress_count = 0;
  my $blocks_imported = 0;
  my $probe_count = 0;

  #Open the probe file and read the neccessary probe data into a hash keyed on probe ID
  open (PROBES, "$probe_file") || die "\nCould not open input probe file: $probe_file\n\n";

  my $first_line = 1;

  print BLUE, "\nImporting Region probe records from: $probe_file\n", RESET;

  while (<PROBES>){
    $progress_count++;
    if ($progress_count == 10000){
      $blocks_imported++;
      print BLUE, "\n\tImported $blocks_imported blocks of 10,000 probes", RESET;
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
    my $unit1_start = $line[$region_columns{Unit1_start}{column_pos}];
    my $unit1_end = $line[$region_columns{Unit1_end}{column_pos}];

    if ($unit1_end < $unit1_start){
      print RED, "\nEnd coordinate of probe $probe_id is smaller than start coordinate!\n\n", RESET;
      exit();
    }

    if ($regions_ref->{$region_id}->{probes}){
      my $probes_ref = $regions_ref->{$region_id}->{probes};
      $probes_ref->{$probe_id}->{probeset_id} = $probeset_id;
      $probes_ref->{$probe_id}->{unit1_start} = $unit1_start;
      $probes_ref->{$probe_id}->{unit1_end} = $unit1_end;
    }else{
      my %probes;
      $probes{$probe_id}{probeset_id} = $probeset_id;
      $probes{$probe_id}{unit1_start} = $unit1_start;
      $probes{$probe_id}{unit1_end} = $unit1_end;
      $regions_ref->{$region_id}->{probes} = \%probes;
    }
  }

  return();
}
