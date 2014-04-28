#!/usr/bin/perl -w
#Written by Malachi Griffith
#Copyright 2009 Malachi Griffith
#This file is part of 'ALEXA-Seq'
#ALEXA-Seq is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ALEXA-Seq is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ALEXA-Seq (COPYING.txt).  If not, see <http://www.gnu.org/licenses/>.

#The purpose of this script is to divide mapped reads into separate files representing genome partitions
#This will be done for reads that map within transcripts, junctions, boundaries, introns and intergenic regions

#Input:
#Partition Coords File, Analysis Dir, Library ID, Read Assignment Class

#Steps:
#1.) Check and form neccessary paths
#2.) Import genome partition coordinates
#3.) Build a hash of read IDs that have been assigned to the Read Class Assignment specified -> get from read record files
#4.) Delete all partition files for this read class to ensure they are created cleanly
#5.) Write header to each partition file to be created (one for each region)
#    - Files should be named: COMBINED_Lane0_X_Y_ENST_v53.txt.gz (Where X and Y are chromosome and region number)
#6.) Go through the mapping results file and identify reads that have been assigned to the current class
#    - If only one read of a read pair has been assigned to the current class, change it's status to 'Ambiguous'
#    - Determine which chromosome and region each read corresponds to and write the record to the corresponding partition file
#    - If read1 and read2 correspond to different regions, write the record to both partition files
#7.) Compress all results files

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use File::Basename;
use BerkeleyDB;
use IO::File;

#Load the ALEXA modules
BEGIN {
  use Cwd 'abs_path';
  if (abs_path($0) =~ /(.*)\/.*\/.*\.pl/){
    push (@INC, $1);
  }
}
use utilities::utility qw(:all);

my $analysis_dir = '';
my $library_id = '';
my $genome_regions_file = '';
my $ensembl_version = '';
my $read_class = '';
my $seq_size = '';
my $annotation_dir = '';

GetOptions ('analysis_dir=s'=>\$analysis_dir,'library_id=s'=>\$library_id, 'genome_regions_file=s'=>\$genome_regions_file, 'ensembl_version=i'=>\$ensembl_version, 'read_class=s'=>\$read_class,
            'seq_size=i'=>\$seq_size, 'annotation_dir=s'=>\$annotation_dir);

print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the analysis dir using: --analysis_dir", RESET;
print GREEN, "\n\tSpecify the annotation dir using: --annotation_dir", RESET;
print GREEN, "\n\tSpecify the junction/boundary size for this build using: --seq_size", RESET;
print GREEN, "\n\tSpecify the library ID using: --library_id", RESET;
print GREEN, "\n\tSpecify the genome regions file using: --genome_regions_file", RESET;
print GREEN, "\n\tSpecify the ensembl version number used for the analysis using: --ensembl_version", RESET;
print GREEN, "\n\tSpecify the read class to process using: --read_class (ENST NOVEL_JUNCTION NOVEL_BOUNDARY INTRON INTERGENIC)", RESET;
print GREEN, "\n\nExample: partitionMapFiles.pl  --analysis_dir=/projects/malachig/alexa_seq/  --annotation_dir=/projects/alexa/sequence_databases/hs_53_36o/  --seq_size=62  --library_id=HS1361   --genome_regions_file=/projects/alexa/sequence_databases/hs_53_36o/Regions_250_Genes.txt  --ensembl_version=53  --read_class=ENST\n\n", RESET;

unless ($analysis_dir && $annotation_dir && $seq_size && $library_id && $genome_regions_file && $ensembl_version && $read_class){
  print RED, "\nRequired input parameter(s) missing\n\n", RESET;
  exit();
}

unless ($read_class =~ /^ENST$|^NOVEL\_JUNCTION$|^NOVEL\_BOUNDARY$|^INTRON$|^INTERGENIC$/){
  print RED, "\nRead class (--read_class) must be one of: (ENST NOVEL_JUNCTION NOVEL_BOUNDARY INTRON INTERGENIC)\n\n", RESET;
  exit();
}

#Define class to dirname mappings - note that for junctions and boundaries, hits to known junctions and boundaries will have been assigned to ENST, so we need to grab those when summarizing the junctions/boundaries for expression
my %class;
$class{'ENST'}{dir_name} = "ENST_v$ensembl_version";
$class{'ENST'}{test_class} = "ENST_U";
$class{'NOVEL_JUNCTION'}{dir_name} = "Junctions_v$ensembl_version";
$class{'NOVEL_JUNCTION'}{test_class} = "NOVEL_JUNCTION_U ENST_U";
$class{'NOVEL_BOUNDARY'}{dir_name} = "Boundaries_v$ensembl_version";
$class{'NOVEL_BOUNDARY'}{test_class} = "NOVEL_BOUNDARY_U ENST_U";
$class{'INTRON'}{dir_name} = "Introns_v$ensembl_version";
$class{'INTRON'}{test_class} = "INTRON_U";
$class{'INTERGENIC'}{dir_name} = "Intergenics_v$ensembl_version";
$class{'INTERGENIC'}{test_class} = "INTERGENIC_U";

#1.) Check and form neccessary paths
$analysis_dir = &checkDir('-dir'=>$analysis_dir, '-clear'=>"no");
$annotation_dir = &checkDir('-dir'=>$annotation_dir, '-clear'=>"no");

unless (-e $genome_regions_file){
  print RED, "\n\nCould not find genome regions file: $genome_regions_file - check path\n\n", RESET;
  exit();
}
my $read_records_dir = "$analysis_dir"."read_records/$library_id/";
$read_records_dir = &checkDir('-dir'=>$read_records_dir, '-clear'=>"no");
my $map_files_dir = "$read_records_dir"."$class{$read_class}{dir_name}/";
$map_files_dir = &checkDir('-dir'=>$map_files_dir, '-clear'=>"no");
my $part_dir = &createNewDir('-path'=>$map_files_dir, '-new_dir_name'=>"part", '-force'=>"yes");
my $working_dir = "$map_files_dir"."temp/";
$working_dir = &checkDir('-dir'=>$working_dir, '-clear'=>"no");

#Open a log file
my $logfile = "$analysis_dir"."logs/$library_id/partitionMapFiles/partitionMapFiles_"."$read_class"."_LOG.txt";
open(LOG, ">$logfile") || die "\n\nCould not open log file: $logfile\n\n";

print BLUE, "\n\nPartitioning map files for library: $library_id and read class: $read_class", RESET;
print LOG "\n\nPartitioning map files for library: $library_id and read class: $read_class";

#2-a.) Import genome partition coordinates and create a sub directory for each region
my $regions_ref = &importRegions('-regions_file'=>$genome_regions_file, '-part_dir'=>$part_dir);
my $mem_message = &memoryUsage();
print YELLOW, "\n$mem_message", RESET;
print LOG "\n$mem_message";


#2-b.) Import features annotations for this class and determine which genome region each corresponds to... 
#Organize the genome regions by chromosome
#For each Gene, Junction, Boundary, Intron or Intergenic ... Get the chromosome, start and end coordinates from the annotation file
#Figure out which region the feature corresponds to and store only the feature ID and region ID
my $features_ref = &importFeatures('-annotation_dir'=>$annotation_dir, '-seq_size'=>$seq_size, '-regions_ref'=>$regions_ref);
$mem_message = &memoryUsage();
print YELLOW, "\n$mem_message", RESET;
print LOG "\n$mem_message";


#3.) Delete all partition files for this read class to ensure they are created cleanly
# - Not neccessary as the 'part' directory was created cleanly...


#4.) Write header to each partition file to be created (one for each region)
#    - Get the map files to be processed
#    - Get the column positions for each map file
#    - Grab the headers from each lane file (make sure they are all the same)
#    - Initialize each map region file by writing the header to it
#    - Files should be named: COMBINED_Lane0_X_Y_ENST_v53.txt.gz (Where X and Y are chromosome and region number)
my $common_suffix = '';
my $header_line = '';
my $map_files_ref = &getMapFiles('-map_files_dir'=>$map_files_dir, '-regions_ref'=>$regions_ref);
$mem_message = &memoryUsage();
print YELLOW, "\n$mem_message", RESET;
print LOG "\n$mem_message";


#5.) Build a hash of read IDs that have been assigned to the Read Class Assignment specified -> get from read record files
my %reads;
my $reads_ref = \%reads;
my $reads_db_file = "$working_dir"."$library_id"."$read_class"."_ReadPairs.btree";
system ("rm -f $reads_db_file");
tie(%reads, 'BerkeleyDB::Btree', -Cachesize => 256000000, -Filename=> $reads_db_file, -Flags => DB_CREATE) or die "can't open file $reads_db_file: $! $BerkeleyDB::Error\n";
my $test_class = $class{$read_class}{test_class};
&importTargetReads('-read_records_dir'=>$read_records_dir, '-read_class'=>$test_class);
$mem_message = &memoryUsage();
print YELLOW, "\n$mem_message", RESET;
print LOG "\n$mem_message";


#6.) Go through the mapping results file and identify reads that have been assigned to the current class
#    - If only one read of a read pair has been assigned to the current class, change it's status to 'Ambiguous'
#    - Determine which chromosome and region each read corresponds to and write the record to the corresponding partition file
#    - If read1 and read2 correspond to different regions, write the record to both partition files
&parseMapFiles('-map_files_ref'=>$map_files_ref, '-features_ref'=>$features_ref, '-reads_ref'=>$reads_ref);
$mem_message = &memoryUsage();
print YELLOW, "\n$mem_message", RESET;
print LOG "\n$mem_message";

#Untie and delete Berkeley DB files
untie(%reads);
system ("rm -f $reads_db_file");


#7.) Close all open file handles and compress all results files
print BLUE, "\n\nCompressing genome region map files:", RESET;
print LOG "\n\nCompressing genome region map files:";

foreach my $r (sort {$a <=> $b} keys %{$regions_ref}){
  my $fh = $regions_ref->{$r}->{region_map_file_handle};
  $fh->close;
  my $region_map_file_path = $regions_ref->{$r}->{region_map_file_path};
  my $cmd_gz = "gzip $region_map_file_path";
  print BLUE, "\n\t$cmd_gz", RESET;
  system($cmd_gz);
}

$mem_message = &memoryUsage();
print YELLOW, "\n$mem_message", RESET;
print "\n\nSCRIPT COMPLETE\n\n";
print LOG "\n$mem_message";
print LOG "\n\nSCRIPT COMPLETE\n\n";

close(LOG);

exit();



####################################################################################################################################
#Import the genome regions file                                                                                                    #
####################################################################################################################################
sub importRegions{
  my %args = @_;
  my $regions_file = $args{'-regions_file'};
  my $part_dir = $args{'-part_dir'};

  print BLUE, "\n\nImporting genome region coordinates", RESET;
  print LOG "\n\nImporting genome region coordinates";

  my %regions;
  open(REGIONS, "$regions_file") || die "\n\nCould not open region file: $regions_file\n\n";
  my $header = 1;
  my $r = 0;
  while(<REGIONS>){
    if ($header == 1){$header = 0; next();}
    $r++;
    chomp($_);
    my @line = split("\t", $_);
    $regions{$r}{chromosome} = $line[0];
    $regions{$r}{region} = $line[1];
    $regions{$r}{start_chr} = $line[2];
    $regions{$r}{end_chr} = $line[3];
    $regions{$r}{size} = $line[4];
    $regions{$r}{gene_count} = $line[5];
  }
  close(REGIONS);

  #Create sub dir for each region within the partition dir
  foreach my $r (sort {$a <=> $b} keys %regions){
    my $sub_dir = "$part_dir"."$regions{$r}{chromosome}"."_"."$regions{$r}{region}/";
    mkdir($sub_dir);
    $regions{$r}{sub_dir} = $sub_dir;
  }

  return(\%regions);
}


####################################################################################################################################
#Import the annotation for the specified class                                                                                     #
####################################################################################################################################
sub importFeatures{
  my %args = @_;
  my $annotation_dir = $args{'-annotation_dir'};
  my $seq_size = $args{'-seq_size'};
  my $regions_ref = $args{'-regions_ref'};

  print BLUE, "\n\nImporting feature annotations and mapping each to a genome region partition", RESET;
  print LOG "\n\nImporting feature annotations and mapping each to a genome region partition";

  #Organize the genome regions by chromosome
  my %regions_chr;
  foreach my $r (sort {$a <=> $b} keys %{$regions_ref}){
    my $chr = $regions_ref->{$r}->{chromosome};
    if ($regions_chr{$chr}){
      my $r_ref = $regions_chr{$chr}{list};
      $r_ref->{$r}=1;
    }else{
      my %temp;
      $temp{$r}=1;
      $regions_chr{$chr}{list} = \%temp;
    }
  }

  #Import features annotations for this class and determine which genome region each corresponds to... 
  #For each Gene, Junction, Boundary, Intron or Intergenic ... Get the chromosome, start and end coordinates from the annotation file
  #Figure out which region the feature corresponds to and store only the feature ID and region ID

  my %annotation_files;
  my $chromosome_col_name = "Chromosome";
  my @coord_col_names;
  $annotation_files{'ENST'}{path} = "$annotation_dir"."genes/genes_annotated.txt.gz";
  @coord_col_names = qw(Unit1_start_chr Unit1_end_chr);
  $annotation_files{'ENST'}{coord_col_names} = \@coord_col_names;

  $annotation_files{'NOVEL_JUNCTION'}{path} = "$annotation_dir"."exonJunctions/exonJunctions_"."$seq_size"."mers_annotated.txt.gz";
  @coord_col_names = qw(Unit1_start_chr Unit1_end_chr Unit2_start_chr Unit2_end_chr);
  $annotation_files{'NOVEL_JUNCTION'}{coord_col_names} = \@coord_col_names;

  $annotation_files{'NOVEL_BOUNDARY'}{path} = "$annotation_dir"."exonBoundaries/exonBoundaries_"."$seq_size"."mers_annotated.txt.gz";
  @coord_col_names = qw(Unit1_start_chr Unit1_end_chr);
  $annotation_files{'NOVEL_BOUNDARY'}{coord_col_names} = \@coord_col_names;

  $annotation_files{'INTRON'}{path} = "$annotation_dir"."introns/introns_annotated.txt.gz";
  @coord_col_names = qw(Unit1_start_chr Unit1_end_chr);
  $annotation_files{'INTRON'}{coord_col_names} = \@coord_col_names;

  $annotation_files{'INTERGENIC'}{path} = "$annotation_dir"."intergenics/intergenics_annotated.txt.gz";
  @coord_col_names = qw(Unit1_start_chr Unit1_end_chr);
  $annotation_files{'INTERGENIC'}{coord_col_names} = \@coord_col_names;

  my %features;
  my $annotation_file = $annotation_files{$read_class}{path};
  @coord_col_names = @{$annotation_files{$read_class}{coord_col_names}};

  my $header = 1;
  my %columns;
  my $features_checked = 0;
  my $features_mapped = 0;
  open (ANN, "zcat $annotation_file |") || die "\n\nCould not open annotation file: $annotation_file\n\n";
  while (<ANN>){
    chomp($_);
    my @line = split("\t", $_);
    if ($header == 1){
      my $column_count = 0;
      foreach my $column (@line){
        $columns{$column}{column_pos} = $column_count;
        $column_count++;
      }
      $header = 0;
      next();
    }
    $features_checked++;
    my @coords;
    foreach my $col_name (@coord_col_names){
      my $coord = $line[($columns{$col_name}{column_pos})];
      push(@coords, $coord);
    }
    my @coords_sort = sort {$a <=> $b} @coords;
    my $lower_coord = $coords_sort[0];
    my $upper_coord = $coords_sort[(scalar(@coords_sort)-1)];

    my $feature_id = $line[0];
    my $chromosome = $line[($columns{$chromosome_col_name}{column_pos})];

    if ($chromosome eq "M"){
      $chromosome = "MT";
    }

    #Which region does this feature correspond to?
    unless($regions_chr{$chromosome}){
      print RED, "\n\nCould not find regions list for chromosome: $chromosome", RESET;
      print LOG "\n\nCould not find regions list for chromosome: $chromosome";
      exit();
    }
    my $r_list_ref = $regions_chr{$chromosome}{list};

    my $found_feature = 0;
    foreach my $r (keys %{$r_list_ref}){
      my $r_chr_start = $regions_ref->{$r}->{start_chr};
      my $r_chr_end = $regions_ref->{$r}->{end_chr};

      if ($lower_coord >= $r_chr_start && $lower_coord <= $r_chr_end && $upper_coord >= $r_chr_start && $upper_coord <= $r_chr_end){
        $features{$feature_id}=$r;
        $features_mapped++;
        $found_feature = 1;
        last();
      }
    }
    unless($found_feature){
      #print YELLOW, "\n\nFailed to assign a feature to a region.  Feature: $chromosome: $lower_coord - $upper_coord", RESET;
      foreach my $r (keys %{$r_list_ref}){
        my $r_chr_start = $regions_ref->{$r}->{start_chr};
        my $r_chr_end = $regions_ref->{$r}->{end_chr};

        if (($lower_coord >= $r_chr_start && $lower_coord <= $r_chr_end) || ($upper_coord >= $r_chr_start && $upper_coord <= $r_chr_end)){
          $features{$feature_id}=$r;
          $features_mapped++;
          $found_feature = 1;
          last();
        }
      }
    }
  }
  close(ANN);

  if ($features_checked == $features_mapped){
    print BLUE, "\n\tChecked $features_checked features and successfully mapped $features_mapped to a genome region", RESET;
    print LOG "\n\tChecked $features_checked features and successfully mapped $features_mapped to a genome region";
  }else{
    print RED, "\n\nFailed to assign ALL features to a genome region!! - (Checked $features_checked features and successfully mapped $features_mapped) - Aborting", RESET;
    print LOG "\n\nFailed to assign ALL features to a genome region!! - (Checked $features_checked features and successfully mapped $features_mapped) - Aborting";
    exit();
  }
  return(\%features);
}


####################################################################################################################################
#Import the target read IDs for this class                                                                                         #
####################################################################################################################################
sub importTargetReads{
  my %args = @_;
  my $read_records_dir = $args{'-read_records_dir'};
  my $read_class_u = $args{'-read_class'};

  print BLUE, "\n\nImporting read IDs assigned to the class: $read_class_u for the library: $library_id", RESET;
  print LOG "\n\nImporting read IDs assigned to the class: $read_class_u for the library: $library_id";

  my %files;
  opendir(DIRHANDLE, "$read_records_dir") || die "\nCannot open directory: $read_records_dir\n\n";
  my @test_files = readdir(DIRHANDLE);
  my $file_count = 0;
  my $stored_reads = 0;

  foreach my $test_file (sort @test_files){
    my $file_path = "$read_records_dir"."$test_file";

    #Skip directories within the specified directory
    if (-e $file_path && -d $file_path){
      next();
    }
    #If the results file is compressed uncompress it
    unless ($file_path =~ /(.*)\.gz$/){
      print RED, "\nFound an uncompressed file: $file_path - make sure no files in progress are in this directory!!\n\n", RESET;
      print LOG "\nFound an uncompressed file: $file_path - make sure no files in progress are in this directory!!\n\n";
      exit();
    }
    $file_count++;
    print BLUE, "\n\t$file_path was added to the list of files to be processed", RESET;
    print LOG "\n\t$file_path was added to the list of files to be processed";
    
    $files{$file_count}{name} = $test_file;
    $files{$file_count}{path} = $file_path;
  }
  closedir(DIRHANDLE);

  #Get reads matching the current class for all lanes of data
  print BLUE, "\n\nProcessing these files and storing read IDs matching the target class", RESET;
  print LOG "\n\nProcessing these files and storing read IDs matching the target class";
  
  foreach my $fc (sort {$a <=> $b} keys %files){
    my $file_name = $files{$fc}{name};
    my $file_path = $files{$fc}{path};
    my %columns;
    print BLUE, "\n\nProcessing: $file_name\n", RESET;
    print LOG "\n\nProcessing: $file_name\n";

    open (READ, "zcat $file_path |") || die "\nCould not open file: $file_path\n\n";
    my $header = 1;
    my $counter = 0;
    while(<READ>){
      $counter++;
      chomp($_);
      my @line = split("\t", $_);
      if ($header == 1){
        my $column_count = 0;
        foreach my $column (@line){
          $columns{$column}{column_pos} = $column_count;
          $column_count++;
        }
        $header = 0;
        next();
      }

      if ($counter == 100000){
        $counter = 0;
        $| = 1; print BLUE, ".", RESET; $| = 0;
      }
      my $r1_id = $line[$columns{'Read1_ID'}{column_pos}];
      my $r2_id = $line[$columns{'Read2_ID'}{column_pos}];
      my $r1_status = $line[$columns{'Read1_Status'}{column_pos}];
      my $r2_status = $line[$columns{'Read2_Status'}{column_pos}];

      if ($read_class_u =~ /($r1_status)/){
        $reads{$r1_id}=1;
        $stored_reads++;
      }
      if ($read_class_u =~ /($r2_status)/){
        $reads{$r2_id}=1;
        $stored_reads++;
      }
    }
    close(READ);
    my $mem_message = &memoryUsage();
    print YELLOW, "\n$mem_message", RESET;
    print LOG "\n$mem_message";
  }

  print BLUE, "\n\nStored a grand total of $stored_reads reads", RESET;
  print LOG "\n\nStored a grand total of $stored_reads reads";

  return();
}


####################################################################################################################################
#Import the mapping results files                                                                                                  #
####################################################################################################################################
sub getMapFiles{
  my %args = @_;
  my $map_files_dir = $args{'-map_files_dir'};
  my $regions_ref = $args{'-regions_ref'};

  print BLUE, "\n\nGet map files to be processed and initializing genome region partition files", RESET;
  print LOG "\n\nGet map files to be processed and initializing genome region partition files";

  #Get the map files to be processed
  my %map_files;
  opendir(DIRHANDLE, "$map_files_dir") || die "\nCannot open directory: $map_files_dir\n\n";
  my @test_files = readdir(DIRHANDLE);
  my $file_count = 0;
  my %suffix_list;
  foreach my $test_file (sort @test_files){
    my $file_path = "$map_files_dir"."$test_file";

    #Skip directories within the specified directory
    if (-e $file_path && -d $file_path){
      next();
    }
    #If the results file is compressed uncompress it
    unless ($file_path =~ /(.*)\.gz$/){
      print RED, "\nFound an uncompressed file: $file_path - make sure no files in progress are in this directory!!\n\n", RESET;
      print LOG "\nFound an uncompressed file: $file_path - make sure no files in progress are in this directory!!\n\n";
      exit();
    }
    my $suffix = '';
    if ($test_file =~ /Lane\d+\_(.*)\.txt\.gz/){
      $suffix = $1;
    }
    $common_suffix = $suffix;
    $suffix_list{$suffix} = 1;

    $file_count++;
    print BLUE, "\n\t$file_path was added to the list of files to be processed", RESET;
    print LOG "\n\t$file_path was added to the list of files to be processed";
    $map_files{$file_count}{name} = $test_file;
    $map_files{$file_count}{path} = $file_path;
    $map_files{$file_count}{suffix} = $suffix;
  }
  closedir(DIRHANDLE);

  my $unique_suffixes = keys %suffix_list;
  if ($unique_suffixes > 1){
    print RED, "\n\nFound more than one distinct suffix for the map files!!!  Aborting...", RESET;
    print LOG "\n\nFound more than one distinct suffix for the map files!!!  Aborting...";
    exit();
  }

  #Get the column positions for each map file
  #Grab the headers from each lane file (make sure they are all the same)
  my %header_list;
  foreach my $fc (sort {$a <=> $b} keys %map_files){
    my %columns;
    my $path = $map_files{$fc}{path};
    open (IN, "zcat $path | head -n 1 | ") || die "\n\nCould not open map file for header check\n\n";
    while(<IN>){
      chomp($_);
      $header_line = $_;
      $header_list{$_}=1;
      my $col_pos = 0;
      my @line = split("\t", $_);
      foreach my $value (@line){
        $columns{$value}{column_pos} = $col_pos;
        $col_pos++;
      }
    }
    close(IN);
    $map_files{$fc}{columns} = \%columns;
    $map_files{$fc}{header_line} = $header_line;
  }
  my $unique_headers = keys %header_list;
  if ($unique_headers > 1){
    print RED, "\n\nFound more than one distinct header line for the map files!!!  Aborting...", RESET;
    print LOG "\n\nFound more than one distinct header line for the map files!!!  Aborting...";
    exit();
  }

  #Initialize each map region file by writing the header to it
  #Files should be named: COMBINED_Lane0_X_Y_ENST_v53.txt.gz (Where X and Y are chromosome and region number)
  foreach my $r (sort {$a <=> $b} keys %{$regions_ref}){
    my $chromosome = $regions_ref->{$r}->{chromosome};
    my $region_number = $regions_ref->{$r}->{region};
    my $sub_dir = $regions_ref->{$r}->{sub_dir};

    my $region_map_file_name = "COMBINED_Lane0_"."$chromosome"."_"."$region_number"."_"."$common_suffix".".txt";
    my $region_map_file_path = "$sub_dir"."$region_map_file_name";
    $regions_ref->{$r}->{region_map_file_name} = $region_map_file_name;
    $regions_ref->{$r}->{region_map_file_path} = $region_map_file_path;
    my $fh = IO::File->new(">$region_map_file_path") || die "\nCould not open output file handle: $region_map_file_path\n\n";

    if (defined($fh)){
      #Do nothing
    }else{
      print YELLOW, "\n\nCould not open file handle, sleeping and attempting one more time", RESET;
      print LOG "\n\nCould not open file handle, sleeping and attempting one more time";
      system("sleep 10");
      $fh = IO::File->new(">$region_map_file_path") || die "\nCould not open output file handle: $region_map_file_path\n\n";
    }

    unless(defined($fh)){
      print RED, "\n\ndid not retrieve valid file handle for: $region_map_file_path\n\n", RESET;
      print LOG "\n\ndid not retrieve valid file handle for: $region_map_file_path\n\n";
      exit();
    }
    $regions_ref->{$r}->{region_map_file_handle} = $fh; 
    print $fh "$header_line\n";
  }

  return(\%map_files);
}


####################################################################################################################################
#Parse the map files into genome region partitions                                                                                 #
####################################################################################################################################
sub parseMapFiles{
  my %args = @_;
  my $map_files_ref = $args{'-map_files_ref'};
  my $features_ref = $args{'-features_ref'};
  my $reads_ref = $args{'-reads_ref'};

  print BLUE, "\n\nParse map files and divide into genome regions for class: $read_class", RESET;
  print LOG "\n\nParse map files and divide into genome regions for class: $read_class";

  #Go through the mapping results file and identify reads that have been assigned to the current class
  #If only one read of a read pair has been assigned to the current class, change the other read's status to 'Ambiguous'
  #Determine which chromosome and region each read corresponds to and write the record to the corresponding partition file
  #If read1 and read2 correspond to different regions, write the record to both partition files

  #Define the column containing start/end coordinates for each class
  my $r1_id_col_name;
  my $r2_id_col_name;
  my $r1_feature_id_col_name;
  my $r2_feature_id_col_name;
  my $r1_hit_type_col_name;
  my $r2_hit_type_col_name;
  my $paired;

  if ($read_class eq "ENST"){
    $paired = 1;
    $r1_id_col_name = "R1_ID";
    $r2_id_col_name = "R2_ID";
    $r1_feature_id_col_name = "R1_GeneID";
    $r2_feature_id_col_name = "R2_GeneID";
    $r1_hit_type_col_name = "R1_HitType";
    $r2_hit_type_col_name = "R2_HitType";
  }elsif($read_class eq "NOVEL_JUNCTION"){
    $paired = 0;
    $r1_id_col_name = "Read_ID";
    $r1_feature_id_col_name = "Junction_ID";
    $r1_hit_type_col_name = "HitType";
  }elsif($read_class eq "NOVEL_BOUNDARY"){
    $paired = 0;
    $r1_id_col_name = "Read_ID";
    $r1_feature_id_col_name = "Boundary_ID";
    $r1_hit_type_col_name = "HitType";
  }elsif($read_class eq "INTRON"){
    $paired = 1;
    $r1_id_col_name = "R1_ID";
    $r2_id_col_name = "R2_ID";
    $r1_feature_id_col_name = "R1_IntronName";
    $r2_feature_id_col_name = "R2_IntronName";
    $r1_hit_type_col_name = "R1_HitType";
    $r2_hit_type_col_name = "R2_HitType";
  }elsif($read_class eq "INTERGENIC"){
    $paired = 1;
    $r1_id_col_name = "R1_ID";
    $r2_id_col_name = "R2_ID";
    $r1_feature_id_col_name = "R1_IntergenicName";
    $r2_feature_id_col_name = "R2_IntergenicName";
    $r1_hit_type_col_name = "R1_HitType";
    $r2_hit_type_col_name = "R2_HitType";
  }else{
    print RED, "\n\nRead class: $read_class not understood\n\n", RESET;
    print LOG "\n\nRead class: $read_class not understood\n\n";
  }


  foreach my $fc (sort {$a <=> $b} keys %{$map_files_ref}){
    my %columns = %{$map_files_ref->{$fc}->{columns}};
    my $file_path = $map_files_ref->{$fc}->{path};

    print BLUE, "\n\nProcessing: $file_path", RESET;
    print LOG "\n\nProcessing: $file_path";
    open(IN, "zcat $file_path | ") || die "\n\nCould not open map file: $file_path\n\n";
    my $header = 1;
    my $counter = 0;
    while(<IN>){
      if ($header == 1){
        $header = 0;
        next();
      }
      chomp($_);
      my @line = split("\t", $_);

      if ($counter == 100000){
        $counter = 0;
        $| = 1; print BLUE, ".", RESET; $| = 0;
      }

      #READ1
      my $r1_test = 0;
      my $r1_id = $line[$columns{$r1_id_col_name}{column_pos}];
      my $r1_feature_id = $line[$columns{$r1_feature_id_col_name}{column_pos}];
      my $r1_hit_type = $line[$columns{$r1_hit_type_col_name}{column_pos}];
      my $r1_region;
      my $r1_region_fh;

      if (($reads_ref->{$r1_id}) && ($r1_hit_type eq "Top_Hit")){
        $r1_test = 1;
        $r1_region = $features_ref->{$r1_feature_id};
        $r1_region_fh = $regions_ref->{$r1_region}->{region_map_file_handle};
        #print YELLOW, "\n\nMATCH:\nLINE: @line\nr1_test: $r1_test\tr1_id: $r1_id\tr1_feature_id: $r1_feature_id\tr1_region: $r1_region", RESET;

      }else{
        #print YELLOW, "\n\nNO-MATCH:\nLINE: @line\nr1_test: $r1_test\tr1_id: $r1_id\tr1_feature_id: $r1_feature_id", RESET;
        $line[$columns{$r1_hit_type_col_name}{column_pos}] = "Ambiguous";
      }

      #READ2
      my $r2_test = 0;
      my $r2_id = '';
      my $r2_feature_id = '';
      my $r2_hit_type;
      my $r2_region;
      my $r2_region_fh;

      if ($paired == 1){
        $r2_id = $line[$columns{$r2_id_col_name}{column_pos}];
        $r2_feature_id = $line[$columns{$r2_feature_id_col_name}{column_pos}];
        $r2_hit_type = $line[$columns{$r2_hit_type_col_name}{column_pos}];

        if (($reads_ref->{$r2_id}) && ($r2_hit_type eq "Top_Hit")){
          $r2_test = 1;
          $r2_region = $features_ref->{$r2_feature_id};
          $r2_region_fh = $regions_ref->{$r2_region}->{region_map_file_handle};
        }else{
          $line[$columns{$r2_hit_type_col_name}{column_pos}] = "Ambiguous";
        }
      }

      #If neither read is a member of the current class, skip
      unless ($r1_test || $r2_test){
        next();
      }
      my $tmp_seperator = $";
      $" = "\t";

      if ($paired == 0){
        #Only Read1 records
        if ($r1_test){
          print $r1_region_fh "@line\n";
        }

      }else{
        #Paired Read1+Read2 records
        if ($r1_test && $r2_test){
          if ($r1_region == $r2_region){
            print $r1_region_fh "@line\n";
          }else{
            print $r1_region_fh "@line\n";
            print $r2_region_fh "@line\n";
          }
        }elsif($r1_test){
          print $r1_region_fh "@line\n";
        }elsif($r2_test){
          print $r2_region_fh "@line\n";
        }
      }
      $" = $tmp_seperator;
    }
    close(IN);
  }
  return();
}

