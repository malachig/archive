#!/usr/bin/perl -w

#Written by Malachi Griffith
#This script takes a file path as input and attempts to concatenate all probe files present in that directory
#Multiple probe files arise for various reasons in the ArrayDesign pipeline.  In particular from jobs done on the cluster 
#where a single probe file is divided into multiple smaller pieces.  I have also seperated the process of extracting different
#probe types (exon-exon versus exon-intron versus exon-internal)

#During this process basic quality checks will be conducted.  Multiple header lines will be removed but will also be checked to 
#ensure that you are not trying to join files with different numbers of columns, etc.
#It will also make sure that the probes are ordered by their design ID number and ensure that no duplicate IDs exist in any of the files

#Two options from this point:
#   A.) Merge the data collected onto a master probe file from which they were originally generated
#       - For example, say we have files with (Probe_Count, simfold_score, pairfold_score)
#       - We can append these probe records to those in the specified master probe file
#   B.) Write out a joined file based on only the files in the specified directory

#NOTE1: All probe files must be tab-delimited and contain numeric unique probe IDs in the first column
#NOTE2: Header lines must start with "Probe_Count"
#NOTE3: Because of the sanity checks this script may be very memory intensive!

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);

my $source_directory = '';
my $probes_expected = '';
my $master_probe_file = '';
my $output_file = '';

GetOptions ('source_dir=s'=>\$source_directory, 'probes_expected=i'=>\$probes_expected, 'master_probe_file=s'=>\$master_probe_file, 'output_file=s'=>\$output_file);

#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify a directory containing the tab-delimited input probe files to be joined using: --source_dir", RESET;
print GREEN, "\n\tSpecify the total number of probes expected in all of these files using: --probes_expected", RESET;
print GREEN, "\n\tSpecify the out probe file name using: --output_file", RESET;
print GREEN, "\n\n\tIF the probe files being joined contain only partial probe info, merge them with a master probe file using:  --master_probe_file", RESET;
print GREEN, "\n\t\tIF the probe files being joined already contain all columns, skip this option", RESET;
print GREEN, "\n\nExample: joinProbeFiles  --source_dir=/home/user/alexa/ALEXA_version/fold/results  --probes_expected=12345  --master_probe_file=/home/user/alexa/ALEXA_version/unfiltered_probes/exonJunctionProbes.txt  --output_file=/home/user/alexa/ALEXA_version/unfiltered_probes/exonJunctionProbes_joined.txt\n\n", RESET;

unless ($source_directory && $probes_expected && $output_file){
  print RED, "\nRequired input parameter(s) missing!\n\n", RESET;
  exit();
}

#First make sure the source path exists and is a directory
unless ($source_directory =~ /.*\/$/){
  $source_directory = "$source_directory"."/";
}

unless (-e $source_directory && -d $source_directory){
  print RED, "\nSpecified directory: $source_directory does not appear valid!\n\n", RESET;
  exit();
}


#First get all files from the specified directory
print BLUE, "\nSearching $source_directory for probe files", RESET;

my @files;
opendir(DIRHANDLE, "$source_directory") || die "\nCannot open directory: $source_directory\n\n";
my @test_files = readdir(DIRHANDLE);

foreach my $test_file (@test_files){
  my $file_path = "$source_directory"."$test_file";

    #Skip directories within the specified directory
    if (-e $file_path && -d $file_path){
      print "\n\t$file_path  is a directory - skipping";
      next();
    }
  push(@files, $file_path);
}

my $file_num = scalar(@files);

print BLUE, "\n\nFound $file_num probe files in the specified directory\n", RESET;

#Check the numbers of each file
my %file_nums;
my $file_count = 0;
foreach my $file (@files){
  $file_count++;
  $file_nums{$file_count}{file_name} = $file;
  if ($file =~ /(\d+)$/){
    $file_nums{$file_count}{file_number} = $1;
  }else{
    $file_nums{$file_count}{file_number} = $file_count;
  }
}

#Display the to be processed files to the user in order
print BLUE, "\nThe files to be processed are as follows:", RESET;
my $loop_count = 0;
foreach my $file (sort {$file_nums{$a}->{file_number} <=> $file_nums{$b}->{file_number}} keys %file_nums){
  $loop_count++;
  print YELLOW, "\n\t$loop_count\tFile: $file_nums{$file}{file_number}\tName: $file_nums{$file}{file_name}", RESET;
}

my %probes;
my $master_probe_count = 0;

#Grab records from these files
my $header_line;
print BLUE, "\n\nBegin processing $file_count files:", RESET;

foreach my $file (sort {$file_nums{$a}->{file_number} <=> $file_nums{$b}->{file_number}} keys %file_nums){

  my $probe_file = $file_nums{$file}{file_name};
  chomp ($probe_file);

  print BLUE, "\n\n\tProcessing $probe_file", RESET;

  open (PROBE_FILE, "$probe_file") || die "\nCould not open probe file: $probe_file\n\n";

  my $probe_count = 0;

  while (<PROBE_FILE>){
    my $line_record = $_;
    chomp ($line_record);

    my @line = split ("\t", $line_record);

    #Watch for the header line which is assumed to start with "Probe_Count"
    if ($line[0] =~ /^Probe_Count$/){
      my $new_header = $line_record;

      #Make sure this header is the same as any encountered before
      if ($header_line){
	unless ($new_header eq $header_line){
	  print RED, "\nYou are attempting to merge files with different header lines!  Bad user!\n\n", RESET;
	  exit();
	}
      }else{
	$header_line = $new_header;
      }
      next();
    }

    #If the line is a data line - make sure it starts with a probe ID
    if ($line[0] =~ /(\d+)/){
      my $probe_id = $1;

      #Check for duplicate probe count IDs
      if ($probes{$probe_id}){
	print RED, "\nFound a probe ID that has already been used previously: $probe_id !!\n\n", RESET;
	exit();	
      }else{
	$probes{$probe_id}{line_record} = $line_record;
	$probe_count++;
	$master_probe_count++;
      }
    }
  }
  close PROBE_FILE;
  print BLUE, "\n\tFound $probe_count probes in this file", RESET;
}

print BLUE, "\n\nFound a grand total of $master_probe_count probes\n\n", RESET;

#Make sure this is the correct number according to those expected by the user
unless ($master_probe_count == $probes_expected){
  print RED, "\nDid not find the expected ($probes_expected) number of probes!!", RESET;
  print RED, "\nCheck files in input directory for problems\n\n", RESET;
  exit();
}

#Two options from this point:
#   A.) Merge the data collected onto a master probe file from which they were originally generated
#       - For example, say we have files with (Probe_Count, simfold_score, pairfold_score)
#       - We can append these probe records to those in the specified master probe file
#   B.) Write out a joined file based on only the files in the specified directory

#Now write all probes out to a single master probe file as specified
open (NEW_PROBE_FILE, ">$output_file") || die "\nCould not open new probe file: $output_file\n\n";

if ($master_probe_file){
  #A.) Merge to master probe file
  open (MASTER_PROBE_FILE, "$master_probe_file") || die "\nCould not open master probe file: $master_probe_file\n\n";

  #get the master header line
  my $master_header_line;
  my $first_line = 1;
  my $master_probes_found = 0;

  while(<MASTER_PROBE_FILE>){
    chomp($_);
    my $current_line_record = $_;
    my @line = split ("\t", $_);

    if ($first_line == 1){
      $master_header_line = $_;
      $first_line = 0;

      #Watch for the header line which is assumed to start with "Probe_Count"
      unless ($line[0] =~ /^Probe_Count$/){
	print RED, "\nMaster probe file: $master_probe_file does not appear to be a valid probe file with a Probe_Count column!\n\n", RESET;
	exit();
      }

      #Print out the master_header_line and appended columns
      my @header = split ("\t", $header_line);
      my $new_header = '';
      shift(@header); #Remove the redundant 'Probe_Count' column
      foreach my $header (@header){
	$new_header = "$new_header"."\t$header";
      }
      print NEW_PROBE_FILE "$master_header_line$new_header\n";
      next();
    }

    my $master_probe_id = $line[0];

    unless ($master_probe_id =~ /^\d+/){
      print RED, "\nProbe_ID: $master_probe_id from master_probe_file does not appear to be valid!\n\n", RESET;
      exit();
    }

    #Make sure this probe ID was found in the files to be merged
    unless ($probes{$master_probe_id}){
      print RED, "\nProbe_ID: $master_probe_id from the master_probe_file was not found in any of the probe files to be merged!\n\n", RESET;
      exit();
    }
    $master_probes_found++;

    #Print out the current line from the master probe file with the additional probe data found for the current probe_ID
    my $data_line = $probes{$master_probe_id}{line_record};
    my @data_line = split ("\t", $data_line);
    my $new_data_line = '';
    shift(@data_line); #Remove the redundant 'Probe_Count' column
    foreach my $line_value (@data_line){
      $new_data_line = "$new_data_line"."\t$line_value";
    }
    print NEW_PROBE_FILE "$current_line_record$new_data_line\n";
  }
  close(MASTER_PROBE_FILE);

  #last check of probe counts
  my $new_probe_data_found = keys %probes;
  unless ($master_probes_found == $new_probe_data_found){
    print RED, "\nDid not find data for all probes in master_probe_file!!\n\n", RESET;
    exit();
  }

}else{
  #B.) Write out simple joined file

  #Write the header line
  print NEW_PROBE_FILE "$header_line\n";

  foreach my $probe_id (sort {$a <=> $b} keys %probes){
    print NEW_PROBE_FILE "$probes{$probe_id}{line_record}\n";
  }

}
close (NEW_PROBE_FILE);

exit();
