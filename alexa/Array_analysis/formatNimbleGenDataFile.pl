#!/usr/bin/perl -w
#Written by Malachi Griffith
#This script takes a NimbleGen datafile containing probe intensity values for several hybridizations and combines this data
#with additional probe info creating a new summary file with probe intensities, Tm, probe type, exon skips, etc.

#To do this will involve the following steps
#1.) Import the NimbleGen probe intensity data file
#2.) Import the filtered (a) exon-junction probes, (b) exon-intron junction probes (c) exon probes, (d) intron probes, and (e) negative control probes.
#    - Get all these probes from a single probe dir
#    - Only do this for probes that are actually in the design file
#3.) Combine data from (1) and (2)
#4.) Generate an output file, tab-delimited, with one experimental probe per line
#    - Info from probe files is simply appended to the end of the file

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);

#Initialize command line options
my $nimblegen_data_file = '';
my $probe_id_column = '';
my $filtered_probes_dir = '';
my $exons_skipped_column = '';
my $out_file = '';

GetOptions ('nimblegen_data_file=s'=>\$nimblegen_data_file, 'probe_id_column=i'=>\$probe_id_column, 'filtered_probes_dir=s'=>\$filtered_probes_dir, 'exons_skipped_column=i'=>\$exons_skipped_column, 'out_file=s'=>\$out_file);

#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the file containing NimbleGen probe intensity data using: --nimblegen_data_file", RESET;
print GREEN, "\n\tSpecify the column in this file containing probe IDs using: --probe_id_column", RESET;
print GREEN, "\n\tNOTE: this script expects a 'All Pair' NimbleGen summary file which contains intensities for multiple experiments", RESET;
print GREEN, "\n\tbut limited other info about the probe and its position on the array", RESET;
print GREEN, "\n\tSpecify the directory containing filtered probe files using: --filtered_probes_dir", RESET;
print GREEN, "\n\tSpecify the column containing exons_skipped values for exon-exon probes using: --exons_skipped_column", RESET;
print GREEN, "\n\tSpecify the output file for submission to NimbleGen using: --out_file", RESET;
print "\n\nExample: formatNimbleGenDataFile.pl --nimblegen_data_file=/home/malachig/AlternativeSplicing/Array_analysis/NimbleGen_Data/NimbleGenData_04Dec2006/AuxillaryData/PairData_formatted/All_hybes.txt  --probe_id_column=2  --filtered_probes_dir=/home/malachig/AlternativeSplicing/Array_design/ALEXA_hs_35_35h/probe_dev/filtered_probes/  --exons_skipped_column=34  --out_file=dataSummary.txt\n\n", RESET;

#Make sure all options were specified
unless ($nimblegen_data_file && $probe_id_column && $filtered_probes_dir && $exons_skipped_column && $out_file){
  print RED, "\nOptions missing!\n\n", RESET;
  exit();
}

#1.) Import the NimblGen probe intensity data file
my %probes;
my $data_header;
&importNimbleGenDataFile('-nimblegen_data_file'=>$nimblegen_data_file, '-probe_id_column'=>$probe_id_column);
print BLUE, "\n\nImported data for the following data columns: $data_header\n\n", RESET;

#2.) Import the filtered junction and exon probes
&importFilteredProbes('-filtered_probes_dir'=>$filtered_probes_dir, '-exons_skipped_column'=>$exons_skipped_column);

#3.) Generate an output file, tab-delimited, with one experimental probe per line
&printSummaryDataFile('-out_file'=>$out_file);

print BLUE, "\n\n", RESET;

exit();


###########################################################################################################################
#1.) Import the NimbleGen probe intensity data file                                                                       #
###########################################################################################################################
sub importNimbleGenDataFile{
  my %args = @_;
  my $data_file = $args{'-nimblegen_data_file'};
  my $probe_id_column = $args{'-probe_id_column'};

  #Process the NimbleGen probe intensity data file
  open (DATA, "$data_file") || die "\nCould not open NimbleGen probe intensity data file: $data_file\n\n";
  print BLUE, "\n\nImporting probes intensities from: $data_file\n\n", RESET;

  my $header = 1;
  my @columns;
  my $column_count;
  my $line_count = 0;

  while (<DATA>){
    $line_count++;
    chomp ($_);

    #Skip the header
    if ($header == 1){
      $header = 0;
      #Get the names of the columns so that the hybridization/experiment ID can be used
      $data_header = $_;
      next();
    }
    my @line = split ("\t", $_);
    my $line = $_;
    my $probe_id_uf = $line[$probe_id_column-1];
    my $probe_id;

    #Format the probe ID, skip lines corresponding to NimbleGen control probes
    if ($probe_id_uf =~ /^ALEXA_(\d+)/){
      $probe_id = $1;
    }else{
      print YELLOW, "\rSkipping NimbleGen probe: $probe_id_uf", RESET;
      next();
    }

    #Create a probe record for this line
    $probes{$probe_id}{line} = $line;
    $probes{$probe_id}{line_count} = $line_count;

  }
  close (DATA);

  my $probe_count = keys %probes;
  print BLUE, "\nFound $probe_count probe data lines in the NimbleGen input file (not including nimblegen controls)\n\n", RESET;

  return();
}


###################################################################################################################
#2.) Import the probes from the filtered exon and junction probes files.                                          #
###################################################################################################################
sub importFilteredProbes{
  my %args = @_;
  my $filtered_probe_dir = $args{'-filtered_probes_dir'};
  my $exons_skipped_column = $args{'-exons_skipped_column'};

  #Information we want to append to the data file: alexa_gene_id, probeset_id, probe_type, probe_length, probe_sequence, probe_tm, exons_skipped
  print BLUE, "\nProcessing input probe files from $filtered_probe_dir to get probe counts for each gene\n\n", RESET;

  my @files = `ls $filtered_probe_dir`;

  print BLUE, "\nFound the following filtered probe files:\n @files", RESET;

  my $data_found = 0;
  my %genes_list;

  foreach my $probe_file (@files){

    chomp($probe_file);
    my $probe_file_path = "$filtered_probe_dir"."/"."$probe_file";

    open (PROBE, "$probe_file_path") || die "\nCould not open probe file: $probe_file_path\n\n";

    print BLUE, "\n\tCounting probes from: $probe_file", RESET;
    my $header = 1;

    while (<PROBE>){
      #Skip the header
      if ($header == 1){
	$header = 0;
	next();
      }
      chomp($_);
      my @line = split ("\t", $_);

      my $probe_id = $line[0];
      my $gene_id = $line[2];

      #Skip this probe unless it was found in the data file
      unless ($probes{$probe_id}){
	next();
      }
      $data_found++;

      $probes{$probe_id}{probeset_id} = $line[1];
      $probes{$probe_id}{gene_id} =  $line[2];
      $probes{$probe_id}{sequence} = $line[3];
      $probes{$probe_id}{length} = $line[4];
      $probes{$probe_id}{tm} = $line[5];
      $probes{$probe_id}{probe_type} = $line[6];

      #Count genes found
      $genes_list{$line[2]}{tmp} = '';

      #For exon-exon probes, get the number of exons skipped
      if ($probes{$probe_id}{probe_type} eq "Exon-Exon"){
	$probes{$probe_id}{exons_skipped} = $line[$exons_skipped_column-1];
      }else{
	$probes{$probe_id}{exons_skipped} = 'na';
      }
    }
    close (PROBE);
  }

  my $genes_with_probes = keys %genes_list;

  print BLUE, "\nFound a probe info for a total of $data_found probes corresponding to $genes_with_probes gene targets\n\n", RESET;

  return();
}


###########################################################################################################################
#3.) Generate an appended output file, tab-delimited, with one probe per line                                             #
###########################################################################################################################
sub printSummaryDataFile{
  my %args = @_;
  my $out_file = $args{'-out_file'};

  print BLUE, "\nPrinting appended output data file: $out_file\n\n", RESET;

  my $records_printed = 0;

  #Open the output file
  open (OUTFILE, ">$out_file") || die "\nCould not open output file: $out_file\n\n";

  #Print the header line
  print OUTFILE "AlexaGene_ID\tProbe_ID\tProbeSet_ID\t$data_header\tSequence\tProbe_length\tProbe_Tm\tProbe_Type\tExons_Skipped\n";

  foreach my $probe_id (sort {$probes{$a}->{line_count} <=> $probes{$b}->{line_count}} keys %probes){
    $records_printed++;

    print OUTFILE "$probes{$probe_id}{gene_id}\t$probe_id\t$probes{$probe_id}{probeset_id}\t$probes{$probe_id}{line}\t$probes{$probe_id}{sequence}\t$probes{$probe_id}{length}\t$probes{$probe_id}{tm}\t$probes{$probe_id}{probe_type}\t$probes{$probe_id}{exons_skipped}\n";

  }

  close OUTFILE;

  print BLUE, "\nPrinted $records_printed probe records to the output file\n\n", RESET;

  return();
}
