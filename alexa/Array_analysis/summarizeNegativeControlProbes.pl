#!/usr/local/bin/perl58 -w
#Written by Malachi Griffith
#This script takes a NimbleGen datafile containing probe intensity values for several hybridizations and combines this data
#with additional probe info creating a new summary file with probe intensities, Tm, probe type, exon skips, etc.
#This script does this only for negative control probes!!
#Creates a seperate line for each negative control probe, instead of placing 'junction and reverse-junction' probes on the same line

#To do this will involve the following steps
#1.) Import the filtered (a) negative control probes.
#2.) Import the NimbleGen probe intensity data file
#3.) Combine data from (1) and (2) 
#4.) Generate an output file, tab-delimited, with one negative control probe per line

use strict;
use Data::Dumper;
use Getopt::Long;

#Initialize command line options
my $negative_probes = '';
my $nimblegen_data_file = '';
my $out_file = '';

GetOptions ('negative_probes=s'=>\$negative_probes,'nimblegen_data=s'=>\$nimblegen_data_file, 'out_file=s'=>\$out_file);

#Provide instruction to the user
print "\n\nUsage:";
print "\n\tSpecify the file containing the filtered negative probes using: --negative_probes";
print "\n\tSpecify the file containing NimbleGen probe intensity data using: --nimblegen_data";
print "\n\tNOTE: this script expects a 'All Pair' NimbleGen summary file which contains intensities for multiple experiments";
print "\n\tbut limited other info about the probe and its position on the array";
print "\n\tSpecify the output file for submission to NimbleGen using: --out_file";
print "\n\nExample: summarizeNegativeControlProbes.pl --negative_probes=filteredNegativeControlProbes.txt --nimblegen_data=All_pair.txt --out=dataSummary.txt\n\n";

#Make sure all options were specified
unless ($negative_probes && $nimblegen_data_file && $out_file){
  print "\nOptions missing!\n\n";
  exit();
}

#1.) Import the filtered junction and exon probes
my %probes;
&importFilteredProbes('-negative_probe_file'=>$negative_probes);
my $probe_count = keys %probes;
print "\nFound $probe_count probes in the input probe files";

#2.) Import the NimblGen probe intensity data file 
my %probe_data;
my @experiment_ids;
&importNimbleGenDataFile('-nimblegen_data_file'=>$nimblegen_data_file);
print "\n\nImported data for the following experiment ids: @experiment_ids\n\n";


#3.) Combine data from (1) and (2) and generate an output file, tab-delimited, with one probe per line
my %probe_summary;
&combineProbeData();

#4.) Generate an output file, tab-delimited, with one experimental probe per line
&printSummaryDataFile('-out_file'=>$out_file);

#print Dumper %probe_summary;

print "\n\n";
exit();


###################################################################################################################
#1.) Import the probes from the filtered exon and junction probes files.                                          #
###################################################################################################################
sub importFilteredProbes{
  my %args = @_;
  my $negative_file = $args{'-negative_probe_file'};

  #Process the negative probe file
  open (NEGATIVE, "$negative_file") || die "\nCould not open negative probe file: $negative_file\n\n";
  print "\nImporting negative probes from: $negative_file";

  my $header = 1;
  while (<NEGATIVE>){
    #Skip the header
    if ($header == 1){
      $header = 0;
      next();
    }
    my @line = split ("\t", $_);
    my $probe_id = $line[0];
    $probes{$probe_id}{gene_id} = $line[1];
    $probes{$probe_id}{sequence} = $line[2];
    $probes{$probe_id}{probe_type} = $line[3];
    $probes{$probe_id}{tm} = $line[12];
    $probes{$probe_id}{probe_pair_id} = $line[10];
    $probes{$probe_id}{exons_skipped} = "na";
  }
  close (NEGATIVE);

  my $probe_records = keys %probes;
  print "\nFound $probe_records probe records from the input junction,exon and control probe files\n\n";

  return();
}


###########################################################################################################################
#2.) Import the NimbleGen probe intensity data file                                                                       #
###########################################################################################################################
sub importNimbleGenDataFile{
  my %args = @_;
  my $data_file = $args{'-nimblegen_data_file'};

  #Process the NimbleGen probe intensity data file
  open (DATA, "$data_file") || die "\nCould not open NimbleGen probe intensity data file: $data_file\n\n";
  print "\n\nImporting probes intensities from: $data_file\n\n";

  my $header = 1;
  my @columns;
  my $column_count;

  while (<DATA>){
    #Skip the header
    if ($header == 1){
      $header = 0;
      #Get the names of the columns so that the hybridization/experiment ID can be used
      chomp ($_);
      @columns = split ("\t", $_);
      $column_count = @columns;  
      @experiment_ids = @columns[4 .. ($column_count-1)];
      next();
    }
    my @line = split ("\t", $_);
    my $gene_id = $line[1];
    my $probe_id_uf = $line[2];
    my $probe_id;

    #Format the probe ID, skip lines corresponding to NimbleGen control probes
    if ($probe_id_uf =~ /^ALEXA_(\d+)/){
    	$probe_id = $1;
    }else{
    	print "\rSkipping NimbleGen probe: $probe_id_uf";
    	next();
    }

    #Create a probe record for this line
    $probe_data{$probe_id}{target_id} = $gene_id;

    #Get the intensity values for each column - skip the first 4 columns of probe info
    my @intensities;
    for (my $i = 4; $i <= $column_count-1; $i++){
    	my $intensity = $line[$i];
    	chomp ($intensity);
	
	push(@intensities, $intensity);
    }
    $probe_data{$probe_id}{intensities} = \@intensities;
  }
  close (DATA);

  my $probe_count = keys %probe_data;
  print "\nFound $probe_count probe data lines in the NimbleGen input file (not including nimblegen controls)\n\n";

  return();
}


##########################################################################################################################
#3.) Combine data from (1) and (2) and generate an output file, tab-delimited, with one probe per line                   #
##########################################################################################################################
sub combineProbeData{
  my $probe_count = 0;
  my $record_count = 0;

  #Go through all the NimbleGen probe intensity data and create a list of probes
  foreach my $probe_id (keys %probe_data){
    #Make sure this probe id was found in the probe info files
    unless ($probes{$probe_id}){
      next();
    }

    #Get the probe type for this probe
    my $probe_type = $probes{$probe_id}{probe_type};

    unless ($probe_type eq "Control-Negative" || $probe_type eq "Control-Negative_Rev"){
      next();
    }

    $probe_summary{$probe_id}{gene_id} = $probes{$probe_id}{gene_id};
    $probe_summary{$probe_id}{target_id} = $probe_data{$probe_id}{target_id};
    $probe_summary{$probe_id}{probe_type} = $probe_type;
    $probe_summary{$probe_id}{sequence} = $probes{$probe_id}{sequence};
    $probe_summary{$probe_id}{tm} = $probes{$probe_id}{tm};
    $probe_summary{$probe_id}{exons_skipped} = $probes{$probe_id}{exons_skipped};
    $probe_summary{$probe_id}{intensities_experimental} = $probe_data{$probe_id}{intensities};

    $probe_count++;
    $record_count++;
  }
  print "\n\nProbe pair summary created to account for $probe_count probes total";
  print "\n\nThese probes correspond to a total of $record_count records\n\n";

  return();
}


###########################################################################################################################
#3.) Combine data from (1) and (2) and generate an output file, tab-delimited, with one probe per line                    #
###########################################################################################################################
sub printSummaryDataFile{
  my %args = @_;
  my $data_file = $args{'-out_file'};

  #Open the output file
  open (OUTFILE, ">$out_file") || die "\nCould not open output file: $out_file\n\n";

  #Print the header line
  print OUTFILE "Probe_ID\tProbe_type\tSequence\tTm";
  foreach my $experiment_id (@experiment_ids){
    print OUTFILE "\t$experiment_id";
  }
  print OUTFILE "\n";  

  foreach my $probe_id (sort {$a <=> $b} keys %probe_summary){

    print OUTFILE "$probe_id\t$probe_summary{$probe_id}{probe_type}\t$probe_summary{$probe_id}{sequence}\t$probe_summary{$probe_id}{tm}";

    my @experimental_intensities = @{$probe_summary{$probe_id}{intensities_experimental}};
    foreach my $experimental_intensity (@experimental_intensities){
      print OUTFILE "\t$experimental_intensity";
    }
    print OUTFILE "\n";
  }

  close OUTFILE;

  return();
}
