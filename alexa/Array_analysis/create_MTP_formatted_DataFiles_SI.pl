#!/usr/bin/perl

#Written by Malachi Griffith
#This script takes a data file containing gene-normalized values as input
#NOTE: these values are already on a log2 scale!
#This script is used in conjuction with the output of processIntensityData_SI.pl, but it had to be modified slightly to be in the right format


#It then creates output files to be used with the R MTP package to test for differences between two conditions
# e.g. compare observations for sensitive cell lines versus resistant cell lines
#The Input file may contain normalized and background corrected values instead of raw values

#Values will be grouped to allow statistical testing of differences between conditions at two levels:
#   A.) The exon level, where each probeset represents a single exon (or portion of an exon)
#   B.) The gene level, where all exon and canonical junction probesets for a single gene are grouped

#1.) Import probe data file
#    - Note: first three columns are expected to contain: AlexaGene_ID\tProbe_ID\tProbeSet_ID
#2.) Create a hash of genes -> probesets -> probes
#3.) Calculate log2 values for all intensities
#4.) Arrange probe-level values and output to exon-level_MTP.txt
#5.) Arrange gene-level values and output to gene-level_MTP.txt

#Note: (sens1, ..., sens3 are three replicate experiments for the sensitive cell line)
#Note: (res1, ..., res3 are three replicate experiments for the resistant cell line)
#Note: (probe1, ..., probe3 are three probes for a single exon)

#Desired format:
#exon1, probe1_sens1, probe1_sens2, probe1_sens3, ... , probe3_sens3, ..., probe1_res1, probe1_res2, probe1_res3, ... , probe3_res3

#Each line (row) will contain all the observations for a particular exon (probeset) or gene (gene_id)
#This script can also be run with files using Affy/NimbleGen common probeset IDs

#Note that before writing the output file, we need to know the max number of observations for a single exon/gene in each condition
#- This will define the number of columns.
#- Since some probesets are smaller than others, there will have to be 'NA' values for some exons/genes in each condition

#Note: Process negative control probes for 'exon'-level analysis but not for gene-level (results in simply too many columns)


use strict;
use Data::Dumper;
use Getopt::Long;
use Math::Complex;
use Term::ANSIColor qw(:constants);

use lib '/home/malachig/AlternativeSplicing/perl_bin';
use utilities::utility qw(:all);
use utilities::ALEXA_DB qw(:all);
use Statistics::Descriptive;

#Initialize command line options
my $data_file = '';
my $out_dir = '';
my $file_type = '';

GetOptions ('data_file=s'=>\$data_file, 'out_dir=s'=>\$out_dir, 'file_type=s'=>\$file_type);

#Provide instruction to the user
print GREEN, "\n\nUsage:", RESET;
print GREEN, "\n\tSpecify the input probe microarray data file using: --data_file", RESET;
print GREEN, "\n\tSpecify a directory for output files using: --out_dir", RESET;
print GREEN, "\n\tSpecify a file type using: --file_type (AA_affy, AA_alexa, or Standard)", RESET;
print GREEN, "\n\nExample: create_MTP_formatted_DataFiles_SI.pl  --data_file=/home/malachig/AlternativeSplicing/Array_analysis/NimbleGen_Data/MIPvs5FUR_04Dec2006/AuxillaryData/PairData_formatted/Quantiles_TmBGC/All_hybes_withProbeInfo_BGC_Norm_GNI.txt  --out_dir=/home/malachig/AlternativeSplicing/Array_analysis/NimbleGen_Data/MIPvs5FUR_04Dec2006/AuxillaryData/Processed_data/Quantiles_TmBGC/MTP_formatted  --file_type=Standard\n\n", RESET;

#Make sure all options were specified
unless ($data_file && $out_dir && $file_type){
  print RED, "\nOptions missing!\n\n", RESET;
  exit();
}
unless ($file_type eq "AA_affy" || $file_type eq "AA_alexa" || $file_type eq "Standard"){
  print RED, "\n\t--file_type must be 'AA_affy', 'AA_alexa', or 'Standard'\n\n",RESET;
  exit();
}

#Define sample pairs that will be used to calculate DE and groups of replicates that will be used to calculate means
#Only data columns defined here will be imported
my %samples;
if ($file_type eq "AA_affy"){
  $samples{MIP_C}{group} = "mip";
  $samples{MIP_C}{pair} = "c";
  $samples{MIP_C}{order} = 1;
}else{
  $samples{MIP_EF_A}{group} = "mip";
  $samples{MIP_EF_A}{pair} = "ef_a";
  $samples{MIP_EF_A}{order} = 1;
}
$samples{MIP_EF_B}{group} = "mip";
$samples{MIP_EF_B}{pair} = "ef_b";
$samples{MIP_EF_B}{order} = 2;
$samples{MIP_GH_A}{group} = "mip";
$samples{MIP_GH_A}{pair} = "gh_a";
$samples{MIP_GH_A}{order} = 3;

if ($file_type eq "AA_affy"){
  $samples{FUR_C}{group} = "fur";
  $samples{FUR_C}{pair} = "c";
  $samples{FUR_C}{order} = 4;
}else{
  $samples{FUR_EF_A}{group} = "fur";
  $samples{FUR_EF_A}{pair} = "ef_a";
  $samples{FUR_EF_A}{order} = 4;
}
$samples{FUR_EF_B}{group} = "fur";
$samples{FUR_EF_B}{pair} = "ef_b";
$samples{FUR_EF_B}{order} = 5;
$samples{FUR_GH_A}{group} = "fur";
$samples{FUR_GH_A}{pair} = "gh_a";
$samples{FUR_GH_A}{order} = 6;

my %pairs;
if ($file_type eq "AA_affy"){
  my @arr1 = qw (MIP_C FUR_C); $pairs{c}{pair} = \@arr1; $pairs{c}{order} = 1;
}else{
  my @arr1 = qw (MIP_EF_A FUR_EF_A); $pairs{ef_a}{pair} = \@arr1; $pairs{ef_a}{order} = 1;
}
my @arr2 = qw (MIP_EF_B FUR_EF_B); $pairs{ef_b}{pair} = \@arr2; $pairs{ef_b}{order} = 2;
my @arr3 = qw (MIP_GH_A FUR_GH_A); $pairs{gh_a}{pair} = \@arr3; $pairs{gh_a}{order} = 3;

my %replicates;
if ($file_type eq "AA_affy"){
  my @rep1 = qw (MIP_C MIP_EF_B MIP_GH_A); $replicates{mip}{replicate} = \@rep1; $replicates{mip}{order} = 1;
  my @rep2 = qw (FUR_C FUR_EF_B FUR_GH_A); $replicates{fur}{replicate} = \@rep2; $replicates{fur}{order} = 2;
}else{
  my @rep1 = qw (MIP_EF_A MIP_EF_B MIP_GH_A); $replicates{mip}{replicate} = \@rep1; $replicates{mip}{order} = 1;
  my @rep2 = qw (FUR_EF_A FUR_EF_B FUR_GH_A); $replicates{fur}{replicate} = \@rep2; $replicates{fur}{order} = 2;
}

#1.) Import probe data file
#    - Note: first three columns are expected to contain: AlexaGene_ID\tProbe_ID\tProbeSet_ID
#2.) Create a hash of genes -> probesets -> probes
my $largest_exon_probeset = 0;
my $largest_gene_probeset = 0;
my $data_ref = &parseProbeFile('-data_file'=>$data_file);

#Open output files
unless ($out_dir =~ /.*\/$/){
  $out_dir = "$out_dir"."/";
}
my $exon_log2_out = "$out_dir"."$file_type"."_exon_log2_data.txt";
my $exon_log2_de_out = "$out_dir"."$file_type"."_exon_log2_de_data.txt";
my $gene_log2_out = "$out_dir"."$file_type"."_gene_log2_data.txt";
my $gene_log2_de_out = "$out_dir"."$file_type"."_gene_log2_de_data.txt";

open (EXON_LOG2, ">$exon_log2_out") || die "\nCould not open output file: $exon_log2_out\n\n";
open (EXON_DE, ">$exon_log2_de_out") || die "\nCould not open output file: $exon_log2_de_out\n\n";
open (GENE_LOG2, ">$gene_log2_out") || die "\nCould not open output file: $gene_log2_out\n\n";
open (GENE_DE, ">$gene_log2_de_out") || die "\nCould not open output file: $gene_log2_de_out\n\n";

print BLUE, "\nPrinting exon-level and gene-level, log2 and log2 DE values to:\n\t$exon_log2_out\n\t$exon_log2_de_out\n\t$gene_log2_out\n\t$gene_log2_de_out\n\n", RESET;

#Print output file headers

#A.) Exon-level (probeset level) - log2 intensity measures
#    ProbeSet_ID, AlexaGene_ID, Probe_Count, Probe_Type, Exons_Skipped, [all log2 observations for each condition for this exon]
print EXON_LOG2 "ProbeSet_ID\tAlexaGene_ID\tProbe_Count\tProbe_Type\tExons_Skipped\t";
foreach my $rep (sort {$replicates{$a}->{order} <=> $replicates{$b}->{order}} keys %replicates){
  my @samples = @{$replicates{$rep}{replicate}};

  foreach my $sample (@samples){

    #Print out a title for every possible column
    for (my $i = 1; $i <= $replicates{$rep}{largest_exon_probeset}; $i++){
      print EXON_LOG2 "$sample","_$i\t";
    }
  }
}
print EXON_LOG2 "\n";


#B.) Exon-level (probeset level) - log2DE intensity measures
#    ProbeSet_ID, AlexaGene_ID, Probe_Count, Probe_Type, Exons_Skipped, [all log2 observations for each condition for this exon]
print EXON_DE "ProbeSet_ID\tAlexaGene_ID\tProbe_Count\tProbe_Type\tExons_Skipped\t";
foreach my $pair (sort {$pairs{$a}->{order} <=> $pairs{$b}->{order}} keys %pairs){

  #Print out a title for every possible column
  for (my $i = 1; $i <= $pairs{$pair}{largest_exon_probeset}; $i++){
    print EXON_DE "$pair","_$i\t";
  }
}
print EXON_DE "\n";


#C.) Gene-level - log2 intensity measures
#    AlexaGene_ID, Probe_Count, [all log2 observations for each condition for this gene]
print GENE_LOG2 "AlexaGene_ID\tProbe_Count\t";
foreach my $rep (sort {$replicates{$a}->{order} <=> $replicates{$b}->{order}} keys %replicates){
  my @samples = @{$replicates{$rep}{replicate}};

  foreach my $sample (@samples){

    #Print out a title for every possible column
    for (my $i = 1; $i <= $replicates{$rep}{largest_gene_probeset}; $i++){
      print GENE_LOG2 "$sample","_$i\t";
    }
  }
}
print GENE_LOG2 "\n";

#D.) Gene-level - log2DE intensity measures
#    AlexaGene_ID, Probe_Count, [all log2 observations for each condition for this gene]
print GENE_DE "AlexaGene_ID\tProbe_Count\t";
foreach my $pair (sort {$pairs{$a}->{order} <=> $pairs{$b}->{order}} keys %pairs){

  #Print out a title for every possible column
  for (my $i = 1; $i <= $pairs{$pair}{largest_gene_probeset}; $i++){
    print GENE_DE "$pair","_$i\t";
  }
}
print GENE_DE "\n";





#3.) Calculate probe-level values and output to probe-level.txt

#Go through each gene
foreach my $gene_id (sort {$a <=> $b} keys %{$data_ref}){

  print BLUE, "\n\tProcessing gene id: $gene_id", RESET;

  my $single_gene_ref = $data_ref->{$gene_id};

  #Each probeset
  my $probesets_ref = $single_gene_ref->{probesets};
  foreach my $probeset_id (sort {$a <=> $b} keys %{$probesets_ref}){

    my $probeset_values_ref = $probesets_ref->{$probeset_id}->{probeset_values};
    my $probeset_pair_values_ref = $probesets_ref->{$probeset_id}->{probeset_pair_values};

    my $probe_type = $probesets_ref->{$probeset_id}->{probe_type};
    my $exons_skipped = $probesets_ref->{$probeset_id}->{exons_skipped};

    #Each probe
    my $probes_ref = $probesets_ref->{$probeset_id}->{probes};
    foreach my $probe_id (sort {$a <=> $b} keys %{$probes_ref}){

      #Each probe/sample value
      my $values_ref = $probes_ref->{$probe_id}->{raw_values};

      #i.) Calculate log2 values for individual probes
      my %log2_values;
      foreach my $sample (sort keys %{$values_ref}){
	my $x = $values_ref->{$sample};
	my $rounded = sprintf("%.4f", $x);
	$log2_values{$sample} = $rounded;

	#Gather log2 values for each probe of this probeset for each sample
	push(@{$probeset_values_ref->{$sample}->{log2_values}}, $rounded);

	#Gather gene level values
	#If the probe type is correct (exon probe or canonical exon junction probe) add to the list of probe used to estimate gene expression
	#Also group all the negative control probes and come up with a single estimate for these (for interests sake)
	if ($probe_type eq "Exon" || (($probe_type eq "Exon-Exon") && ($exons_skipped eq "0")) || $probe_type eq "Control-Negative"){
	  push(@{$single_gene_ref->{gene_values}->{$sample}->{log2_values}}, $rounded);
	}
      }
      $probes_ref->{$probe_id}->{log2_values} = \%log2_values;

      #ii.) Calculate de values for individual probes
      my %log2_de_values;
      foreach my $pair (sort keys %pairs){
	my @pair = @{$pairs{$pair}{pair}};

	my $de = ($values_ref->{$pair[0]}) - ($values_ref->{$pair[1]});
	my $rounded = sprintf("%.4f", $de);
	$log2_de_values{$pair} = $rounded;

	#Gather log2 de values for each probe sample pair of this probeset
	push(@{$probeset_pair_values_ref->{$pair}->{log2_de_values}}, $rounded);

	#Gather gene level values
	#If the probe type is correct (exon probe or canonical exon junction probe) add to the list of probe used to estimate gene expression
	#Also group all the negative control probes and come up with a single estimate for these (for interests sake)
	if ($probe_type eq "Exon" || (($probe_type eq "Exon-Exon") && ($exons_skipped eq "0")) || $probe_type eq "Control-Negative"){
	  push(@{$single_gene_ref->{gene_pair_values}->{$pair}->{log2_de_values}}, $rounded);
	}

      }
      $probes_ref->{$probe_id}->{log2_de_values} = \%log2_de_values;

    }#probe loop

    #Now that values have been grouped to the probeset level - perform calculations on these

  }#probesets loop

  ####################################################################################################################################
  #Calculations complete for this gene.  Now print results to output files

  #A.) Exon-level (probeset level) - log2 intensity measures
  #    ProbeSet_ID, AlexaGene_ID, Probe_Count, Probe_Type, Exons_Skipped, [all log2 observations for each condition for this exon]
  foreach my $probeset_id (sort {$a <=> $b} keys %{$probesets_ref}){

    my $probe_type = $probesets_ref->{$probeset_id}->{probe_type};
    my $exons_skipped = $probesets_ref->{$probeset_id}->{exons_skipped};
    my $probe_count = $probesets_ref->{$probeset_id}->{probe_count};

    print EXON_LOG2 "$probeset_id\t$gene_id\t$probe_count\t$probe_type\t$exons_skipped\t";

    #For each condition (e.g. sensitive and resistant)
    foreach my $rep (sort {$replicates{$a}->{order} <=> $replicates{$b}->{order}} keys %replicates){
      #For each probe
      my $probes_ref = $probesets_ref->{$probeset_id}->{probes};

      #For each sample measured with this probe in this condition
      my @samples = @{$replicates{$rep}{replicate}};
      foreach my $sample (@samples){
	my $exon_values_printed = 0;

	foreach my $probe_id (sort {$a <=> $b} keys %{$probes_ref}){
	  my $log2_values_ref = $probes_ref->{$probe_id}->{log2_values};

	  print EXON_LOG2 "$log2_values_ref->{$sample}\t";

	  $exon_values_printed++;
	}

	#Now that all values for this sample-condition have been printed, fill in leftover spaces with 'NA'
	if ($exon_values_printed < $replicates{$rep}{largest_exon_probeset}){
	  for (my $i = $exon_values_printed+1; $i <= $replicates{$rep}{largest_exon_probeset}; $i++){
	    print EXON_LOG2 "NA\t";
	  }
	}
      }
    }
    print EXON_LOG2 "\n";
  }

  #B.) Exon-level (probeset level) - log2DE intensity measures
  #    ProbeSet_ID, AlexaGene_ID, Probe_Count, Probe_Type, Exons_Skipped, [all log2 observations for each condition for this exon]

  foreach my $probeset_id (sort {$a <=> $b} keys %{$probesets_ref}){

    my $probe_type = $probesets_ref->{$probeset_id}->{probe_type};
    my $exons_skipped = $probesets_ref->{$probeset_id}->{exons_skipped};
    my $probe_count = $probesets_ref->{$probeset_id}->{probe_count};

    print EXON_DE "$probeset_id\t$gene_id\t$probe_count\t$probe_type\t$exons_skipped\t";

    my $probeset_pair_values_ref = $probesets_ref->{$probeset_id}->{probeset_pair_values};

    #For each pair (e.g. sensitive vs resistant - times the number of replicates)
    foreach my $pair (sort {$pairs{$a}->{order} <=> $pairs{$b}->{order}} keys %pairs){
      my $exon_values_printed = 0;

      my @probe_de_values = @{$probeset_pair_values_ref->{$pair}->{log2_de_values}};

      foreach my $value (@probe_de_values){
	print EXON_DE "$value\t";
	 $exon_values_printed++;
      }

      #Now that all values for this condition pair have been printed, fill in leftover spaces with 'NA'
      if ($exon_values_printed < $pairs{$pair}{exon_column_count}){
	for (my $i = $exon_values_printed+1; $i <= $pairs{$pair}{largest_exon_probeset}; $i++){
	  print EXON_DE "NA\t";
	}
      }
    }
    print EXON_DE "\n";
  }

  #C.) Gene-level - log2 intensity measures
  #    AlexaGene_ID, Probe_Count, [all log2 observations for each condition for this gene]

  #Skip the negative control probes (which have gene_id == 0)
  unless ($gene_id == 0){
    my $gene_probe_count = $single_gene_ref->{gene_probe_count};
    print GENE_LOG2 "$gene_id\t$gene_probe_count\t";

    #For each condition (e.g. sensitive and resistant)
    foreach my $rep (sort {$replicates{$a}->{order} <=> $replicates{$b}->{order}} keys %replicates){

      my @samples = @{$replicates{$rep}{replicate}};

      #For each sample measured with this probe in this condition
      foreach my $sample (@samples){

	my @log2_values = @{$single_gene_ref->{gene_values}->{$sample}->{log2_values}};

	my $gene_log2_values_printed = 0;

	foreach my $value (@log2_values){
	  print GENE_LOG2 "$value\t";
	  $gene_log2_values_printed++;
	}

	#Now that all values for this sample condition have been printed, fill in leftover spaces with 'NA'
	if ($gene_log2_values_printed < $replicates{$rep}{largest_gene_probeset}){
	  for (my $i = $gene_log2_values_printed+1; $i <= $replicates{$rep}{largest_gene_probeset}; $i++){
	    print GENE_LOG2 "NA\t";
	  }
	}
      }

    }
    print GENE_LOG2 "\n";
  }

  #D.) Gene-level - log2DE intensity measures
  #    AlexaGene_ID, Probe_Count, [all log2 observations for each condition for this gene]
  #Skip the negative control probes (which have gene_id == 0)
  unless ($gene_id == 0){
    my $gene_probe_count = $single_gene_ref->{gene_probe_count};
    print GENE_DE "$gene_id\t$gene_probe_count\t";

    #For each condition (e.g. sensitive and resistant)
    foreach my $pair (sort {$pairs{$a}->{order} <=> $pairs{$b}->{order}} keys %pairs){

      my @de_values = @{$single_gene_ref->{gene_pair_values}->{$pair}->{log2_de_values}};
      my $gene_de_values_printed = 0;

      #For each sample measured with this probe in this condition
      foreach my $value (@de_values){
	print GENE_DE "$value\t";
	$gene_de_values_printed++;
      }

      #Now that all values for this sample condition have been printed, fill in leftover spaces with 'NA'
      if ($gene_de_values_printed < $pairs{$pair}{largest_gene_probeset}){
	for (my $i = $gene_de_values_printed+1; $i <= $pairs{$pair}{largest_gene_probeset}; $i++){
	  print GENE_DE "NA\t";
	}
      }
    }
    print GENE_DE "\n";
  }



  #print Dumper $single_gene_ref;

}#gene loop

print BLUE "\n\n", RESET;

close (EXON_LOG2);
close (EXON_DE);
close (GENE_LOG2);
close (GENE_DE);

exit();


#################################################################################################################
#Parse data file
#################################################################################################################
sub parseProbeFile{
  my %args = @_;
  my $data_file = $args{'-data_file'};
  my %data;

  open (DATA, "$data_file") || die "\nCould not open data file: $data_file\n\n";

  my $first_line = 1;
  my %columns;

  while(<DATA>){
    chomp($_);
    my @line = split("\t", $_);

    #Watch for the header line
    if ($first_line == 1){
      $first_line = 0;
      my $column_count = 0;
      foreach my $column (@line){
	$columns{$column}{column_pos} = $column_count;
	$column_count++;
      }

      #Check for required columns
      unless($columns{'AlexaGene_ID'} && $columns{'Probe_ID'} && $columns{'ProbeSet_ID'} && $columns{'Probe_Type'} && $columns{'Exons_Skipped'}){
	print RED, "\nExpected column missing\n\n", RESET;
	exit();
      }

      if ($file_type eq "AA_affy" || $file_type eq "AA_alexa"){
	unless ($columns{'AA_probeset_id'}){
	  print RED, "\nAA_probeset_id column missing\n\n", RESET;
	  exit();
	}
      }

      next();
    }

    my $gene_id = $line[$columns{'AlexaGene_ID'}{column_pos}];

    #################################################################################
    #Skip gene: H19 which has unusually large probsets??
    if ($gene_id == 6617){
      print YELLOW, "\nSkipping probe for gene: H19", RESET;
      next();
    }

    #Process H19 seperately
    #unless ($gene_id == 6617){
    #  next();
    #}


    #Change gene_id for NC probes from 'na' to '0'
    if ($gene_id eq "na") {$gene_id = 0;}

    #Get the probeset ID (if this is a common probeset datafile, replace the normal probeset ID with the AA_probeset_ID
    my $probeset_id;

    if ($file_type eq "AA_affy" || $file_type eq "AA_alexa"){
      $probeset_id = $line[$columns{'AA_probeset_id'}{column_pos}];
    }else{
      $probeset_id = $line[$columns{'ProbeSet_ID'}{column_pos}];
    }

    my $probe_id = $line[$columns{'Probe_ID'}{column_pos}];
    my $probe_type = $line[$columns{'Probe_Type'}{column_pos}];
    my $exons_skipped = $line[$columns{'Exons_Skipped'}{column_pos}];

    if ($data{$gene_id}){
      #If this gene was previously observed 
      my $probesets_ref = $data{$gene_id}{probesets};

      if ($probe_type eq "Exon" || (($probe_type eq "Exon-Exon") && ($exons_skipped eq "0"))){
	$data{$gene_id}{gene_probe_count}++;  #Number of probes that will be used to estimate gene expression
      }

      if ($probesets_ref->{$probeset_id}){
	#If this probeset was previously observed
	my $probes_ref = $probesets_ref->{$probeset_id}->{probes};

	my %raw_values;
	foreach my $sample (sort keys %samples){
	  if ($columns{$sample}){
	    $raw_values{$sample} = $line[$columns{$sample}{column_pos}];
	  }
	}
	$probes_ref->{$probe_id}->{raw_values} = \%raw_values;
	$probesets_ref->{$probeset_id}->{probe_count}++;

      }else{
	#First time this probeset has been observed
	my %probes;

	my %raw_values;
	my %probeset_values;

	foreach my $sample (sort keys %samples){
	  if ($columns{$sample}){
	    $raw_values{$sample} = $line[$columns{$sample}{column_pos}];

	    my @tmp1;
	    $probeset_values{$sample}{log2_values} = \@tmp1;
	  }
	}
	$probes{$probe_id}{raw_values} = \%raw_values;

	$probesets_ref->{$probeset_id}->{probeset_values} = \%probeset_values;
	$probesets_ref->{$probeset_id}->{probes} = \%probes;
	$probesets_ref->{$probeset_id}->{probe_type} = $probe_type;
	$probesets_ref->{$probeset_id}->{exons_skipped} = $exons_skipped;
	$probesets_ref->{$probeset_id}->{probe_count} = 1;

	#Initialize arrays to store de values for all samples pairs for the entire probeset
	my %probeset_pair_values;
	foreach my $pair (sort keys %pairs){
	  my @tmp1;
	  $probeset_pair_values{$pair}{log2_de_values} = \@tmp1;
	}
	$probesets_ref->{$probeset_id}->{probeset_pair_values} = \%probeset_pair_values;

      }
    }else{
      #First time gene was observed
      my %probes;
      my %probesets;

      my %raw_values;
      my %probeset_values;
      my %gene_values;

      foreach my $sample (sort keys %samples){
	if ($columns{$sample}){
	  $raw_values{$sample} = $line[$columns{$sample}{column_pos}];

	  my @tmp1;
	  my @tmp2;
	  $probeset_values{$sample}{log2_values} = \@tmp1;
	  $gene_values{$sample}{log2_values} = \@tmp2;

	}
      }
      $probes{$probe_id}{raw_values} = \%raw_values;

      $probesets{$probeset_id}{probeset_values} = \%probeset_values;
      $probesets{$probeset_id}{probes} = \%probes;
      $probesets{$probeset_id}{probe_type} = $probe_type;
      $probesets{$probeset_id}{exons_skipped} = $exons_skipped;
      $probesets{$probeset_id}{probe_count} = 1;

      #Initialize arrays to store de values for all samples pairs for the entire probeset
      my %probeset_pair_values;
      my %gene_pair_values;
      foreach my $pair (sort keys %pairs){
	my @tmp1;
	my @tmp2;
	$probeset_pair_values{$pair}{log2_de_values} = \@tmp1;
	$gene_pair_values{$pair}{log2_de_values} = \@tmp2;

      }
      $probesets{$probeset_id}{probeset_pair_values} = \%probeset_pair_values;

      $data{$gene_id}{probesets} = \%probesets;
      $data{$gene_id}{gene_values} = \%gene_values;
      $data{$gene_id}{gene_pair_values} = \%gene_pair_values;

      if ($probe_type eq "Exon" || (($probe_type eq "Exon-Exon") && ($exons_skipped eq "0"))){
	$data{$gene_id}{gene_probe_count} = 1;  #Number of probes that will be used to estimate gene expression
      }else{
	$data{$gene_id}{gene_probe_count} = 0;
      }
    }
  }

  close (DATA);

  #Determine the largest number of observations for a single condition (number of probes in a probeset or for a gene)
  #This will be used to determine the number of columns required for the exon-level and gene-level file
  foreach my $gene_id (keys %data){

    my $probesets_ref = $data{$gene_id}{probesets};

    foreach my $probeset_id (keys %{$probesets_ref}){
      if ($probesets_ref->{$probeset_id}->{probe_count} > $largest_exon_probeset){
	$largest_exon_probeset = $probesets_ref->{$probeset_id}->{probe_count};
      }

      if ($data{$gene_id}{gene_probe_count} > $largest_gene_probeset){
	$largest_gene_probeset = $data{$gene_id}{gene_probe_count};
      }
    }
  }

  my $rep_count;
  foreach my $rep (keys %replicates){
    my @reps = @{$replicates{$rep}{replicate}};
    $rep_count = @reps;

    $replicates{$rep}{largest_exon_probeset} = $largest_exon_probeset;
    $replicates{$rep}{largest_gene_probeset} = $largest_gene_probeset;
    $replicates{$rep}{exon_column_count} = $rep_count * $largest_exon_probeset;
    $replicates{$rep}{gene_column_count} = $rep_count * $largest_gene_probeset;

    print BLUE, "\n\nCondition: $rep will have $replicates{$rep}{exon_column_count} columns per exon and $replicates{$rep}{gene_column_count} columns per gene\n\n", RESET;
  }

  foreach my $pair (keys %pairs){
    $pairs{$pair}{largest_exon_probeset} = $largest_exon_probeset;
    $pairs{$pair}{largest_gene_probeset} = $largest_gene_probeset;
    $pairs{$pair}{exon_column_count} = $rep_count * $largest_exon_probeset;
    $pairs{$pair}{gene_column_count} = $rep_count * $largest_gene_probeset;
  }

  return(\%data);
}












