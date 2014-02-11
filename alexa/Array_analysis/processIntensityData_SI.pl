#!/usr/bin/perl

#Written by Malachi Griffith
#This script takes a raw data file as input and calculates mean, log2, DE, etc. values
#Input file may contain normalized and background corrected values instead of raw values
#Values will be summarized at three levels and outputed to three different files

#1.) Import probe data file
#    - Note: first three columns are expected to contain: AlexaGene_ID\tProbe_ID\tProbeSet_ID
#2.) Create a hash of genes -> probesets -> probes
#3.) Calculate probe-level values
#4.) Calculate probeset-level values
#5.) Calculate gene-level values

#6.) Calculate GENE-NORMALIZED values.
#    - Individual probe values divided by gene expression estimate
#    - Mean probeset values divided by gene expression estimate

#7.) Calculate SIMPLE SI values for all probes
#    - Try the formula: log2((exon1_mip/gene1_mip) / (exon1_fur/gene1_fur))

#    - Note: log2(exon1_mip/gene1_mip) - log2(exon1_fur/gene1_fur)
#            IS THE SAME AS
#            (log2(exon1_mip)-log2(gene1_mip)) - (log2(exon1_fur) - log2(gene1_fur))

#8.) Calculate COMPLEX SI values for applicable probe types
#    - Calculate an SI value based on inclusion/exclusion evidence
#    - For each exon skipping probe, get the probes that provide evidence for the inclusion of an exon (or series of exons) and compare to
#      probes that provide evidence for the exons being skipped
#    - Exclusion probes: the exon skipping probe and other members of the same probeset (covering the same junction)
#    - Inclusion probes: exon and canonical junction probes within the affected region
#       - Exon probes that are entirely within the unit1_end and unit2_start of the current exon-skip probe
#       - Canonical junction probes that are entirely within the unit1_end and unit2_start of the current exon-skip probe
#       - OR canonical junction probes that have the same unit1_end or unit2_start as the current exon-skip probe

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
print GREEN, "\n\nExample: processIntensityData.pl  --data_file=/home/malachig/AlternativeSplicing/Array_analysis/NimbleGen_Data/MIPvs5FUR_04Dec2006/AuxillaryData/PairData_formatted/Quantiles_TmBGC/All_hybes_withProbeInfo_BGC_Norm.txt  --out_dir=/home/malachig/AlternativeSplicing/Array_analysis/NimbleGen_Data/MIPvs5FUR_04Dec2006/AuxillaryData/Processed_data/Quantiles_TmBGC/  --file_type=Standard\n\n", RESET;

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
my $data_ref = &parseProbeFile('-data_file'=>$data_file);

#Open output files
unless ($out_dir =~ /.*\/$/){
  $out_dir = "$out_dir"."/";
}
my $probe_out = "$out_dir"."$file_type"."_"."probe_data.txt";
my $probeset_out = "$out_dir"."$file_type"."_"."probeset_data.txt";
my $gene_out = "$out_dir"."$file_type"."_"."gene_data.txt";
open (PROBE, ">$probe_out") || die "\nCould not open output file: $probe_out\n\n";
open (PROBESET, ">$probeset_out") || die "\nCould not open output file: $probeset_out\n\n";

print BLUE, "\nPrinting probe-level and probeset-level to:\n\t$probe_out\n\t$probeset_out\n\n", RESET;

#Print output file headers
#Probe-level
#Probe_ID, ProbeSet_ID, AlexaGene_ID, Probe_Type, Exons_Skipped, log2_gni (all reps), log2_si (all_pairs), mean_log2_gni (mean of reps), mean_log2_si (mean of pairs)
print PROBE "Probe_ID\tProbeSet_ID\tAlexaGene_ID\tProbe_Type\tExons_Skipped\t";
foreach my $sample (sort {$samples{$a}->{order} <=> $samples{$b}->{order}} keys %samples){
  print PROBE "$sample","_LOG2_GNI\t";
}
foreach my $pair (sort {$pairs{$a}->{order} <=> $pairs{$b}->{order}} keys %pairs){
  print PROBE "$pair","_LOG2_SI\t";
}
foreach my $rep (sort {$replicates{$a}->{order} <=> $replicates{$b}->{order}} keys %replicates){
  print PROBE "$rep","_LOG2_GNI_U\t";
}
print PROBE "mip_v_fur_SI_U\n";

#Probeset-level
#ProbeSet_ID, AlexaGene_ID, Probe_Count, Probe_Type, Exons_Skipped, log2_gni (all reps), log2_si (all_pairs), mean_log2_gni (mean of reps), mean_log2_si (mean of pairs)
print PROBESET "ProbeSet_ID\tAlexaGene_ID\tProbe_Count\tProbe_Type\tExons_Skipped\t";
foreach my $sample (sort {$samples{$a}->{order} <=> $samples{$b}->{order}} keys %samples){
  print PROBESET "$sample","_LOG2_GNI\t";
}
foreach my $pair (sort {$pairs{$a}->{order} <=> $pairs{$b}->{order}} keys %pairs){
  print PROBESET "$pair","_LOG2_SI\t";
}
foreach my $rep (sort {$replicates{$a}->{order} <=> $replicates{$b}->{order}} keys %replicates){
  print PROBESET "$rep","_LOG2_GNI_U\t";
}
print PROBESET "mip_v_fur_SI_U\n";

#print Dumper $data_ref;

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
      my $raw_values_ref = $probes_ref->{$probe_id}->{raw_values};

      foreach my $sample (sort keys %{$raw_values_ref}){
	my $raw_value = $raw_values_ref->{$sample};

	#Gather probeset level values
	push (@{$probeset_values_ref->{$sample}->{raw_values}}, $raw_value);

	#Gather gene level values
	#If the probe type is correct (exon probe or canonical exon junction probe) add to the list of probe used to estimate gene expression
	#Also group all the negative control probes and come up with a single estimate for these (for interests sake)
	if ($probe_type eq "Exon" || (($probe_type eq "Exon-Exon") && ($exons_skipped eq "0")) || $probe_type eq "Control-Negative"){
	  push(@{$single_gene_ref->{gene_values}->{$sample}->{raw_values}}, $raw_value);
	}
      }

      #a.) Calculate the mean raw intensity across replicates ($rep = mip OR $rep = fur)
      my %mean_values;
      foreach my $rep (sort keys %replicates){
	my @reps = @{$replicates{$rep}{replicate}};

	my @replicate_values;
	foreach my $rep (@reps){
	  push (@replicate_values, $raw_values_ref->{$rep});
	}

	my $stat = Statistics::Descriptive::Full->new();
	$stat->add_data(@replicate_values);
	my $mean = $stat->mean();
	my $rounded = sprintf("%.4f", $mean);
	$mean_values{$rep} = $rounded;
      }
      $probes_ref->{$probe_id}->{mean_values} = \%mean_values;

    }#probe loop

    #Now that values have been grouped to the probeset level - perform calculations on these

    #b.) raw intensity (mean for entire probeset)
    foreach my $sample (sort keys %{$probeset_values_ref}){
      my @raw_values = @{$probeset_values_ref->{$sample}->{raw_values}};
      my $stat = Statistics::Descriptive::Full->new();
      $stat->add_data(@raw_values);
      my $mean = $stat->mean();
      my $rounded = sprintf("%.4f", $mean);
      $probeset_values_ref->{$sample}->{raw_mean} = $rounded;
    }

    #c.) grand mean of raw intensities (mean across replicates of the probeset means)
    my %mean_values;
    foreach my $rep (sort keys %replicates){
      my @samples = @{$replicates{$rep}{replicate}};

      my @raw_values;
      foreach my $sample (@samples){
	push (@raw_values, $probeset_values_ref->{$sample}->{raw_values});
      }
      my $stat = Statistics::Descriptive::Full->new();
      $stat->add_data(@raw_values);
      my $mean = $stat->mean();
      my $rounded = sprintf("%.4f", $mean);
      $mean_values{$rep} = $rounded;
    }
    $probesets_ref->{$probeset_id}->{mean_values} = \%mean_values;

  }#probesets loop

  #Now that values have been grouped to the gene level - perform calculations on these
  #d.) gene mean of raw values (mean across all exon and canonical probes for a single gene in a single experiment)
  my $gene_values_ref = $single_gene_ref->{gene_values};
  foreach my $sample (sort keys %{$gene_values_ref}){
    my @gene_raw_values = @{$gene_values_ref->{$sample}->{raw_values}};
    my $stat = Statistics::Descriptive::Full->new();
    $stat->add_data(@gene_raw_values);
    my $mean = $stat->mean();
    my $rounded = sprintf("%.4f", $mean);
    $gene_values_ref->{$sample}->{raw_mean} = $rounded;
  }

  #e.) grand gene mean of raw values (mean across all exon and canonical probes for a single gene across all experiments)
  my %mean_gene_values;
  foreach my $rep (sort keys %replicates){
    my @samples = @{$replicates{$rep}{replicate}};

    my @values;
    foreach my $sample (@samples){
      push (@values, $gene_values_ref->{$sample}->{raw_mean});
    }
    my $stat = Statistics::Descriptive::Full->new();
    $stat->add_data(@values);
    my $mean = $stat->mean();
    my $rounded = sprintf("%.4f", $mean);
    $mean_gene_values{$rep} = $rounded;
  }
  $single_gene_ref->{mean_values} = \%mean_gene_values;


  ###################################################################################################################################
  #Now that gene-level values are present use these to gene-normalize probe and probeset-level values
  foreach my $probeset_id (sort {$a <=> $b} keys %{$probesets_ref}){

    my $probeset_values_ref = $probesets_ref->{$probeset_id}->{probeset_values};
    my $probeset_pair_values_ref = $probesets_ref->{$probeset_id}->{probeset_pair_values};

    my $probe_type = $probesets_ref->{$probeset_id}->{probe_type};
    my $exons_skipped = $probesets_ref->{$probeset_id}->{exons_skipped};

    #Each PROBE
    my $probes_ref = $probesets_ref->{$probeset_id}->{probes};
    foreach my $probe_id (sort {$a <=> $b} keys %{$probes_ref}){

      #Each probe/sample value
      my $values_ref = $probes_ref->{$probe_id}->{raw_values};

      #i) Calculate gene normalized values for individual probes:  log2(exon1_mip/gene1_mip)
      my %log2_gni_values;
      foreach my $sample (sort keys %{$values_ref}){
	
	#Get the gene expression estimate for this gene for this sample
	
	my $gene_value = $gene_values_ref->{$sample}->{raw_mean};
	if ($gene_value == 0){
	  my $log2_gni_ratio = 'na';
	  $log2_gni_values{$sample} = 'na';
	}else{
	  my $log2_gni_ratio = (($values_ref->{$sample}) / $gene_value);

	  my $x = logn($log2_gni_ratio, 2);
	  my $rounded = sprintf("%.4f", $x);
	  $log2_gni_values{$sample} = $rounded;

	  #GATHER log2 gni values for each probe of this probeset for each sample
	  push(@{$probeset_values_ref->{$sample}->{log2_gni_values}}, $rounded);
	}

      }
      $probes_ref->{$probe_id}->{log2_gni_values} = \%log2_gni_values;

      #ii.) Calculate the mean log2 GNI
      my %mean_log2_gni_values;
      foreach my $rep (sort keys %replicates){
	my @samples = @{$replicates{$rep}{replicate}};

	my @values;
	foreach my $sample (@samples){
	  my $x = $log2_gni_values{$sample};
	  push (@values, $x);
	}
	my $stat = Statistics::Descriptive::Full->new();
	$stat->add_data(@values);
	my $mean = $stat->mean();
	my $rounded = sprintf("%.4f", $mean);
	$mean_log2_gni_values{$rep} = $rounded;
      }
      $probes_ref->{$probe_id}->{mean_log2_gni_values} = \%mean_log2_gni_values;


      #iii.) Calculate the SI value for individual probes:  log2((exon1_mip/gene1_mip) / (exon1_fur/gene1_fur))

      my %log2_si_values;
      foreach my $pair (sort keys %pairs){
	my @pair = @{$pairs{$pair}{pair}};

	my $sample1_gene_value = $gene_values_ref->{$pair[0]}->{raw_mean};
	my $sample2_gene_value = $gene_values_ref->{$pair[1]}->{raw_mean};

	if ($sample1_gene_value == 0 || $sample2_gene_value == 0){
	  $log2_si_values{$pair} = 'na';
	}else{

	  my $si_ratio = (($values_ref->{$pair[0]} / $sample1_gene_value) / ($values_ref->{$pair[1]} / $sample2_gene_value));
	  my $x = logn($si_ratio, 2);
	  my $rounded = sprintf("%.4f", $x);
	  $log2_si_values{$pair} = $rounded;

	  #GATHER log2 si values for each probe sample pair of this probeset
	  push(@{$probeset_pair_values_ref->{$pair}->{log2_si_values}}, $rounded);
	}
      }
      $probes_ref->{$probe_id}->{log2_si_values} = \%log2_si_values;

      #iv.) Calculate the mean log2 SI
      my @si_values;
      foreach my $pair (keys %log2_si_values){
	push(@si_values, $log2_si_values{$pair});
      }
      my $stat = Statistics::Descriptive::Full->new();
      $stat->add_data(@si_values);
      my $mean = $stat->mean();
      my $rounded = sprintf("%.4f", $mean);

      $probes_ref->{$probe_id}->{mean_log2_si_value} = $rounded;

    }#probe loop

    #Each PROBESET
    #Now that values have been grouped to the probeset level - perform calculations on these

    #v.) log2 gni (mean for entire probeset)
    foreach my $sample (sort keys %{$probeset_values_ref}){
      my @log2_gni_values = @{$probeset_values_ref->{$sample}->{log2_gni_values}};
      my $stat = Statistics::Descriptive::Full->new();
      $stat->add_data(@log2_gni_values);
      my $mean = $stat->mean();
      my $rounded = sprintf("%.4f", $mean);

      $probeset_values_ref->{$sample}->{log2_gni_mean} = $rounded;
    }

    #vi.) log2 SI (mean for entire probeset)
    foreach my $pair (sort keys %{$probeset_pair_values_ref}){
      my @log2_si_values = @{$probeset_pair_values_ref->{$pair}->{log2_si_values}};
      my $stat = Statistics::Descriptive::Full->new();
      $stat->add_data(@log2_si_values);
      my $mean = $stat->mean();
      my $rounded = sprintf("%.4f", $mean);

      $probeset_pair_values_ref->{$pair}->{log2_si_mean} = $rounded;
    }

    #vii.) grand mean of log2 gni (mean across replicates of the probeset means)
    my %mean_log2_gni_values;
    foreach my $rep (sort keys %replicates){
      my @samples = @{$replicates{$rep}{replicate}};

      my @values;
      foreach my $sample (@samples){
	push (@values, @{$probeset_values_ref->{$sample}->{log2_gni_values}});
      }
      my $stat = Statistics::Descriptive::Full->new();
      $stat->add_data(@values);
      my $mean = $stat->mean();
      my $rounded = sprintf("%.4f", $mean);
      $mean_log2_gni_values{$rep} = $rounded;
    }
    $probesets_ref->{$probeset_id}->{log2_gni_mean} = \%mean_log2_gni_values;

    #viii.) grand mean of log2 de (mean across replicates of probeset means)
    my @probeset_si_values;
    foreach my $pair (sort keys %{$probeset_pair_values_ref}){
      push (@probeset_si_values, @{$probeset_pair_values_ref->{$pair}->{log2_si_values}});
    }
    my $stat = Statistics::Descriptive::Full->new();
    $stat->add_data(@probeset_si_values);
    my $mean = $stat->mean();
    my $rounded = sprintf("%.4f", $mean);
    $probesets_ref->{$probeset_id}->{log2_si_mean} = $rounded;

  }#probesets loop



  ####################################################################################################################################
  #PRINT OUTPUT FILES
  #Calculations complete for this gene.  Now print results to output files

  foreach my $probeset_id (sort {$a <=> $b} keys %{$probesets_ref}){

    my $probe_type = $probesets_ref->{$probeset_id}->{probe_type};
    my $exons_skipped = $probesets_ref->{$probeset_id}->{exons_skipped};

    #NOTE: 'gni' stands for gene-normalized intensity

    #Each probe
    my $probes_ref = $probesets_ref->{$probeset_id}->{probes};
    foreach my $probe_id (sort {$a <=> $b} keys %{$probes_ref}){

      #Probe-level
      #Probe_ID, ProbeSet_ID, AlexaGene_ID, Probe_Type, Exons_Skipped, log2_gni (all reps), log2_si (all_pairs), mean_log2_gni (mean of reps), mean_log2_si (mean of pairs)
      print PROBE "$probe_id\t$probeset_id\t$gene_id\t$probe_type\t$exons_skipped\t";

      my $log2_gni_values_ref = $probes_ref->{$probe_id}->{log2_gni_values};
      foreach my $sample (sort {$samples{$a}->{order} <=> $samples{$b}->{order}} keys %samples){
	print PROBE "$log2_gni_values_ref->{$sample}\t";
      }
      my $log2_si_values_ref = $probes_ref->{$probe_id}->{log2_si_values};
      foreach my $pair (sort {$pairs{$a}->{order} <=> $pairs{$b}->{order}} keys %pairs){
	print PROBE "$log2_si_values_ref->{$pair}\t";
      }

      my $mean_log2_gni_values_ref = $probes_ref->{$probe_id}->{mean_log2_gni_values};
      foreach my $rep (sort {$replicates{$a}->{order} <=> $replicates{$b}->{order}} keys %replicates){
	print PROBE "$mean_log2_gni_values_ref->{$rep}\t";
      }
      print PROBE "$probes_ref->{$probe_id}->{mean_log2_si_value}\n";
    }
  }

  #Probeset-level
  #ProbeSet_ID, AlexaGene_ID, Probe_Count, Probe_Type, Exons_Skipped, log2_gni (all reps), log2_si (all_pairs), mean_log2_gni (mean of reps), mean_log2_si (mean of pairs)
  foreach my $probeset_id (sort {$a <=> $b} keys %{$probesets_ref}){

    my $probe_type = $probesets_ref->{$probeset_id}->{probe_type};
    my $probe_count = $probesets_ref->{$probeset_id}->{probe_count};
    my $exons_skipped = $probesets_ref->{$probeset_id}->{exons_skipped};

    my $probeset_values_ref = $probesets_ref->{$probeset_id}->{probeset_values};
    my $probeset_pair_values_ref = $probesets_ref->{$probeset_id}->{probeset_pair_values};

    print PROBESET "$probeset_id\t$gene_id\t$probe_count\t$probe_type\t$exons_skipped\t";

    foreach my $sample (sort {$samples{$a}->{order} <=> $samples{$b}->{order}} keys %samples){
      print PROBESET "$probeset_values_ref->{$sample}->{log2_gni_mean}\t";
    }
    foreach my $pair (sort {$pairs{$a}->{order} <=> $pairs{$b}->{order}} keys %pairs){
      print PROBESET "$probeset_pair_values_ref->{$pair}->{log2_si_mean}\t";
    }

    my $mean_log2_gni_values_ref = $probesets_ref->{$probeset_id}->{log2_gni_mean};
    foreach my $rep (sort {$replicates{$a}->{order} <=> $replicates{$b}->{order}} keys %replicates){
      print PROBESET "$mean_log2_gni_values_ref->{$rep}\t";
    }
    print PROBESET "$probesets_ref->{$probeset_id}->{log2_si_mean}\n";
  }

  #print Dumper $single_gene_ref;

}#gene loop

print BLUE "\n\n", RESET;

close (PROBE);
close (PROBESET);
close (GENE);

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

      if ($probe_type eq "Exon" || (($probe_type eq "Exon-Exon") && ($exons_skipped eq "0")) || $probe_type eq "Control-Negative"){
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
	    my @tmp2;
	    $probeset_values{$sample}{raw_values} = \@tmp1;
	    $probeset_values{$sample}{log2_gni_values} = \@tmp2;
	  }
	}
	$probes{$probe_id}{raw_values} = \%raw_values;

	$probesets_ref->{$probeset_id}->{probeset_values} = \%probeset_values;
	$probesets_ref->{$probeset_id}->{probes} = \%probes;
	$probesets_ref->{$probeset_id}->{probe_type} = $probe_type;
	$probesets_ref->{$probeset_id}->{exons_skipped} = $exons_skipped;
	$probesets_ref->{$probeset_id}->{probe_count} = 1;

	#Initialize arrays to store de values for all samples pairs for the (($values_ref->{$pair[0]} / $sample1_gene_value) / ($values_ref->{$pair[1]} / $sample2_gene_value)entire probeset
	my %probeset_pair_values;
	foreach my $pair (sort keys %pairs){
	  my @tmp1;
	  $probeset_pair_values{$pair}{log2_si_values} = \@tmp1;
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
	  my @tmp3;
	  $probeset_values{$sample}{raw_values} = \@tmp1;
	  $probeset_values{$sample}{log2_gni_values} = \@tmp2;
	  $gene_values{$sample}{raw_values} = \@tmp3;

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
	$probeset_pair_values{$pair}{log2_si_values} = \@tmp1;

      }
      $probesets{$probeset_id}{probeset_pair_values} = \%probeset_pair_values;

      $data{$gene_id}{probesets} = \%probesets;
      $data{$gene_id}{gene_values} = \%gene_values;
      $data{$gene_id}{gene_pair_values} = \%gene_pair_values;

      if ($probe_type eq "Exon" || (($probe_type eq "Exon-Exon") && ($exons_skipped eq "0")) || $probe_type eq "Control-Negative"){
	$data{$gene_id}{gene_probe_count} = 1;  #Number of probes that will be used to estimate gene expression
      }
    }
  }

  close (DATA);

  return(\%data);
}












