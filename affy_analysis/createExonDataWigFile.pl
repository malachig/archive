#!/usr/bin/perl -w
#Written by Malachi Griffith

#Create a wiggle custom ucsc file for all HL60 Affymetrix exon array data

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);
use Math::Complex;


#Initialize command line options
my $probeset_data_file = '';  #Contains affymetrix exon array data (expression and DE values for all exon probesets - average across replicates)
                              # - this file should also contain coordinate data for all probesets
my $param_file = '';          #Contains description of each data column to be made into a wig track (name, priority, etc.)
my $output_dir = '';          #Directory where wig files will be written, one per chromosome


my $probeset_data_file = "/home/malachig/AlternativeSplicing/Array_analysis/Affymetrix_Exon_Arrays/HL60/ExpressionConsoleResults/iterplier-exon-all/AllExons_iterplier_DAGB_all_means_withCoords_hg18_positionSorted.txt";

my $output_dir = "/home/malachig/AlternativeSplicing/Array_analysis/Affymetrix_Exon_Arrays/HL60/ExpressionConsoleResults/iterplier-exon-all/wigfiles/";

my %data;

&parseProbeDataFile('-probeset_data_file'=>$probeset_data_file);

my $duplicate_probesets = 0;
my $overlap_probesets = 0;

my $plier_file = "$output_dir"."AffyExonData_BothStrands_Plier_hg18.txt";
open (PLIER, ">$plier_file") || die "\nCould not open PLIER file: $plier_file\n\n";

my $dabg_file = "$output_dir"."AffyExonData_BothStrands_DABG_hg18.txt";
open (DABG, ">$dabg_file") || die "\nCould not open DABG file: $dabg_file\n\n";


foreach my $strand (sort keys %data){

  #One file for each strand
  if ($strand eq '+'){

    print PLIER "\ntrack type=wiggle_0 name=\"AffyExon_Plier_P\" description=\"Affymetrix exon array PLIER expression data - Positive Strand\" visibility=full color=0,100,0 autoScale=off viewLimits=0.0:14.6 yLineMark=8.88 yLineOnOff=on priority=97";

    print DABG "\ntrack type=wiggle_0 name=\"AffyExon_DABG_P\" description=\"Affymetrix exon array DABG P-Value - Log10 Scale - Positive Strand\" visibility=full color=178,34,34 autoScale=off viewLimits=0.0:3.0 yLineMark=1.30103 yLineOnOff=on priority=99";


  }elsif ($strand eq '-'){

    print PLIER "\ntrack type=wiggle_0 name=\"AffyExon_Plier_N\" description=\"Affymetrix exon array PLIER expression data - Negative Strand\" visibility=full color=0,134,139 altColor=0,134,139 autoScale=off viewLimits=0.0:-14.6 yLineMark=-8.88 yLineOnOff=on priority=98";

    print DABG "\ntrack type=wiggle_0 name=\"AffyExon_DABG_N\" description=\"Affymetrix exon array DABG P-Value - Log10 Scale - Negative Strand\" visibility=full color=255,69,0 altColor=255,69,0 autoScale=off viewLimits=0.0:-3.0 yLineMark=-1.30103 yLineOnOff=on priority=100";

  }else{
    print RED, "\nStrand not understood\n\n", RESET;
    exit();
  }

  my $chr_ref = $data{$strand}{chromosomes};

  foreach my $chr (sort {$a cmp $b} keys %{$chr_ref}){
    my $probesets_ref = $chr_ref->{$chr}->{probesets};

    my $current_start = 0;
    my $current_end = 0;

    foreach my $probeset_id (sort {$probesets_ref->{$a}->{start} <=> $probesets_ref->{$b}->{start}} keys %{$probesets_ref}){
      my $start = $probesets_ref->{$probeset_id}->{start};
      my $end = $probesets_ref->{$probeset_id}->{end};

      #Watch for duplicate probesets
      #This is a hack which should be dealt with properly (average of duplicates?)

      if (($start == $current_start) && ($end == $current_end)){
	$duplicate_probesets++;
	next();
      }

      #Watch for overlaping probesets
      #This is a hack which should be dealt with properly (merge overlaping probesets?)

      if ($start < $current_end){
	$overlap_probesets++;
	next();
      }

      if ($strand eq '+'){
	print PLIER "\nchr$chr\t$start\t$end\t$probesets_ref->{$probeset_id}->{mean_plier_log2}";
	print DABG "\nchr$chr\t$start\t$end\t$probesets_ref->{$probeset_id}->{dabg_pvalue_log10}";

      }elsif ($strand eq '-'){
	print PLIER "\nchr$chr\t$start\t$end\t-$probesets_ref->{$probeset_id}->{mean_plier_log2}";
	print DABG "\nchr$chr\t$start\t$end\t-$probesets_ref->{$probeset_id}->{dabg_pvalue_log10}";
      }

      $current_start = $start;
      $current_end = $end;

    }
  }
}

close (PLIER);
close (DABG);

print YELLOW, "\nFound $duplicate_probesets duplicate probesets!  Only using the first occurence when this happens\n\n", RESET;
print YELLOW, "\nFound $overlap_probesets overlap probesets!  Only using the first occurence when this happens\n\n", RESET;

my $cmd = "/usr/bin/gzip $output_dir*.txt";
system ($cmd);
print BLUE, "\n$cmd\n\n", RESET;

exit();


#############################################################################################################################
#Parse probe input data                                                                                                     #
#############################################################################################################################
sub parseProbeDataFile{
  my %args = @_;
  my $probeset_data_file = $args{'-probeset_data_file'};

  print GREEN, "\n\nParsing Affy probeset expression file: $probeset_data_file\n", RESET;

  open (DATA, "$probeset_data_file") || die "\n\nCould not open probe data file: $probeset_data_file\n\n";

  my $probe_record_count = 0;      #Number of probeset records examined

  my $count = 0;


  while(<DATA>){
    chomp($_);
    my $line = $_;

    $count++;

    if ($count == 10000){
      $| = 1;
      print GREEN, ".", RESET;
      $| = 0;
      $count = 0;
    }

    #Remove the quote characters from each line
    my @line = split ("\t", $line);

    #Unless the first element of this line array contains a number, skip the line (header lines, etc.)
    unless ($line[0] =~ /^\d+/){
      next();
    }

    #Grab the following data: probeset_id,transcript_cluster_id,seqname,strand,start,stop,probe_count
    my $probeset_id = $line[0];
    my $chr = $line[1];
    my $start = $line[2];
    my $end = $line[3];
    my $strand = $line[4];
    my $mean_plier_log2 = $line[5];
    my $dabg = $line[6];

    #Format DABG P-value and convert to log10 scale
    if ($dabg == 0){
      $dabg = 0.001;
    }

    my $dabg_log10 = abs(logn($dabg, 10));

    #Strand -> Chr -> ProbeSet
    if ($data{$strand}){
      my $chr_ref = $data{$strand}{chromosomes};

      if ($chr_ref->{$chr}){
	my $probesets_ref = $chr_ref->{$chr}->{probesets};
	$probesets_ref->{$probeset_id}->{start} = $start;
	$probesets_ref->{$probeset_id}->{end} = $end;
	$probesets_ref->{$probeset_id}->{mean_plier_log2} = $mean_plier_log2;
	$probesets_ref->{$probeset_id}->{dabg_pvalue_log10} = $dabg_log10;
	$probe_record_count++;

      }else{
	#First time this chromosome has been observed
	my %probesets;
	$probesets{$probeset_id}{start} = $start;
	$probesets{$probeset_id}{end} = $end;
	$probesets{$probeset_id}{mean_plier_log2} = $mean_plier_log2;
	$probesets{$probeset_id}{dabg_pvalue_log10} = $dabg_log10;

	$probe_record_count++;

	$chr_ref->{$chr}->{probesets} = \%probesets;
      }

    }else{
      #First this time strand has been observed
      my %chromosomes;
      my %probesets;
      $probesets{$probeset_id}{start} = $start;
      $probesets{$probeset_id}{end} = $end;
      $probesets{$probeset_id}{mean_plier_log2} = $mean_plier_log2;
      $probesets{$probeset_id}{dabg_pvalue_log10} = $dabg_log10;

      $probe_record_count++;

      $chromosomes{$chr}{probesets} = \%probesets;
      $data{$strand}{chromosomes} = \%chromosomes;
    }

    #DEBUG
    #if ($probe_record_count == 100000){
    #  last();
    #}

  }

  #Print summary statistics for the mapping of probesets
  print BLUE, "\n\nTotal probeset records processed: $probe_record_count\n\n", RESET;
  close DATA;

  return();
}

