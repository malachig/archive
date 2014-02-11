#!/usr/bin/perl -w
#Written by Malachi Griffith

#This file takes as input a mapping file containing probeset ID to Ensembl Gene mappings and creates a
#meta-probeset file for use with the Affymetrix ExACT or Expression Console software package.  Essentially I am defining meta-probesets for
#my own gene annotations (ALEXA/Ensembl in this case)
#This will allow me to generate gene-level PLIER estimates of expression for each Ensembl Gene

use strict;
use Data::Dumper;
use Getopt::Long;
use Term::ANSIColor qw(:constants);

#Initialize command line options
my $infile = '';
my $outfile = '';

GetOptions ('infile=s'=>\$infile,'outfile=s'=>\$outfile);

#Provide instruction to the user
print BLUE, "\n\nUsage:", RESET;
print BLUE, "\n\tSpecify the input probeset-to-gene mapfile using: --infile", RESET;
print BLUE, "\n\t\tThis file should have already been created by the scripts mapAffyProbesToEnsembl.pl and identifyCoreEnsemblProbesets.pl", RESET;
print BLUE, "\n\n\tSpecify the output file name for the Expression Console metaprobeset file using: --outfile", RESET;
print YELLOW, "\n\n\tNOTE: Genome version info for header of output file is hard coded!! Make sure it is correct!!", RESET;
print YELLOW, "\n\t\tYou can get this info from one of the .mps files provided by Affymetrix", RESET;
print BLUE, "\n\nExample: createEnsemblMetaProbesetFile.pl --infile probesets_mapped_to_ensembl.txt --outfile meta-probeset.ensembl_hs_47_36i.mps\n\n", RESET;

#Make sure all options were specified
unless ($infile && $outfile){
  print RED, "\nOptions missing!\n\n", RESET;
  exit();
}

#1.) Read input file with mapped probset-gene relationships
my %gene_probesets;
&getProbesetGeneData('-mapfile'=>$infile);

#2.) Dump output file in the required Exact format
&writeMetaProbesetFile('-outfile'=>$outfile);

exit();


#####################################################################################################################
#1.) Read input file with mapped probset-gene relationships                                                         #
#####################################################################################################################
sub getProbesetGeneData{
  my %args = @_;
  my $mapfile = $args{'-mapfile'};

  open (MAP, "$mapfile") || die "\nCould not open input mapfile: $mapfile\n\n";

  #Info needed: probeset_id     transcript_cluster_id   probeset_list   probe_count
  #Use:         probeset_id     alexa_id                list probesets  cumulative_probe_count

  my $probe_count;
  while (<MAP>){

    unless ($_ =~ /^\d+/){
      next();
    }
    $probe_count++;

    chomp($_);
    my @line = split ("\t",$_);

    my $probe_set_id = $line[0];
    my $alexa_gene_id = $line[1];
    my $probe_count = $line[8];

    #If this gene has already been observed add this probeset to its record, otherwise create a new gene record
    if ($gene_probesets{$alexa_gene_id}){
      my @probeset_list = @{$gene_probesets{$alexa_gene_id}{probeset_list}};
      push(@probeset_list, $probe_set_id);
      $gene_probesets{$alexa_gene_id}{probeset_list} = \@probeset_list;
      $gene_probesets{$alexa_gene_id}{probe_count} += $probe_count;
    }else{
      my @probeset_list;
      push(@probeset_list, $probe_set_id);
      $gene_probesets{$alexa_gene_id}{probeset_list} = \@probeset_list;
      $gene_probesets{$alexa_gene_id}{probe_count} = $probe_count;
    }
  }
  close (MAP);

  return();
}


#####################################################################################################################
#2.) Dump output file in the required Exact format                                                                  #
#####################################################################################################################
sub writeMetaProbesetFile{
  my %args = @_;
  my $outfile = $args{'-outfile'};

  open (OUT, ">$outfile") || die "\nCould not open outfile: $outfile\n\n";

  #print file header info
  my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
  my @weekDays = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
  my ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
  my $year = 1900 + $yearOffset;
  my $theTime = "$hour:$minute:$second, $weekDays[$dayOfWeek] $months[$month] $dayOfMonth, $year";

  print OUT "\#\%chip_type\=HuEx-1_0-st-v2\n";
  print OUT "\#\%chip_type\=HuEx-1_0-st-v1\n";
  print OUT "\#\%chip_type\=HuEx-1_0-st-ta1\n";
  print OUT "\#\%lib-set-name\=HuEx-1_0-st\n";
  print OUT "\#\%lib-set-version\=v2\n";
  print OUT "\#\%create_date\=$theTime\n";
  print OUT "\#\%genome-species\=Homo sapiens\n";
  print OUT "\#\%genome-version\=hg18\n";
  print OUT "\#\%genome-version-ucsc\=hg18\n";
  print OUT "\#\%genome-version-ncbi\=36\n";
  print OUT "\#\%genome-version-create_date\=2006 March\n";
  print OUT "probeset_id\ttranscript_cluster_id\tprobeset_list\tprobe_count\n";

  foreach my $gene_id (sort {$a <=> $b} keys %gene_probesets){

    my @probeset_list = @{$gene_probesets{$gene_id}{probeset_list}};

    print OUT "$gene_id\t$gene_id\t@probeset_list\t$gene_probesets{$gene_id}{probe_count}\n";

  }

  return();
}







